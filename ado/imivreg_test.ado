// imivreg.ado
program define imivreg, eclass
    version 17
    syntax varlist(numeric min=1) [if] [in] [aw fw iw pw], ///
        CLuster(varlist numeric) [TIMEperiod(varname numeric)] ///
        IV(string asis)
    
    if !regexm("`iv'", "^([^=]+)=(.+)$") {
        di as error "iv() must be specified as endogenous_vars = instrument_vars"
        exit 198
    }
    
    local endog_var = trim(regexs(1))
    local instrument_vars = trim(regexs(2))
    di as text "Endogenous vars: `endog_var'"
    di as text "Instruments: `instrument_vars'"
	
    
    // varlist contains depvar and exog_vars
    tokenize `varlist'
    local depvar "`1'"
    macro shift
    local exog_vars `*'

    // Check if endogenous variable is in exog_vars (it shouldn't be)
    local present : list endog_var in exog_vars
    if `present' {
        di as err "Endogenous variable (`endog_var') cannot also be listed as an exogenous variable."
        exit 198
    }

    tempname __htscalar
    if "`timeperiod'" == "" {
        tempvar __timedummy
        quietly gen double `__timedummy' = 0 `if' `in'
        local timeperiod `__timedummy'
        scalar `__htscalar' = 0
    }
    else {
        quietly summarize `timeperiod' `if' `in'
        scalar `__htscalar' = (r(min) < r(max))
    }

    // Create a list of all structural regressors (endog + exog) for output naming
    local all_struc_regressors `endog_var' `exog_vars'

    mata{
        uniformseed(12345)
        rseed(12345)

        Y_s         = st_data(., "`depvar'")      
        X_e         = st_data(., "`endog_var'")
        if ("`exog_vars'" != "") {
            X_k     = st_data(., tokens("`exog_vars'"))
        } else {
            X_k     = J(rows(Y_s), 0, .) // Empty matrix if no exog_vars
        }
        Z_i         = st_data(., tokens("`instrument_vars'"))
        
        coord       = st_data(.,"`cluster'")
        timePeriod  = st_data(.,"`timeperiod'")
        n           = length(Y_s)

        // Structural regressors: X_s_all = (X_e, X_k, constant)
        X_s_all     = (X_e, X_k, J(n,1,1))
        p_s         = cols(X_s_all) // Number of structural parameters including constant

        // Instruments for structural equation: Z_s_all = (Z_i, X_k, constant)
        Z_s_all     = (Z_i, X_k, J(n,1,1))

        // Regressors for first stage of X_e: X_fs_all = (Z_i, X_k, constant)
        X_fs_all    = (Z_i, X_k, J(n,1,1))
        p_fs        = cols(X_fs_all)

        G_max = ceil(n^(1/3))
        st_numscalar("G_max", G_max)
        G_vec = range(2, G_max, 1)
        l_G = length(G_vec)
        st_numscalar("l_G", l_G)
    }

    gen id0 = _n 
    matrix dissim dis_mat = `cluster',L2 
    matrix dissim time_mat = `timeperiod', L2
    forvalues i =2/`=G_max' { 
        qui clpam group`i', distmat(dis_mat) id(id0) medoids(`i') ga
    }

    mata{
        hasTime = st_numscalar("`__htscalar'")
        data_medoids = st_data(.,.)
        ncols_med = cols(data_medoids)
        clusteringSet =  data_medoids[.,ncols_med-G_max+2..ncols_med]

        // Initial Full Sample Estimations for residuals U_hat and V_hat
        // Structural equation: Y_s on X_s_all with instruments Z_s_all
        beta_s_2sls = invsym(X_s_all'*Z_s_all*invsym(Z_s_all'*Z_s_all)*Z_s_all'*X_s_all) * (X_s_all'*Z_s_all*invsym(Z_s_all'*Z_s_all)*Z_s_all'*Y_s)
        U_hat = Y_s - X_s_all * beta_s_2sls
        
        // First stage: X_e on X_fs_all (OLS)
        pi_fs_ols = invsym(X_fs_all'*X_fs_all)*X_fs_all'*X_e
        V_hat = X_e - X_fs_all * pi_fs_ols
        
        M_ident = I(n)
        Qd_ident=Rd_ident=ex_ident=.
        qrdp(M_ident, Qd_ident, Rd_ident, ex_ident) // ex_ident will be 1..n
        useQML_all = ex_ident[1..n] // Use all observations for QMLE on these residuals

        dis_mat = st_matrix("dis_mat")
        time_mat = st_matrix("time_mat")

        // QMLE for Sigma_U using U_hat
        sig2_0_U = log(mean(U_hat:^2))
        d_med_U  = mm_median(select(colshape(dis_mat,1), colshape(dis_mat:>0,1)))
        rho_0_U  = log(d_med_U/(-ln(0.30)))
        if (hasTime) {
            t_vec_U = select( colshape(time_mat,1) , colshape(time_mat :> 0 , 1) )
            t_med_U  = (rows(t_vec_U)>0 ? mm_median(t_vec_U) : 1)
            tau_0_U  = log(t_med_U/(-ln(0.30)))
        }
        S_U = optimize_init()
        if (hasTime) {
            optimize_init_evaluator(S_U,  &QMLE_new() )
            optimize_init_params(S_U, (sig2_0_U, rho_0_U, tau_0_U))
            optimize_init_argument(S_U, 5, U_hat) // Pass U_hat as resid
        } else {
            optimize_init_evaluator(S_U,  &QMLE_bin() )
            optimize_init_params(S_U, (sig2_0_U, rho_0_U))
            optimize_init_argument(S_U, 4, U_hat) // Pass U_hat as resid
        }
        optimize_init_which(S_U, "min")
        optimize_init_evaluatortype(S_U, "d0")
        optimize_init_argument(S_U, 1, M_ident) // Effectively no projection for QMLE itself
        optimize_init_argument(S_U, 2, useQML_all)
        optimize_init_argument(S_U, 3, dis_mat)
        if (hasTime) optimize_init_argument(S_U, 4, time_mat)
        
        alphaHat_U = optimize(S_U)
        if (hasTime) {
			SigmaHat_U = Sigma_func_DGP(alphaHat_U, dis_mat, time_mat)
			}
        else { 
            SigmaHat_U = Sigma_func_DGP_bin(alphaHat_U, dis_mat)
            SigmaHat_U = SigmaHat_U + I(rows(SigmaHat_U))*mean(U_hat:^2)*1e-2 // Regularization from ccereg
        }

        // QMLE for Sigma_V using V_hat
        sig2_0_V = log(mean(V_hat:^2))
        d_med_V  = d_med_U // Assume same distance metric applies
        rho_0_V  = rho_0_U
        if (hasTime) {
            t_med_V = t_med_U
            tau_0_V = tau_0_U
        }
        S_V = optimize_init()
        if (hasTime) {
            optimize_init_evaluator(S_V,  &QMLE_new() )
            optimize_init_params(S_V, (sig2_0_V, rho_0_V, tau_0_V))
            optimize_init_argument(S_V, 5, V_hat) // Pass V_hat as resid
        } else {
            optimize_init_evaluator(S_V,  &QMLE_bin() )
            optimize_init_params(S_V, (sig2_0_V, rho_0_V))
            optimize_init_argument(S_V, 4, V_hat) // Pass V_hat as resid
        }
        optimize_init_which(S_V, "min")
        optimize_init_evaluatortype(S_V, "d0")
        optimize_init_argument(S_V, 1, M_ident)
        optimize_init_argument(S_V, 2, useQML_all)
        optimize_init_argument(S_V, 3, dis_mat)
        if (hasTime) optimize_init_argument(S_V, 4, time_mat)

        alphaHat_V = optimize(S_V)
        if (hasTime) {
			SigmaHat_V = Sigma_func_DGP(alphaHat_V, dis_mat, time_mat) 
			}
        else { 
            SigmaHat_V = Sigma_func_DGP_bin(alphaHat_V, dis_mat) 
            SigmaHat_V = SigmaHat_V + I(rows(SigmaHat_V))*mean(V_hat:^2)*1e-2 // Regularization
        }
        
        // Estimate rhoHat
        CSHat_U = cholesky(SigmaHat_U)
        CSHat_V = cholesky(SigmaHat_V)
        U_transformed = invsym(CSHat_U) * U_hat // L^-1 * U
        V_transformed = invsym(CSHat_V) * V_hat // L^-1 * V
        rhoHat = correlation(U_transformed, V_transformed)
		display("AAA")
        if (missing(rhoHat)) rhoHat = 0 // Handle potential issues if correlation is undefined

        // Construct Sigma_sim_block for (U,V)
        Sigma_UV = rhoHat * CSHat_U * CSHat_V'
        Sigma_sim_block = J(2*n, 2*n, 0)
        Sigma_sim_block[1..n, 1..n] = SigmaHat_U
        Sigma_sim_block[(n+1)..(2*n), (n+1)..(2*n)] = SigmaHat_V
        Sigma_sim_block[1..n, (n+1)..(2*n)] = Sigma_UV
        Sigma_sim_block[(n+1)..(2*n), 1..n] = Sigma_UV'
        
        // Ensure Sigma_sim_block is positive definite for Cholesky
        // Small ridge if needed, though QMLE regularization might help
        eigval = symeigensystem(Sigma_sim_block, ., .)
        if (min(eigval) <= 1e-8) { // Check for positive definiteness
            Sigma_sim_block = Sigma_sim_block + I(2*n)*max((1e-6, -min(eigval)+1e-6))
        }

        sigLevel = .05
        Bboot = 1000 // Reduce for speed, paper uses 10000 for example, 1000 for sim
        CSHat_sim_block = cholesky(Sigma_sim_block)
        UVbootMat = CSHat_sim_block' * rnormal(2*n, Bboot, 0, 1)
        
        resultsMat = J(p_s, 7, 0) // p_s = number of structural parameters (X_e, X_k, const)
        GstarVec = J(p_s, 1, 0)

        // Loop for each structural parameter
        for(iCov = 1; iCov <= p_s; iCov++ ){
            // beta_s_H0: structural parameters under H0 for current iCov
            beta_s_H0 = beta_s_2sls // Full sample 2SLS estimates
            beta_s_H0[iCov] = 0     // Set coeff of interest to 0 for H0
            
            simPowerVec = J(l_G,1,0)
            pValSim = J(Bboot,l_G,0)
            sigLevelAdjVec = J(l_G,1,0)

            for(kk = 1; kk <= l_G; kk++){
                clustering = clusteringSet[.,kk]
                G = G_vec[kk]
                
                theta_C_vec_sim = J(Bboot, 1, 0) // Stores mean(theta_g) for each sim
                se_of_mean_theta_C_vec_sim = J(Bboot, 1, 0) // Stores SE(mean(theta_g))

                for(rr = 1; rr <= Bboot; rr++){
                    U_boot_iter = UVbootMat[1..n, rr]
                    V_boot_iter = UVbootMat[(n+1)..(2*n), rr]
                    
                    X_e_boot = X_fs_all * pi_fs_ols + V_boot_iter
                    X_s_all_boot = X_s_all // Start with original
                    X_s_all_boot[.,1] = X_e_boot // Replace endogenous var with simulated one (assuming X_e is first col of X_s_all)
                                                // This needs care if X_e is not first. Let's assume X_s_all = (X_e, X_k, const)
                                                // And iCov refers to index in this X_s_all.

                    Y_s_boot = X_s_all_boot * beta_s_H0 + U_boot_iter
                    
                    // Perform FamaMacbethIV for this simulated dataset
                    // It should return G x p_s matrix of coeffs
                    b_coeffs_clusterwise_sim = J(G, p_s, .)
                    FamaMacbethIV(Y_s_boot, X_s_all_boot, Z_s_all, clustering, b_coeffs_clusterwise_sim)
                    
                    current_theta_C_sim = b_coeffs_clusterwise_sim[.,iCov] // Coeff of interest for this simulation
                    
                    theta_C_vec_sim[rr] = mean(current_theta_C_sim)
                    if (G > 1) {
                        se_of_mean_theta_C_vec_sim[rr] = sqrt(variance(current_theta_C_sim)/G)
                    } else {
                        se_of_mean_theta_C_vec_sim[rr] = . // Undefined for G=1
                    }
                    
                    if (se_of_mean_theta_C_vec_sim[rr] != . & se_of_mean_theta_C_vec_sim[rr] > 1e-9) {
                         pValSim[rr,kk] = 2*t(G-1, -abs(theta_C_vec_sim[rr]/se_of_mean_theta_C_vec_sim[rr]))
                    } else {
                         pValSim[rr,kk] = 1
                    }
                }

                if (G > 1) {
                    sigLevelAdjVec[kk] = min((mm_quantile(pValSim[.,kk],1,0.05),0.05))
                    if (missing(sigLevelAdjVec[kk])) sigLevelAdjVec[kk] = 0.05 // Fallback
                } else {
                    sigLevelAdjVec[kk] = 0.05 // Cannot adjust for G=1
                }
                
                // Power calculation
                alternatives = (-10::1)' * (1::10)/sqrt(n)
                nalt = rows(alternatives)
                power_sum_for_G = 0
                if (G > 1) {
                    for (alt_idx = 1; alt_idx <= nalt; alt_idx++) {
                        current_alt_val = alternatives[alt_idx]
                        // theta_C_vec_sim are means of (theta_g - 0). For power, it's mean(theta_g - alt_val)
                        // Or rather, t_stat is (mean(theta_g) - alt_val) / SE(mean(theta_g))
                        t_stats_alt = (theta_C_vec_sim :- current_alt_val) :/ se_of_mean_theta_C_vec_sim
                        pvals_for_alt = 2*t(G-1, -abs(t_stats_alt))
                        power_sum_for_G = power_sum_for_G + mean(pvals_for_alt :< sigLevelAdjVec[kk])
                    }
                    simPowerVec[kk] = power_sum_for_G / nalt
                } else {
                    simPowerVec[kk] = 0 // No power for G=1
                }
            } // end kk loop (G_vec)

            indStar = windStar = .
            maxindex(simPowerVec,1,indStar,windStar)
            if (rows(indStar) > 1) indStar = indStar[1]
            if (missing(indStar)) indStar = l_G // Fallback if all power is 0 or missing
            
            Gstar = G_vec[indStar]
            clusteringStar = clusteringSet[.,indStar]
            GstarVec[iCov] = Gstar

            // Final estimation on original data with Gstar
            b_coeffs_clusterwise_final = J(Gstar, p_s, .)
            FamaMacbethIV(Y_s, X_s_all, Z_s_all, clusteringStar, b_coeffs_clusterwise_final)
            
            theta_C_final = b_coeffs_clusterwise_final[., iCov]
            Coef = mean(theta_C_final)
            SE = .
            if (Gstar > 1) SE = sqrt(variance(theta_C_final)/Gstar)
            
            tVal = . ; pStar = . ; pValAdj = . ; CI_lower = . ; CI_upper = .
            if (SE != . & SE > 1e-9) {
                tVal = Coef/SE
                pStar = 2*t(Gstar-1,-abs(tVal))
                if (Gstar > 1) {
                     pValAdj = mean(pStar :>= pValSim[.,indStar])
                     if(missing(pValAdj)) pValAdj = pStar // Fallback
                } else {
                     pValAdj = pStar
                }
                gap = -invt(Gstar-1,sigLevelAdjVec[indStar]/2)*SE
                CI_lower = Coef-gap
                CI_upper = Coef+gap
            } else { // Gstar=1 or SE is zero/missing
                pValAdj = (Gstar==1 ? 1 : (abs(Coef)>1e-9 ? 0 : 1) ) // crude pval for G=1
            }
            
            resultsMat[iCov,.] = (Coef,SE,tVal,pValAdj,CI_lower,CI_upper,Gstar)
        } // end iCov loop

        st_numscalar("n_obs", n)
        st_matrix("resultsMat",resultsMat)
    }

    mata{ // For R2 and RootMSE - use full sample 2SLS results
        betaHat_final_2sls = beta_s_2sls // from earlier full sample
        Yhat_final   = X_s_all * betaHat_final_2sls
        resid_final  = Y_s - Yhat_final
        RSS_final    = sum((resid_final:^2))
        TSS_final    = sum((Y_s :- mean(Y_s)):^2)
        R2_final     = 1 - RSS_final/TSS_final
        // df for IV: n - number of parameters. Here p_s.
        RootMSE_final= sqrt(RSS_final/(n-p_s))
        st_numscalar("R2", R2_final)
        st_numscalar("RootMSE", RootMSE_final)
    }

    mat colnames resultsMat = Coefficient Std_err t PValue CI_lower CI_upper Gstar
    // Construct row names: endog_var, exog_vars, _cons
    local rnames `all_struc_regressors' _cons
    mat rownames resultsMat = `rnames'


    mata: st_local("R2_val", strofreal(R2_final))
    mata: st_local("RootMSE_val", strofreal(RootMSE_final))

    local obs_text "Number of obs  = "
    local obs_num = string(n_obs, "%9.0g")
    local r2_text  "R^2 (2SLS)     = "
    local r2_num = string(`R2_val', "%9.4f")
    local rmse_text "Root MSE (2SLS)= "
    local rmse_num = string(`RootMSE_val', "%9.4f")

    di _n as text "Ibragimov and Muller (IV) with learned cluster"
    di _col(65) as text "`obs_text'" _col(5) as result "`obs_num'"
    di _col(65) as text "`r2_text'" _col(5) as result "`r2_num'"
    di _col(65) as text "`rmse_text'" _col(5) as result "`rmse_num'"
    di as text "{hline 85}"
    // Header for output table, depvar is `depvar'
    di as text %12s abbrev("`depvar'",12) _col(14) " {c |} Coefficient  Std. err.      t    P>|t|     [95% conf. interval]   Clusters"
    di as text "{hline 85}"

    // Loop to display results for each structural parameter
    local list_of_struc_param_names : list rnames
    forvalues i = 1/`=p_s' { // p_s is number of structural parameters
        local name : word `i' of `list_of_struc_param_names'
        local coef = resultsMat[`i',1]
        local se   = resultsMat[`i',2]
        local t    = resultsMat[`i',3]
        local p    = resultsMat[`i',4]
        local lci  = resultsMat[`i',5]
        local uci  = resultsMat[`i',6]
        local gval = resultsMat[`i',7]

        di as text %12s abbrev("`name'",12) ///
            _col(14) " {c |} " as res %8.6f `coef'  ///
            as res _col(28) %8.6f `se'  ///
            as res _col(38) %8.3f `t'    ///
            as res _col(48) %8.6f `p'    ///
            as res _col(58) %8.6f `lci' "   " %8.6f `uci' ///
            as res _col(78) %5.0f `gval'
    }
    di as text "{hline 85}"
end



mata:
// New FamaMacbethIV function
// Inputs: Y_s (outcome), X_s_all (structural regressors), Z_s_all (instruments for structural), 
//         clustering_labels, b_coeffs_clusterwise (output matrix G x p_s)
void FamaMacbethIV(Y_s_arg, X_s_all_arg, Z_s_all_arg, clustering_labels, b_coeffs_clusterwise_out) {
    real scalar G_fm, p_s_fm, n_fm
    G_fm = max(clustering_labels)
    p_s_fm = cols(X_s_all_arg)
    n_fm = rows(Y_s_arg)
    
    // Ensure output matrix has correct dimensions if passed pre-sized
    if (rows(b_coeffs_clusterwise_out) != G_fm | cols(b_coeffs_clusterwise_out) != p_s_fm) {
        b_coeffs_clusterwise_out = J(G_fm, p_s_fm, .)
    }

    for (ii=1; ii<=G_fm; ii++ ){
        fii = selectindex(clustering_labels :== ii)
        if (rows(fii) == 0) continue // Skip if cluster is empty

        Y_c = Y_s_arg[fii,.]
        X_s_c = X_s_all_arg[fii,.]
        Z_s_c = Z_s_all_arg[fii,.]
        
        // Check for sufficient observations and rank
        if (rows(Y_c) < cols(X_s_c)) { // Not enough obs
            b_coeffs_clusterwise_out[ii,.] = J(1, p_s_fm, .) // missing
            continue
        }
        
        // Use error handling with invsym()
        ZZ_inv = invsym(Z_s_c'*Z_s_c)
        if (ZZ_inv[1,1] == .) { // Check if inversion failed
            b_coeffs_clusterwise_out[ii,.] = J(1, p_s_fm, .)
            continue
        }
        
        XZ = X_s_c'*Z_s_c
        ZX = Z_s_c'*X_s_c
        ZY = Z_s_c'*Y_c
        
        // Compute the IV matrix components
        middle = invsym(XZ * ZZ_inv * ZX)
        if (middle[1,1] == .) { // Check if inversion failed
            b_coeffs_clusterwise_out[ii,.] = J(1, p_s_fm, .)
            continue
        }
        
        beta_c = middle * (XZ * ZZ_inv * ZY)
        b_coeffs_clusterwise_out[ii,.] = beta_c'
    }
}
end