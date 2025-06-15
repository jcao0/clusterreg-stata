// cceivreg.ado
program define cceivreg_test, eclass
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
    
    tokenize `varlist'
    local depvar "`1'"
    macro shift
    local exog_vars `*'

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
    local all_struc_regressors `endog_var' `exog_vars'

    mata{
        uniformseed(12345)
        rseed(12345)

        Y_s         = st_data(., "`depvar'")      
        X_e         = st_data(., "`endog_var'")
        if ("`exog_vars'" != "") {
            X_k     = st_data(., tokens("`exog_vars'"))
        } else {
            X_k     = J(rows(Y_s), 0, .)
        }
        Z_i         = st_data(., tokens("`instrument_vars'"))
        
        coord       = st_data(.,"`cluster'")
        timePeriod  = st_data(.,"`timeperiod'")
        n           = length(Y_s)

        X_s_all     = (X_e, X_k, J(n,1,1))
        p_s         = cols(X_s_all)

        Z_s_all     = (Z_i, X_k, J(n,1,1))
        X_fs_all    = (Z_i, X_k, J(n,1,1))

        G_max = ceil(n^(1/3))
        if (G_max < 2) G_max = 2
        st_numscalar("G_max", G_max)
        
        if (G_max >= 2) {
            G_vec = range(2, G_max, 1)
            l_G = length(G_vec)
        } else {
            G_vec = J(0,1,.)
            l_G = 0
        }
        st_numscalar("l_G", l_G)
    }

    local G_max_val = G_max 
    if (`G_max_val' >= 2) {
        gen id0 = _n 
        matrix dissim dis_mat = `cluster',L2 
        matrix dissim time_mat = `timeperiod', L2
        forvalues i = 2/`G_max_val' { 
            qui clpam group`i', distmat(dis_mat) id(id0) medoids(`i') ga
        }
    }

    mata{
        
        hasTime = st_numscalar("`__htscalar'")
        l_G = st_numscalar("l_G")
        Y_s = st_data(.,"Y"); X_e=st_data(.,"`endog_var'"); X_k = st_data(.,tokens("`exog_vars'")); Z_i = st_data(.,tokens("`instrument_vars'"));
        X_s_all = (X_e, X_k, J(rows(Y_s),1,1)); Z_s_all = (Z_i, X_k, J(rows(Y_s),1,1)); X_fs_all = (Z_i, X_k, J(rows(Y_s),1,1));
        n = rows(Y_s); p_s = cols(X_s_all); G_max = st_numscalar("G_max");
        if(l_G > 0) G_vec = range(2,G_max,1);
        
        if (l_G > 0) {
            data_medoids = st_data(.,.)
            ncols_med = cols(data_medoids)
            clusteringSet =  data_medoids[., (ncols_med-l_G+1)..ncols_med]
        } else {
            clusteringSet = J(n,0,.)
        }

        // Initial Full Sample Estimations with error handling
        beta_s_2sls = J(p_s, 1, .)
        pi_fs_ols = J(cols(X_fs_all), 1, .)
        U_hat = J(n,1,.)
        V_hat = J(n,1,.)
        
        temp_ZtZ = Z_s_all' * Z_s_all
        if (rows(temp_ZtZ) > 0 && issquare(temp_ZtZ) && cond(temp_ZtZ) < 1e12) {
            temp_inv_ZtZ = invsym(temp_ZtZ)
            if (!(isscalar(temp_inv_ZtZ) && missing(temp_inv_ZtZ))) {
                temp_Xprime_PZX = X_s_all' * Z_s_all * temp_inv_ZtZ * Z_s_all' * X_s_all
                if (rows(temp_Xprime_PZX) > 0 && issquare(temp_Xprime_PZX) && cond(temp_Xprime_PZX) < 1e12) {
                    beta_s_2sls_inv = invsym(temp_Xprime_PZX)
                    if (!(isscalar(beta_s_2sls_inv) && missing(beta_s_2sls_inv))) {
                         beta_s_2sls = beta_s_2sls_inv * (X_s_all' * Z_s_all * temp_inv_ZtZ * Z_s_all' * Y_s)
                         U_hat = Y_s - X_s_all * beta_s_2sls
                    }
                }
            }
        }
        if (any(missing(beta_s_2sls))) {
            errprintf("Error: Could not compute initial 2SLS estimates. Matrix singularity or rank deficiency.\n")
            exit(3300)
        }
        
        temp_Xfs_prime_Xfs = X_fs_all' * X_fs_all
        if (rows(temp_Xfs_prime_Xfs) > 0 && issquare(temp_Xfs_prime_Xfs) && cond(temp_Xfs_prime_Xfs) < 1e12) {
            pi_fs_ols_inv = invsym(temp_Xfs_prime_Xfs)
            if (!(isscalar(pi_fs_ols_inv) && missing(pi_fs_ols_inv))) {
                pi_fs_ols = pi_fs_ols_inv * X_fs_all' * X_e
                V_hat = X_e - X_fs_all * pi_fs_ols
            }
        }
         if (any(missing(pi_fs_ols))) {
            errprintf("Error: Could not compute initial first-stage OLS estimates. Matrix singularity or rank deficiency.\n")
            exit(3300)
        }

        M_ident = I(n)
        Qd_ident=Rd_ident=ex_ident=.
        qrdp(M_ident, Qd_ident, Rd_ident, ex_ident) 
        useQML_all = ex_ident[1..n] 

        dis_mat = st_matrix("dis_mat")
        time_mat = st_matrix("time_mat")
        
        // --- QMLE for Sigma_U and Sigma_V ---
        // QMLE for Sigma_U
        sig2_0_U = log(max((1e-8, mean(U_hat:^2))))
        d_vec_U_select = select(colshape(dis_mat,1), colshape(dis_mat:>0,1))
        if(rows(d_vec_U_select)==0){
			d_med_U = 1
		} else {
			d_med_U  = mm_median(d_vec_U_select)
		}
        if(missing(d_med_U) || d_med_U <=0) d_med_U = 1; 
        rho_0_U  = log(max((1e-8, d_med_U/(-ln(0.30)))))
        if (hasTime) {
            t_vec_U_select = select(colshape(time_mat,1), colshape(time_mat:>0,1))
            if(rows(t_vec_U_select)==0) {
				t_med_U = 1
			} else {
				t_med_U = mm_median(t_vec_U_select)
			}
            if(missing(t_med_U) || t_med_U <=0) t_med_U = 1;
            tau_0_U  = log(max((1e-8, t_med_U/(-ln(0.30)))))
        }
        S_U = optimize_init()
        if (hasTime) {
            optimize_init_evaluator(S_U,  &QMLE_new() )
            optimize_init_params(S_U, (sig2_0_U, rho_0_U, tau_0_U))
            optimize_init_argument(S_U, 5, U_hat) 
        } else {
            optimize_init_evaluator(S_U,  &QMLE_bin() )
            optimize_init_params(S_U, (sig2_0_U, rho_0_U))
            optimize_init_argument(S_U, 4, U_hat) 
        }
        optimize_init_which(S_U, "min"); optimize_init_evaluatortype(S_U, "d0")
        optimize_init_argument(S_U, 1, M_ident); optimize_init_argument(S_U, 2, useQML_all)
        optimize_init_argument(S_U, 3, dis_mat)
        if (hasTime) optimize_init_argument(S_U, 4, time_mat)
        alphaHat_U = optimize(S_U)
        
        // QMLE for Sigma_V
        sig2_0_V = log(max((1e-8, mean(V_hat:^2))))
        rho_0_V = rho_0_U
        if (hasTime) tau_0_V = tau_0_U
        S_V = optimize_init()
        if (hasTime) {
            optimize_init_evaluator(S_V,  &QMLE_new() )
            optimize_init_params(S_V, (sig2_0_V, rho_0_V, tau_0_V))
            optimize_init_argument(S_V, 5, V_hat) 
        } else {
            optimize_init_evaluator(S_V,  &QMLE_bin() )
            optimize_init_params(S_V, (sig2_0_V, rho_0_V))
            optimize_init_argument(S_V, 4, V_hat) 
        }
        optimize_init_which(S_V, "min"); optimize_init_evaluatortype(S_V, "d0")
        optimize_init_argument(S_V, 1, M_ident); optimize_init_argument(S_V, 2, useQML_all)
        optimize_init_argument(S_V, 3, dis_mat)
        if (hasTime) optimize_init_argument(S_V, 4, time_mat)
        alphaHat_V = optimize(S_V)
        
        if (hasTime) { 
            SigmaHat_U = Sigma_func_DGP(alphaHat_U, dis_mat, time_mat)
            SigmaHat_V = Sigma_func_DGP(alphaHat_V, dis_mat, time_mat) 
        } else { 
            SigmaHat_U = Sigma_func_DGP_bin(alphaHat_U, dis_mat)
            SigmaHat_U = SigmaHat_U + I(rows(SigmaHat_U))*max((1e-8, mean(U_hat:^2)*1e-2))
            SigmaHat_V = Sigma_func_DGP_bin(alphaHat_V, dis_mat) 
            SigmaHat_V = SigmaHat_V + I(rows(SigmaHat_V))*max((1e-8, mean(V_hat:^2)*1e-2))
        }

        CSHat_U = cholesky(SigmaHat_U)
        CSHat_V = cholesky(SigmaHat_V)
        U_transformed = lusolve(CSHat_U, U_hat)
        V_transformed = lusolve(CSHat_V, V_hat)
        rhoHat = correlation(U_transformed, V_transformed)
        if (missing(rhoHat)) rhoHat = 0 

        Sigma_UV = rhoHat * CSHat_U * CSHat_V'
        Sigma_sim_block = J(2*n, 2*n, 0)
        Sigma_sim_block[1..n, 1..n] = SigmaHat_U
        Sigma_sim_block[(n+1)..(2*n), (n+1)..(2*n)] = SigmaHat_V
        Sigma_sim_block[1..n, (n+1)..(2*n)] = Sigma_UV
        Sigma_sim_block[(n+1)..(2*n), 1..n] = Sigma_UV'
        
        eigval = symeigensystem(Sigma_sim_block, ., .)
        min_eigval = min(eigval)
        if (missing(min_eigval) || min_eigval <= 1e-8) { 
            Sigma_sim_block = Sigma_sim_block + I(2*n)*max((1e-6, (min_eigval <= 1e-8 ? -min_eigval + 1e-6 : 1e-6)))
        }
        
        // --- Simulation setup ---
        sigLevel = .05
        Bboot = 1000 
        CSHat_sim_block = cholesky(Sigma_sim_block)
        UVbootMat = CSHat_sim_block' * rnormal(2*n, Bboot, 0, 1)
        
        resultsMat = J(p_s, 7, .)
        GstarVec = J(p_s, 1, .)

        // =========================================================
        // === CORE LOGIC CHANGE: CCE Method Implementation      ===
        // =========================================================
        for(iCov = 1; iCov <= p_s; iCov++ ){
            beta_s_H0 = beta_s_2sls 
            beta_s_H0[iCov] = 0     
            
            simPowerVec = J(l_G,1,0)
            pValSim = J(Bboot,l_G,1) 
            sigLevelAdjVec = J(l_G,1,sigLevel)

            if (l_G == 0) {
                resultsMat[iCov,.] = (beta_s_2sls[iCov], ., ., 1, ., ., G_max)
                continue
            }
            
            for(kk = 1; kk <= l_G; kk++){
                clustering = clusteringSet[.,kk]
                G = G_vec[kk]
                
                abootVec = J(Bboot,1,.)
                sbootVec = J(Bboot,1,.)

                for(rr = 1; rr <= Bboot; rr++){
                    U_boot_iter = UVbootMat[1..n, rr]
                    V_boot_iter = UVbootMat[(n+1)..(2*n), rr]
                    
                    X_e_boot = X_fs_all * pi_fs_ols + V_boot_iter
                    X_s_all_boot = X_s_all 
                    X_s_all_boot[.,1] = X_e_boot 
                    Y_s_boot = X_s_all_boot * beta_s_H0 + U_boot_iter
                    
                    // CCE uses full-sample estimates
                    bh_boot = invsym(X_s_all_boot'*Z_s_all*invsym(Z_s_all'*Z_s_all)*Z_s_all'*X_s_all_boot) * (X_s_all_boot'*Z_s_all*invsym(Z_s_all'*Z_s_all)*Z_s_all'*Y_s_boot)
                    resid_boot = Y_s_boot - X_s_all_boot * bh_boot
                    
                    // Get cluster robust SE for the simulated full-sample 2SLS
                    se_boot_vec = cluster_se_iv(X_s_all_boot, Z_s_all, resid_boot, clustering)

                    abootVec[rr] = bh_boot[iCov]
                    sbootVec[rr] = se_boot_vec[iCov]

                    if (missing(sbootVec[rr]) || sbootVec[rr] <= 1e-9) {
                        pValSim[rr,kk] = 1
                    } else {
                        t_boot = abootVec[rr]/sbootVec[rr]
                        t_boot_adj = t_boot * sqrt((G-1)/G)
                        pValSim[rr,kk] = 2*t(G-1, -abs(t_boot_adj))
                    }
                } 

                if (G > 1 && rows(selectindex(!missing(pValSim[.,kk]))) > 0) {
                    sigLevelAdjVec[kk] = min((mm_quantile(pValSim[.,kk],1,0.05),0.05))
                    if (missing(sigLevelAdjVec[kk])) sigLevelAdjVec[kk] = 0.05 
                } else {
                    sigLevelAdjVec[kk] = 0.05 
                }
                
                alt_neg_range = range(-10, -1, 1); alt_pos_range = range(1, 10, 1)   
                alternatives = (alt_neg_range \ alt_pos_range) :/ sqrt(n)
                nalt = rows(alternatives)
                
                valid_power_idx = selectindex(!missing(abootVec) :& !missing(sbootVec) :& sbootVec :> 1e-9)
                
                if (rows(valid_power_idx) > 0) {
                     aboot_valid = abootVec[valid_power_idx, .]
                     sboot_valid = sbootVec[valid_power_idx, .]
                     power_sum_for_G = 0
                     for (alt_idx = 1; alt_idx <= nalt; alt_idx++) {
                         t_stats_alt = (aboot_valid :- alternatives[alt_idx]) :/ sboot_valid
                         t_stats_alt_adj = t_stats_alt * sqrt((G-1)/G)
                         pvals_for_alt = 2*t(max((1,G-1)), -abs(t_stats_alt_adj))
                         power_sum_for_G = power_sum_for_G + mean(pvals_for_alt :< sigLevelAdjVec[kk])
                     }
                     simPowerVec[kk] = power_sum_for_G / nalt
                } else {
                    simPowerVec[kk] = 0
                }
            } // end kk loop

            indStar = windStar = .
            if (l_G > 0) {
                maxindex(simPowerVec,1,indStar,windStar)
                if (rows(indStar) > 1) indStar = indStar[1]
                if (missing(indStar)) indStar = l_G
            } else {
                indStar = . 
            }
            Gstar = missing(indStar) ? G_max : G_vec[indStar]
            if (missing(Gstar) || Gstar < 2) Gstar = G_max
            if (Gstar < 2 && G_max >=2) Gstar = 2; else if (Gstar < 2 && G_max < 2) Gstar = G_max;
            
            if (l_G > 0 && !missing(indStar)) {
                 clusteringStar = clusteringSet[.,indStar]
            } else { 
                 if (cols(clusteringSet) > 0) clusteringStar = clusteringSet[., cols(clusteringSet)];
                 else clusteringStar = J(n,1,1);
                 if (Gstar < 2 && cols(clusteringSet)==0) Gstar = 1;
            }
            GstarVec[iCov] = Gstar

            // --- Final CCE results calculation ---
            Coef_full = beta_s_2sls[iCov]
            // Calculate final full-sample residuals for SE calculation
            resid_full = Y_s - X_s_all * beta_s_2sls
            se_final_vec = cluster_se_iv(X_s_all, Z_s_all, resid_full, clusteringStar)
            SE_final = se_final_vec[iCov]

            tVal = . ; pStar = . ; pValAdj = . ; CI_lower = . ; CI_upper = .
            if (!missing(SE_final) && SE_final > 1e-9) {
                tVal = Coef_full/SE_final
                tVal_adj = tVal * sqrt((Gstar-1)/Gstar)
                pStar = 2*t(Gstar-1, -abs(tVal_adj))

                currentSigLevelAdj = (l_G > 0 && !missing(indStar) && indStar <= rows(sigLevelAdjVec)) ? sigLevelAdjVec[indStar] : sigLevel
                if(missing(currentSigLevelAdj)) currentSigLevelAdj = sigLevel

                if (l_G > 0 && !missing(indStar)) {
                    pValAdj = mean(pStar :>= pValSim[.,indStar])
                } else {
                    pValAdj = pStar
                }
                if (missing(pValAdj)) pValAdj = pStar

                gap = -invt(Gstar-1, currentSigLevelAdj/2) * SE_final * sqrt(Gstar/(Gstar-1)) // Adjust CI for t_dist
                CI_lower = Coef_full - gap
                CI_upper = Coef_full + gap
            }
            
            resultsMat[iCov,.] = (Coef_full, SE_final, tVal, pValAdj, CI_lower, CI_upper, Gstar)
        } // end iCov loop

        st_numscalar("n_obs", n)
        st_matrix("resultsMat",resultsMat)
    }

    // --- Output section (identical to crsivreg.ado, with title change) ---
    mata{ 
        if(!any(missing(beta_s_2sls))) {
            betaHat_final_2sls = st_matrix("beta_s_2sls")
            Yhat_final   = st_matrix("X_s_all") * betaHat_final_2sls
            resid_final  = st_matrix("Y_s") - Yhat_final
            RSS_final    = sum((resid_final:^2))
            TSS_final    = sum((st_matrix("Y_s") :- mean(st_matrix("Y_s"))):^2)
            if (TSS_final > 1e-9) {
				R2_final = 1 - RSS_final/TSS_final
			} else {
				R2_final = .
			}
            if ((n-p_s) > 0) {
				RootMSE_final= sqrt(RSS_final/(n-p_s))
			} else {
				RootMSE_final = .
			}
        } else {
            R2_final = .
            RootMSE_final = .
        }
        st_numscalar("R2", R2_final)
        st_numscalar("RootMSE", RootMSE_final)
    }

    mat colnames resultsMat = Coefficient Std_err t PValue CI_lower CI_upper Gstar
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

    di _n as text "CCE method (IV) with learned clusters" // Title changed
    di _col(65) as text "`obs_text'" _col(5) as result "`obs_num'"
    di _col(65) as text "`r2_text'" _col(5) as result "`r2_num'"
    di _col(65) as text "`rmse_text'" _col(5) as result "`rmse_num'"
    di as text "{hline 85}"
    di as text %12s abbrev("`depvar'",12) _col(14) " {c |} Coefficient  Std. err.      t    P>|t|     [95% conf. interval]   Clusters"
    di as text "{hline 85}"

    local list_of_struc_param_names : rownames resultsMat
    local p_s_count : word count `list_of_struc_param_names'
    forvalues i = 1/`p_s_count' { 
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


// =============================================================================
// Mata Helper Functions Block
// =============================================================================
mata:

function cluster_se_iv(X, Z, e, group) {
    real scalar k, n, G
    n = rows(e)
    k = cols(X)
    G = max(group)
    if (missing(G) || G==0) return(J(k,1,.))

    // 1. (X'PzX)^-1
    real matrix bread_inv, bread
    bread_inv = X'*Z*invsym(Z'*Z)*Z'*X
    bread = invsym(bread_inv)

    // 2. sum_g(S_g * S_g') where S_g = sum_{i in g} Z_i*e_i
    real matrix meat_sum
    meat_sum = J(cols(Z), cols(Z), 0)
    
    real matrix Z_g, e_g, S_g
    real colvector fii

    for (ii=1; ii<=G; ii++) {
        fii = selectindex(group:==ii)
        if (rows(fii)==0) continue
        
        Z_g = Z[fii,.]
        e_g = e[fii,.]
        S_g = Z_g' * e_g 
        meat_sum = meat_sum + S_g * S_g'
    }
    
    // 3. (V)
    real matrix middle_term, vcluster_mat
    middle_term = X'*Z*invsym(Z'*Z) * meat_sum * invsym(Z'*Z)*Z'*X
    vcluster_mat = bread * middle_term * bread
    

    if (G > 1) {
        vcluster_mat = vcluster_mat * (G/(G-1)) * ((n-1)/(n-k))
    }
    
    real colvector se
    se = sqrt(diagonal(vcluster_mat))
    return(se)
}


// FamaMacbethIV as refined in previous interactions
void FamaMacbethIV(Y_s_arg, X_s_all_arg, Z_s_all_arg, clustering_labels, b_coeffs_clusterwise_out) {
    real scalar G_fm, p_s_fm 
    G_fm = max(clustering_labels)
    // Robust handling of G_fm if no valid clusters
    if (missing(G_fm) || G_fm < 1) {
        G_fm = 0 // No loops will run if G_fm is 0
    }
    p_s_fm = cols(X_s_all_arg)
    
    b_coeffs_clusterwise_out = J(G_fm, p_s_fm, .) // Initialize with missings

    for (ii=1; ii<=G_fm; ii++ ){
        fii = selectindex(clustering_labels :== ii)
        if (rows(fii) == 0) { 
            // b_coeffs_clusterwise_out[ii,.] is already missing
            continue
        }

        Y_c = Y_s_arg[fii,.]
        X_s_c = X_s_all_arg[fii,.]
        Z_s_c = Z_s_all_arg[fii,.]
        
        if (rows(Y_c) < cols(X_s_c) || rows(Z_s_c) < cols(Z_s_c) || rows(Y_c) < cols(Z_s_c) /*added this last one too*/) { 
            // b_coeffs_clusterwise_out[ii,.] is already missing
            continue
        }
        // Ensure Z_s_c is not empty if it's expected to have columns (e.g. for constant)
        if (cols(Z_s_c) == 0) { // If no instruments defined for the cluster (e.g. all Z_i and X_k are empty after selection)
             // b_coeffs_clusterwise_out[ii,.] is already missing
            continue
        }
        
        real matrix ZZ, ZZ_inv, XZ, ZX, ZY, middle_term, XZ_ZZinv, Xprime_PZ_X, beta_c

        ZZ = Z_s_c'*Z_s_c
        // Check if ZZ is square and not empty before invsym
        if (rows(ZZ) == 0 || !issquare(ZZ)) {
             // b_coeffs_clusterwise_out[ii,.] is already missing
            continue
        }
        ZZ_inv = invsym(ZZ)
        if (isscalar(ZZ_inv) && ZZ_inv[1,1] == .) { 
            // b_coeffs_clusterwise_out[ii,.] is already missing
            continue
        }
        
        XZ = X_s_c'*Z_s_c
        ZX = Z_s_c'*X_s_c // This is XZ'
        ZY = Z_s_c'*Y_c
        
        XZ_ZZinv = XZ * ZZ_inv
        Xprime_PZ_X = XZ_ZZinv * ZX 

        // Check if Xprime_PZ_X is square and not empty
        if (rows(Xprime_PZ_X) == 0 || !issquare(Xprime_PZ_X)) {
            // b_coeffs_clusterwise_out[ii,.] is already missing
            continue
        }
        middle_term = invsym(Xprime_PZ_X) // Renamed from 'middle' for clarity
        if (isscalar(middle_term) && middle_term[1,1] == .) { 
            // b_coeffs_clusterwise_out[ii,.] is already missing
            continue
        }
        
        beta_c = middle_term * (XZ_ZZinv * ZY)
        
        if (isscalar(beta_c) && missing(beta_c)) { // Final check on result
            // b_coeffs_clusterwise_out[ii,.] is already missing
            continue;
        }
        b_coeffs_clusterwise_out[ii,.] = beta_c'
    }
}

// QMLE_new with internal stability ridge
void QMLE_new(todo, w, M, useQML, dis_mat, time_mat, resid, Q, grad, hessian) {
    real matrix Sigma_func, R_chol, invSigma_resid_vec
    Sigma_func = exp(w[1]) :* exp(-dis_mat/exp(w[2])) :* exp(-time_mat/exp(w[3]))
    
    // Stability modification
    Sigma_func = Sigma_func + I(rows(Sigma_func)) * 1e-9 
    
    real colvector resid_vec
    resid_vec = resid[useQML,1]

    R_chol = cholesky(Sigma_func)
    if (isscalar(R_chol) && missing(R_chol)) { 
        Q = 1e100 
        return
    }
    
    invSigma_resid_vec = lusolve(R_chol, lusolve(R_chol', resid_vec))
    real scalar log_det_Sigma
    log_det_Sigma = 2*sum(log(diagonal(R_chol))) 

    Q = 0.5*log_det_Sigma + 0.5*quadcolsum(resid_vec :* invSigma_resid_vec)
    if (missing(Q)) Q = 1e100;
}

// QMLE_bin with internal stability ridge
void QMLE_bin(todo, w, M, useQML, dis_mat, resid, Q, grad, H) {
    real scalar s2, rho, pen
    s2  = exp(w[1])
    rho = exp(w[2])
    real matrix Kfull, Sigma_for_QMLE, R_chol, invSigma_resid_vec
    Kfull      = exp(-dis_mat / rho)
    Sigma_for_QMLE = s2 :* Kfull
    
    // Stability modification
    Sigma_for_QMLE = Sigma_for_QMLE + I(rows(Sigma_for_QMLE)) * 1e-9 
    
    pen = 0
    if (abs(w[1]) > 20) pen = pen + (abs(w[1]) - 20)^2
    if (abs(w[2]) > 10) pen = pen + (abs(w[2]) - 10)^2

    real colvector resid_vec
    resid_vec = resid[useQML,1]

    R_chol = cholesky(Sigma_for_QMLE)
    if (isscalar(R_chol) && missing(R_chol)) {
        Q = 1e100 + pen
        return
    }

    invSigma_resid_vec = lusolve(R_chol, lusolve(R_chol', resid_vec))
    real scalar log_det_Sigma
    log_det_Sigma = 2*sum(log(diagonal(R_chol)))

    Q = 0.5*log_det_Sigma + 0.5*quadcolsum(resid_vec :* invSigma_resid_vec) + pen
    if (missing(Q)) Q = 1e100 + pen;
}

// Other standard helper functions
function logdet(A){
    if (rows(A)==0 || cols(A)==0) return(.)
    L=U=P=.
    lud(A,L,U,P)
    if (isdiagonal(U) && sum(diagonal(U):==0) > 0 ) return(.)
    c = det(P)
    v = sum(log(abs(diagonal(U))))
    if (c == -1) return(.);
    v = log(c) + v;
    return(v)
}

function Sigma_func_DGP(w,dis_mat,time_mat){
    return(exp(w[1]):*exp(-dis_mat/exp(w[2])-time_mat/exp(w[3])))
}

function Sigma_func_DGP_bin(w, dis_mat){
    return( exp(w[1]) :* exp(-dis_mat/exp(w[2])) )
}

real scalar sd_total(matrix X) {
    if (rows(X)==0 || cols(X)==0) return(.)
    n_elem = rows(X) * cols(X)
    if (n_elem <= 1) return(.)
    return(sqrt(variance(colshape(X, 1), 1)))  
}

real scalar issquare(real matrix A)
{
    return(rows(A)==cols(A) & rows(A)>0)
}

real scalar isscalar(real matrix X)
{
    return(rows(X)==1 & cols(X)==1)
}

end