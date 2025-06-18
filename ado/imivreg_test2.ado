*! imivreg version 2.0.0
program define imivreg_test2, eclass
    version 17

    syntax varlist(numeric min=1) [if] [in] [aw fw iw pw], ///
           CLuster(varlist numeric) [TIMEperiod(varname numeric)] ///
           IV(string asis)_test
    
    // (endog1 endog2 = inst1 inst2)
    local iv_content "`iv'"
    gettoken endog_vars iv_content : iv_content, parse("=")
    gettoken eq_sign iv_content : iv_content, parse("=")
    local instrument_vars "`iv_content'"

    if `"`eq_sign'"' != "=" | `"`endog_vars'"' == "" | `"`instrument_vars'"' == "" {
        di as error "iv() must be specified as endogenous_vars = instrument_vars"
        exit 198
    }
    di as text "Endogenous vars: `endog_vars'"
    di as text "Instruments: `instrument_vars'"
    
    tokenize `varlist'
    local depvar "`1'"
    macro shift
    local exog_vars `*'

    foreach v of local endog_vars {
        local present : list v in exog_vars
        if `present' {
            di as err "Endogenous variable (`v') cannot also be listed as an exogenous variable."
            exit 198
        }
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

    local all_struc_regressors `endog_vars' `exog_vars'

    mata{
        uniformseed(12345)
        rseed(12345)

        Y_s         = st_data(., "`depvar'")      
        X_e_mat     = st_data(., tokens("`endog_vars'"))
        if ("`exog_vars'" != "") {
            X_k     = st_data(., tokens("`exog_vars'"))
        } else {
            X_k     = J(rows(Y_s), 0, .)
        }
        Z_i         = st_data(., tokens("`instrument_vars'"))
        
        coord       = st_data(.,"`cluster'")
        timePeriod  = st_data(.,"`timeperiod'")
        n           = length(Y_s)
        m_endog     = cols(X_e_mat)
        q_inst      = cols(Z_i)

        X_s_all     = (X_e_mat, X_k, J(n,1,1))
        p_s         = cols(X_s_all)
        Z_s_all     = (Z_i, X_k, J(n,1,1))
        X_fs_all    = (Z_i, X_k, J(n,1,1))
        p_fs        = cols(X_fs_all)
        
        if (q_inst < m_endog) {
            errprintf("Error: Model is underidentified. Number of instruments (%f) must be at least the number of endogenous variables (%f).\n", q_inst, m_endog)
            exit(3498)
        }
        
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
        st_matrix("Y_s", Y_s); st_matrix("X_e_mat", X_e_mat); st_matrix("X_k", X_k); st_matrix("Z_i", Z_i);
        st_matrix("X_s_all", X_s_all); st_matrix("Z_s_all", Z_s_all); st_matrix("X_fs_all", X_fs_all);
        st_numscalar("p_s", p_s); st_numscalar("m_endog", m_endog); st_numscalar("n", n);
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
        Y_s = st_matrix("Y_s"); X_s_all = st_matrix("X_s_all"); Z_s_all = st_matrix("Z_s_all");
        X_e_mat = st_matrix("X_e_mat"); X_k = st_matrix("X_k"); X_fs_all = st_matrix("X_fs_all");
        p_s = st_numscalar("p_s"); m_endog = st_numscalar("m_endog"); n = st_numscalar("n");
        G_max = st_numscalar("G_max");
        if(l_G > 0) G_vec = range(2,G_max,1); else G_vec=J(0,1,.);

        if (l_G > 0) {
            data_medoids = st_data(.,.)
            ncols_med = cols(data_medoids)
            clusteringSet =  data_medoids[., (ncols_med-l_G+1)..ncols_med]
        } else {
            clusteringSet = J(n,0,.)
        }
        

        run_2sls(Y_s, X_s_all, Z_s_all, beta_s_2sls=., U_hat=.)
        
        
        pi_fs_ols_mat = J(cols(X_fs_all), m_endog, .)
        V_hat_mat     = J(n, m_endog, .)
        for (j=1; j<=m_endog; j++) {
            run_ols(X_e_mat[.,j], X_fs_all, pi_fs_ols_mat[.,j], V_hat_mat[.,j])
        }
        
        dis_mat = st_matrix("dis_mat"); time_mat = st_matrix("time_mat");
        M_ident = I(n); Qd_ident=Rd_ident=ex_ident=.; qrdp(M_ident, Qd_ident, Rd_ident, ex_ident); useQML_all = ex_ident[1..n]; 
        
        alphaHat_U = SigmaHat_U = .
        run_qmle(U_hat, hasTime, M_ident, useQML_all, dis_mat, time_mat, alphaHat_U, SigmaHat_U)
        
        SigmaHat_V_array = J(m_endog, 1, NULL) 
        for (j=1; j<=m_endog; j++) {
            alphaHat_V_j = SigmaHat_V_j =.
            run_qmle(V_hat_mat[.,j], hasTime, M_ident, useQML_all, dis_mat, time_mat, alphaHat_V_j, SigmaHat_V_j)
            SigmaHat_V_array[j] = &SigmaHat_V_j
        }
        
        
        all_residuals = (U_hat, V_hat_mat)
        
        whitened_residuals = J(n, m_endog+1, .)
        
        CSHat_U = cholesky(SigmaHat_U)
        whitened_residuals[.,1] = lusolve(CSHat_U, all_residuals[.,1])

        CSHat_V_array = J(m_endog, 1, NULL) 
        for (j=1; j<=m_endog; j++) {
            
            CSHat_V_j = cholesky(*SigmaHat_V_array[j])
            CSHat_V_array[j] = &CSHat_V_j
            whitened_residuals[.,j+1] = lusolve(CSHat_V_j, all_residuals[.,j+1])
        }
        
        rhoHat_mat = correlation(whitened_residuals)
        if (any(missing(rhoHat_mat))) rhoHat_mat = I(m_endog+1)

        
        Sigma_sim_block = J((m_endog+1)*n, (m_endog+1)*n, 0)
        all_CSHats = (&CSHat_U, CSHat_V_array)

        for (i_block=1; i_block<=m_endog+1; i_block++) {
            for (j_block=i_block; j_block<=m_endog+1; j_block++) {
                
                block_cov = rhoHat_mat[i_block, j_block] * (*all_CSHats[i_block]) * (*all_CSHats[j_block])'
                Sigma_sim_block[(i_block-1)*n+1..i_block*n, (j_block-1)*n+1..j_block*n] = block_cov
                if (i_block != j_block) {
                    Sigma_sim_block[(j_block-1)*n+1..j_block*n, (i_block-1)*n+1..i_block*n] = block_cov'
                }
            }
        }
        
        eigval = symeigensystem(Sigma_sim_block, ., .)
        min_eigval = min(eigval)
        if (missing(min_eigval) || min_eigval <= 1e-8) { 
            Sigma_sim_block = Sigma_sim_block + I((m_endog+1)*n)*max((1e-6, (min_eigval <= 1e-8 ? -min_eigval + 1e-6 : 1e-6)))
        }
        
        sigLevel = .05; Bboot = 1000;

        CSHat_sim_block = cholesky(Sigma_sim_block)
        UVbootMat = CSHat_sim_block' * rnormal((m_endog+1)*n, Bboot, 0, 1)
        
        resultsMat = J(p_s, 7, .)
        GstarVec = J(p_s, 1, .)

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
                
                theta_C_vec_sim = J(Bboot, 1, .) 
                se_of_mean_theta_C_vec_sim = J(Bboot, 1, .)

                for(rr = 1; rr <= Bboot; rr++){

                    U_boot_iter = UVbootMat[1..n, rr]
                    V_boot_iter_vec = UVbootMat[(n+1)..rows(UVbootMat), rr]
                    V_boot_iter_mat = J(n, m_endog, .)
                    for(j_v=1; j_v<=m_endog; j_v++) {
                        V_boot_iter_mat[., j_v] = V_boot_iter_vec[(j_v-1)*n+1..j_v*n, 1]
                    }
                    
                    X_e_boot_mat = X_fs_all * pi_fs_ols_mat + V_boot_iter_mat
                    X_s_all_boot = (X_e_boot_mat, X_k, J(n,1,1))
                    Y_s_boot = X_s_all_boot * beta_s_H0 + U_boot_iter
                    
                    b_coeffs_clusterwise_sim = J(G, p_s, .)
                    FamaMacbethIV(Y_s_boot, X_s_all_boot, Z_s_all, clustering, b_coeffs_clusterwise_sim)
                    
                    btemp_raw = b_coeffs_clusterwise_sim[.,iCov]
                    num_valid_clusters_sim = rows(btemp_raw) - missing(btemp_raw)

                    if (num_valid_clusters_sim > 0) theta_C_vec_sim[rr] = mean(btemp_raw)
                    if (num_valid_clusters_sim > 1) se_of_mean_theta_C_vec_sim[rr] = sd_total(btemp_raw)/sqrt(num_valid_clusters_sim)
                    
                    if (missing(se_of_mean_theta_C_vec_sim[rr]) || se_of_mean_theta_C_vec_sim[rr] <= 1e-9) {
                        pValSim[rr,kk] = 1
                    } else {
                        pValSim[rr,kk] = 2*t(G-1, -abs(theta_C_vec_sim[rr]/se_of_mean_theta_C_vec_sim[rr]))
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
                

                valid_power_idx = selectindex(!missing(theta_C_vec_sim) :& !missing(se_of_mean_theta_C_vec_sim) :& se_of_mean_theta_C_vec_sim :> 1e-9)
                
                if (rows(valid_power_idx) > 0) {

                     aboot_valid = theta_C_vec_sim[valid_power_idx, .]
                     sboot_valid = se_of_mean_theta_C_vec_sim[valid_power_idx, .]
                     power_sum_for_G = 0
                     for (alt_idx = 1; alt_idx <= nalt; alt_idx++) {
                         t_stats_alt = (aboot_valid :- alternatives[alt_idx]) :/ sboot_valid
                         pvals_for_alt = 2*t(max((1,G-1)), -abs(t_stats_alt))
                         power_sum_for_G = power_sum_for_G + mean(pvals_for_alt :< sigLevelAdjVec[kk])
                     }
                     simPowerVec[kk] = power_sum_for_G / nalt
                } else {
                    simPowerVec[kk] = 0
                }
            } 

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

            // Final estimation on original data with Gstar
            b_coeffs_clusterwise_final = J(Gstar, p_s, .)
            FamaMacbethIV(Y_s, X_s_all, Z_s_all, clusteringStar, b_coeffs_clusterwise_final)
            
            theta_C_final_raw = b_coeffs_clusterwise_final[., iCov]
            num_valid_clusters = rows(theta_C_final_raw) - missing(theta_C_final_raw)
            Coef = . ; SE = . ; tVal = . ; pStar = . ; pValAdj = . ; CI_lower = . ; CI_upper = .
            
            if (num_valid_clusters > 1) {
                valid_indices = selectindex(theta_C_final_raw :!= .)
                theta_C_final = theta_C_final_raw[valid_indices, .]
                
                Coef = mean(theta_C_final)
                SE   = sd_total(theta_C_final)/sqrt(num_valid_clusters)

                if (!missing(SE) && SE > 1e-9) {
                    tVal = Coef/SE
                    pStar = 2*t(Gstar-1,-abs(tVal))
                    
                    currentSigLevelAdj = (l_G > 0 && !missing(indStar) && indStar <= rows(sigLevelAdjVec)) ? sigLevelAdjVec[indStar] : sigLevel
                    if (missing(currentSigLevelAdj)) currentSigLevelAdj = sigLevel
                    
                    if(l_G > 0 && !missing(indStar)) {
                         pValAdj = mean(pStar :>= pValSim[.,indStar])
                    } else {
                         pValAdj = pStar
                    }
                    if (missing(pValAdj)) pValAdj = pStar
                    
                    gap = -invt(Gstar-1,currentSigLevelAdj/2)*SE
                    CI_lower = Coef-gap
                    CI_upper = Coef+gap
                }
            } else if (num_valid_clusters == 1) {
                Coef = mean(theta_C_final_raw)
            }
            
            resultsMat[iCov,.] = (Coef, SE, tVal, pValAdj, CI_lower, CI_upper, Gstar)
        } 

        st_numscalar("n_obs", n)
        st_matrix("resultsMat",resultsMat)
        st_matrix("beta_s_2sls", beta_s_2sls)
    }

    // --- Output section ---
    mata{ 
        beta_s_2sls = st_matrix("beta_s_2sls")
        if(!any(missing(beta_s_2sls))) {
            Y_s_out = st_data(.,"Y")
            X_s_all_out = st_data(., tokens("`all_struc_regressors' _cons"))
            Yhat_final   = X_s_all_out * beta_s_2sls
            resid_final  = Y_s_out - Yhat_final
            RSS_final    = sum((resid_final:^2))
            TSS_final    = sum((Y_s_out :- mean(Y_s_out)):^2)
            if (TSS_final > 1e-9) R2_final = 1 - RSS_final/TSS_final else R2_final = .
            p_s_out = cols(X_s_all_out)
            n_out = rows(Y_s_out)
            if ((n_out-p_s_out) > 0) RootMSE_final= sqrt(RSS_final/(n_out-p_s_out)) else RootMSE_final = .
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

    di _n as text "Ibragimov and Muller (IV-Multi) with learned cluster" 
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
void run_ols(y, x, b_out, e_out) {
    if (rows(x) < cols(x)) {
        b_out = J(cols(x),1,.)
        e_out = J(rows(y),1,.)
        return
    }
    real matrix bread_inv
    bread_inv = x'*x
    if (cond(bread_inv) > 1e12) { 
		b_out=J(cols(x),1,.); e_out=J(rows(y),1,.); return 
		}
    b_out = invsym(bread_inv)*x'*y
    e_out = y - x*b_out
}

void run_2sls(y, x, z, | b_out, e_out) {
     if (rows(x) < cols(x) || rows(z) < cols(z)) {
        b_out=J(cols(x),1,.); e_out=J(rows(y),1,.); return
    }
    real matrix bread_inv, bread
    bread_inv = x'*z*invsym(z'*z)*z'*x
    if (cond(bread_inv) > 1e12) {
		b_out=J(cols(x),1,.); e_out=J(rows(y),1,.); return 
		}
    bread = invsym(bread_inv)
    b_out = bread * (x'*z*invsym(z'*z)*z'*y)
    e_out = y - x*b_out
}

void run_qmle(resid_vec, hasTime, M, useQML, dis_mat, time_mat, |alphaHat, SigmaHat) {
    real scalar sig2_0, rho_0, tau_0
    real matrix S, d_vec_select, t_vec_select, d_med, t_med
    
    sig2_0 = log(max((1e-8, mean(resid_vec:^2))))
    d_vec_select = select(colshape(dis_mat,1), colshape(dis_mat:>0,1))
    if(rows(d_vec_select)==0) {
		d_med = 1
	} else {
		d_med  = mm_median(d_vec_select)
	}
    if(missing(d_med) || d_med <=0) d_med = 1; 
    rho_0  = log(max((1e-8, d_med/(-ln(0.30)))))
    
    S = optimize_init()
    if (hasTime) {
        t_vec_select = select(colshape(time_mat,1), colshape(time_mat:>0,1))
        if(rows(t_vec_select)==0) {
			t_med = 1
		} else {
			t_med = mm_median(t_vec_select)
		}
        if(missing(t_med) || t_med <=0) t_med = 1;
        tau_0  = log(max((1e-8, t_med/(-ln(0.30)))))
        optimize_init_evaluator(S,  &QMLE_new() )
        optimize_init_params(S, (sig2_0, rho_0, tau_0))
        optimize_init_argument(S, 5, resid_vec) 
    } else {
        optimize_init_evaluator(S,  &QMLE_bin() )
        optimize_init_params(S, (sig2_0, rho_0))
        optimize_init_argument(S, 4, resid_vec) 
    }
    optimize_init_which(S, "min"); optimize_init_evaluatortype(S, "d0")
    optimize_init_argument(S, 1, M); optimize_init_argument(S, 2, useQML)
    optimize_init_argument(S, 3, dis_mat)
    if (hasTime) optimize_init_argument(S, 4, time_mat)
    
    alphaHat = optimize(S)
    
    if (hasTime) { 
        SigmaHat = Sigma_func_DGP(alphaHat, dis_mat, time_mat)
    } else { 
        SigmaHat = Sigma_func_DGP_bin(alphaHat, dis_mat)
        SigmaHat = SigmaHat + I(rows(SigmaHat))*max((1e-8, mean(resid_vec:^2)*1e-2))
    }
}

void FamaMacbethIV(Y_s_arg, X_s_all_arg, Z_s_all_arg, clustering_labels, b_coeffs_clusterwise_out) {
    real scalar G_fm, p_s_fm
    G_fm = max(clustering_labels)
    if (missing(G_fm) || G_fm < 1) {
        G_fm = 0
    }
    p_s_fm = cols(X_s_all_arg)
    
    b_coeffs_clusterwise_out = J(G_fm, p_s_fm, .)

    for (ii=1; ii<=G_fm; ii++ ){
        fii = selectindex(clustering_labels :== ii)
        if (rows(fii) == 0) { 
            continue
        }

        Y_c = Y_s_arg[fii,.]
        X_s_c = X_s_all_arg[fii,.]
        Z_s_c = Z_s_all_arg[fii,.]
        
        if (rows(Y_c) < cols(X_s_c) || rows(Z_s_c) < cols(Z_s_c) || rows(Y_c) < cols(Z_s_c)) { 
            continue
        }
        if (cols(Z_s_c) == 0) {
            continue
        }
        
        real matrix ZZ, ZZ_inv, XZ, ZX, ZY, middle_term, XZ_ZZinv, Xprime_PZ_X, beta_c

        ZZ = Z_s_c'*Z_s_c
        if (rows(ZZ) == 0 || !issquare(ZZ)) {
            continue
        }
        ZZ_inv = invsym(ZZ)
        if (isscalar(ZZ_inv) && ZZ_inv[1,1] == .) { 
            continue
        }
        
        XZ = X_s_c'*Z_s_c
        ZX = Z_s_c'*X_s_c // This is XZ'
        ZY = Z_s_c'*Y_c
        
        XZ_ZZinv = XZ * ZZ_inv
        Xprime_PZ_X = XZ_ZZinv * ZX 

        if (rows(Xprime_PZ_X) == 0 || !issquare(Xprime_PZ_X)) {
            continue
        }
        middle_term = invsym(Xprime_PZ_X)
        if (isscalar(middle_term) && middle_term[1,1] == .) { 
            continue
        }
        
        beta_c = middle_term * (XZ_ZZinv * ZY)
        
        if (isscalar(beta_c) && missing(beta_c)) {
            continue;
        }
        b_coeffs_clusterwise_out[ii,.] = beta_c'
    }
}
function logdet(A){
		Ldecomposition=Udecomposition =pdecomposition=.
		lud(A,Ldecomposition,Udecomposition,pdecomposition)
		Pdecomposition = I(rows(Ldecomposition))[pdecomposition,.]
		du = diagonal(Udecomposition)
		prod = 1
		for (i=1;i<=rows(du);i++){
			prod =prod * sign(du[i])
		}
		c = det(Pdecomposition) * prod
		v = log(c) + sum(log(abs(du)))
		return(v)
	}

void QMLE_new(todo, w, M, useQML, dis_mat, time_mat, resid, Q, grad, hessian) {
    Sigma_func = exp(w[1]) :* exp(-dis_mat/exp(w[2])) :* exp(-time_mat/exp(w[3]))
    if (rows(M) > 0) {
        Sigma_func = M[useQML,.] * Sigma_func * (M[useQML,.])'
    }
    R = cholesky(Sigma_func)
    invSigma_resid = lusolve(R, lusolve(R', resid[useQML,1]))
    Q = 0.5*logdet(Sigma_func) + 0.5*quadcolsum(resid[useQML,1] :* invSigma_resid)
}

function Sigma_func_DGP(w,dis_mat,time_mat){
		SigmaHat = exp(w[1])*exp(-dis_mat/exp(w[2])-time_mat/exp(w[3]))
		return(SigmaHat)
	}

void FamaMacbeth(D,X,Y,Z,index,b,se){
		X_mat = D,X
		Z_mat = Z,X
		k = cols(X_mat)
		G = rows(uniqrows(index))
		btemp = J(G,k,0)
		for (ii=1; ii<=G; ii++ ){
			fii = index:== ii 
			temp = select(Z_mat,fii)'* select(X_mat,fii)
			ktemp= select(Z_mat,fii)'*select(Y,fii)
			btemp[ii,.] = (invsym(temp)*ktemp)'
		}
		b = mean(btemp)
		se = (diagonal(sqrt(variance(btemp)))/sqrt(G))'
	}

void QMLE_bin(todo, w, M, useQML, dis_mat, resid, Q, grad, H)
{
    real scalar s2, rho, eps, pen, g1, g2
    s2  = exp(w[1])
    rho = exp(w[2])
    real matrix Kfull, Sigma_full
    Kfull      = exp(-dis_mat / rho)
    Sigma_full = s2 :* Kfull
    real scalar haveM
    haveM = (rows(M) > 0)
    real matrix Sigma, Msub
    real colvector rsub
    rsub = resid[useQML, 1]
    if (haveM) {
        Msub  = M[useQML, .]
        Sigma = Msub * Sigma_full * Msub'
    }
    else {
        Sigma = Sigma_full[useQML, useQML]
    }
    eps   = mean(rsub:^2) * 1e-2
    Sigma = Sigma + I(rows(Sigma)) * eps
    pen = 0
    if (abs(w[1]) > 20) pen = pen + (abs(w[1]) - 20)^2
    if (abs(w[2]) > 10) pen = pen + (abs(w[2]) - 10)^2
    real matrix R
    R      = cholesky(Sigma)
    rsub   = lusolve(R', rsub)
    rsub   = lusolve(R , rsub)
    Q      = 0.5*logdet(Sigma) + 0.5*quadcolsum(resid[useQML,1] :* rsub) + pen

    if (args() >= 8) {
        real matrix dS1_full, dS2_full, dK_drho
        dS1_full = Sigma_full
        dK_drho  = (dis_mat:/rho) :* Kfull
        dS2_full = s2 :* dK_drho
        real matrix dS1, dS2
        if (haveM) {
            dS1 = Msub * dS1_full * Msub'
            dS2 = Msub * dS2_full * Msub'
        }
        else {
            dS1 = dS1_full[useQML, useQML]
            dS2 = dS2_full[useQML, useQML]
        }
        real matrix invS_d1, invS_d2
        invS_d1 = lusolve(R , lusolve(R', dS1))
        invS_d2 = lusolve(R , lusolve(R', dS2))
        g1 = 0.5*(trace(invS_d1) - quadcolsum(rsub :* (invS_d1 * rsub)))
        g2 = 0.5*(trace(invS_d2) - quadcolsum(rsub :* (invS_d2 * rsub)))
        if (abs(w[1]) > 20) g1 = g1 + 2*sign(w[1])*(abs(w[1]) - 20)
        if (abs(w[2]) > 10) g2 = g2 + 2*sign(w[2])*(abs(w[2]) - 10)
        grad = (g1, g2)
    }
}

function Sigma_func_DGP_bin(w, dis_mat)
{
    return( exp(w[1]) :* exp(-dis_mat/exp(w[2])) )
}

void FamaMacbethCRS(D, X, Y, Z, index, btemp){
    X_mat = D, X
    Z_mat = Z, X
    k = cols(X_mat)
    G = rows(uniqrows(index))
    btemp = J(G, 1, 0)
    for (ii = 1; ii <= G; ii++){
        fii = index :== ii 
        temp = select(Z_mat, fii)' * select(X_mat, fii)
        ktemp = select(Z_mat, fii)' * select(Y, fii)
        btemp[ii, 1] = (invsym(temp) * ktemp)[1,1]
    }
}

real scalar sd_total(matrix X) {
    n = rows(X) * cols(X)
    return(sqrt(variance(colshape(X, 1), 1)))  //
}

function cluster_se(x, e, XpXinv, group, |k){
    n = rows(e)
    if (args() < 5){
        k = rows(XpXinv)
    }
    k = rows(XpXinv)
    V = J(k, k, 0)
    for(ii = 1; ii <= max(group); ii++){
        I = group :!= ii
        V = V + (select(x, I)' * select(e, I)) * (select(x, I)' * select(e, I))'
    }
    vcluster = ((n-1)/(n-k)) * (max(group)/(max(group)-1)) * XpXinv * V * XpXinv'
    se = sqrt(diagonal(vcluster))
    return(se)
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