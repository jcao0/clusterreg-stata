*! cceivreg.ado
*! Version 2.1: Supports multiple endogenous and instrumental variables.
program define cceivreg, eclass
    version 17
    // 1. Syntax Parsing for Multiple IVs (identical to imivreg/crsivreg)
    syntax varlist(numeric min=1) [if] [in] [aw fw iw pw], ///
        CLuster(varlist numeric) [TIMEperiod(varname numeric)] ///
        IV(string asis)

    // --- This entire setup section is identical to imivreg and crsivreg ---
    gettoken depvar main_vars : varlist
    local exog_vars `main_vars'
    local iv_content `iv'
    local eq_pos = strpos("`iv_content'", "=")
    if (`eq_pos' == 0) {
        di as error "iv() must be specified as endogenous_vars = instrument_vars"
        exit 198
    }
    local endog_vars = trim(substr("`iv_content'", 1, `eq_pos' - 1))
    local instrument_vars = trim(substr("`iv_content'", `eq_pos' + 1, .))

    di as text "Dependent var: `depvar'"
    di as text "Exogenous vars: `exog_vars'"
    di as text "Endogenous vars: `endog_vars'"
    di as text "Instruments: `instrument_vars'"

    foreach var of local endog_vars {
        local present : list var in exog_vars
        if `present' {
            di as err "Endogenous variable (`var') cannot also be listed as an exogenous variable."
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

    // Start Mata block for setup
    mata {
        uniformseed(12345)
        rseed(12345)
        Y_s         = st_data(., "`depvar'")
        if ("`endog_vars'" != "") { 
			X_e = st_data(., tokens("`endog_vars'")) 
			} else { 
				X_e = J(rows(Y_s), 0, .) 
				}
        if ("`exog_vars'" != "") { 
			X_k = st_data(., tokens("`exog_vars'")) 
			} else { 
				X_k = J(rows(Y_s), 0, .) 
				}
        if ("`instrument_vars'" != "") { 
			Z_i = st_data(., tokens("`instrument_vars'")) 
			} else { 
				Z_i = J(rows(Y_s), 0, .) 
				}
        
        coord       = st_data(.,"`cluster'")
        timePeriod  = st_data(.,"`timeperiod'")
        n           = length(Y_s)
        X_s_all     = (X_e, X_k, J(n,1,1))
        p_s         = cols(X_s_all) 
        m_e         = cols(X_e)
        Z_s_all     = (Z_i, X_k, J(n,1,1))
        X_fs_all    = Z_s_all
        
        G_max = ceil(n^(1/3))
        if(G_max < 2) G_max = 2
        st_numscalar("G_max", G_max)
        G_vec = range(2, G_max, 1)
        l_G = length(G_vec)
        st_numscalar("l_G", l_G)
    }

    // Cluster Generation
    gen id0 = _n 
    matrix dissim dis_mat = `cluster',L2 
    matrix dissim time_mat = `timeperiod', L2
    forvalues i =2/`=G_max' { 
        qui clpam group`i', distmat(dis_mat) id(id0) medoids(`i') ga
    }

    // Main Mata Block for Estimation
    mata {
        hasTime = st_numscalar("`__htscalar'")
        data_medoids = st_data(.,.)
        ncols_med = cols(data_medoids)
        clusteringSet =  data_medoids[.,ncols_med-G_max+2..ncols_med]

        // --- 3a. Initial Estimations & Residuals (identical to imivreg/crsivreg) ---
        beta_s_2sls = invsym(X_s_all'*Z_s_all*invsym(Z_s_all'*Z_s_all)*Z_s_all'*X_s_all) * (X_s_all'*Z_s_all*invsym(Z_s_all'*Z_s_all)*Z_s_all'*Y_s)
        U_hat = Y_s - X_s_all * beta_s_2sls
        pi_fs_ols = invsym(X_fs_all'*X_fs_all) * X_fs_all' * X_e
        V_hat = X_e - X_fs_all * pi_fs_ols
        all_residuals = (U_hat, V_hat)
        num_resids = cols(all_residuals)

        // --- 3b. Robust QMLE (identical to imivreg/crsivreg) ---
        M_ident = I(n)
        Qd_ident=Rd_ident=ex_ident=.; qrdp(M_ident, Qd_ident, Rd_ident, ex_ident); useQML_all = ex_ident[1..n]
        dis_mat = st_matrix("dis_mat"); time_mat = st_matrix("time_mat")
		SigmaHats = CSHats = J(num_resids, 1, NULL)
        for(j=1; j<=num_resids; j++){
            current_resid = all_residuals[., j]
            sig2_0 = log(max((1e-8, mean(current_resid:^2))))
            d_vec_select = select(colshape(dis_mat,1), colshape(dis_mat:>0,1))
            d_med = rows(d_vec_select)>0 ? mm_median(d_vec_select) : 1
            if(missing(d_med) || d_med<=0) d_med = 1
            rho_0 = log(max((1e-8, d_med/(-ln(0.30)))))
            S_j = optimize_init(); optimize_init_which(S_j, "min"); optimize_init_evaluatortype(S_j, "d0")
            optimize_init_argument(S_j, 1, M_ident); optimize_init_argument(S_j, 2, useQML_all); optimize_init_argument(S_j, 3, dis_mat)
            if (hasTime) {
                t_vec_select = select(colshape(time_mat,1), colshape(time_mat:>0,1))
                t_med = rows(t_vec_select)>0 ? mm_median(t_vec_select) : 1
                if(missing(t_med) || t_med<=0) t_med = 1
                tau_0 = log(max((1e-8, t_med/(-ln(0.30)))))
                optimize_init_evaluator(S_j,  &QMLE_new() ); optimize_init_params(S_j, (sig2_0, rho_0, tau_0))
                optimize_init_argument(S_j, 4, time_mat); optimize_init_argument(S_j, 5, current_resid)
            } else {
                optimize_init_evaluator(S_j,  &QMLE_bin() ); optimize_init_params(S_j, (sig2_0, rho_0))
                optimize_init_argument(S_j, 4, current_resid)
            }
            alphaHat_j = optimize(S_j)

            if (hasTime) { 
				tempSigma = Sigma_func_DGP(alphaHat_j, dis_mat, time_mat) 
				} 
            else { 
                tempSigma = Sigma_func_DGP_bin(alphaHat_j, dis_mat)
                tempSigma = tempSigma + I(n)*max((1e-8, mean(current_resid:^2)*1e-2))
            }
            SigmaHats[j] = &tempSigma; tempCS = cholesky(tempSigma); CSHats[j] = &tempCS
        }

        // --- 3c. Build Full Covariance Matrix (identical to imivreg/crsivreg) ---
        whitened_residuals = J(n, num_resids, .) 
		for(j=1; j<=num_resids; j++){ 
			whitened_residuals[.,j] = lusolve(*CSHats[j], all_residuals[.,j]) 
			}
        rhoHat_matrix = correlation(whitened_residuals); if(any(missing(rhoHat_matrix))) rhoHat_matrix = I(num_resids)
        Sigma_sim_block = J(n*num_resids, n*num_resids, 0)
        for(j=1; j<=num_resids; j++){
            for(k=j; k<=num_resids; k++){
                start_row = (j-1)*n + 1; end_row = j*n; start_col = (k-1)*n + 1; end_col = k*n
                block_jk = rhoHat_matrix[j,k] * ((*CSHats[j])'*(*CSHats[k]))
                Sigma_sim_block[start_row..end_row, start_col..end_col] = block_jk
                if(j!=k) Sigma_sim_block[start_col..end_col, start_row..end_row] = block_jk'
            }
        }
        eigval = symeigensystem(Sigma_sim_block, ., .); min_eigval = min(eigval)
        if (missing(min_eigval) || min_eigval <= 1e-8) { 
            ridge = max((1e-6, (min_eigval <= 1e-8 ? -min_eigval + 1e-6 : 1e-6)))
            Sigma_sim_block = Sigma_sim_block + I(rows(Sigma_sim_block))*ridge
        }

        // --- 3d. Simulation Loop for CCE ---
        sigLevel = .05; Bboot = 500;
        CSHat_sim_block = cholesky(Sigma_sim_block)
        SimulatedErrorMatrix = CSHat_sim_block' * rnormal(n*num_resids, Bboot, 0, 1)
        resultsMat = J(p_s, 7, .); GstarVec = J(p_s, 1, .)

        for(iCov = 1; iCov <= p_s; iCov++ ){
            beta_s_H0 = beta_s_2sls; beta_s_H0[iCov] = 0
            simPowerVec = J(l_G,1,0); pValSim = J(Bboot,l_G,1); sigLevelAdjVec = J(l_G,1,sigLevel)

            for(kk = 1; kk <= l_G; kk++){
                clustering = clusteringSet[.,kk]; G = G_vec[kk]
                
                // ================== CCE MODIFICATION START ==================
                abootVec = J(Bboot,1,.); sbootVec = J(Bboot,1,.)

                for(rr = 1; rr <= Bboot; rr++){
                    U_boot_iter = SimulatedErrorMatrix[1..n, rr]
                    V_boot_iter_mat = J(n, m_e, 0)
                    for (v_idx=1; v_idx<=m_e; v_idx++){
                         V_boot_iter_mat[.,v_idx] = SimulatedErrorMatrix[((v_idx)*n+1)..((v_idx+1)*n), rr]
                    }
                    X_e_boot = X_fs_all * pi_fs_ols + V_boot_iter_mat
                    X_s_all_boot = (X_e_boot, X_k, J(n,1,1))
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
                        t_boot_adj = t_boot * sqrt((G-1)/G) // Adjustment for CCE t-stat
                        pValSim[rr,kk] = 2*t(G-1, -abs(t_boot_adj))
                    }
                }
                // ==================  CCE MODIFICATION END  ==================
                
                valid_pvals = select(pValSim[.,kk], pValSim[.,kk] :!= .);
                if (G > 1 && rows(valid_pvals) > 0) {
                    sigLevelAdjVec[kk] = min((mm_quantile(valid_pvals,1,0.05),0.05))
                    if (missing(sigLevelAdjVec[kk])) sigLevelAdjVec[kk] = 0.05
                } else { 
					sigLevelAdjVec[kk] = 0.05 
					}

                alternatives = (range(-10,-1,1) \ range(1,10,1)) :/ sqrt(n); nalt = rows(alternatives)
                valid_power_idx = selectindex((abootVec:!=.) :& (sbootVec:!=.) :& (sbootVec :> 1e-9))
                if (rows(valid_power_idx) > 0 && G > 1) {
                     aboot_valid = abootVec[valid_power_idx, .]; sboot_valid = sbootVec[valid_power_idx, .]
                     power_sum_for_G = 0
                     for (alt_idx = 1; alt_idx <= nalt; alt_idx++) {
                         t_stats_alt = (aboot_valid :- alternatives[alt_idx]) :/ sboot_valid
                         t_stats_alt_adj = t_stats_alt * sqrt((G-1)/G)
                         pvals_for_alt = 2*t(G-1, -abs(t_stats_alt_adj))
                         power_sum_for_G = power_sum_for_G + mean(pvals_for_alt :< sigLevelAdjVec[kk])
                     }
                     simPowerVec[kk] = power_sum_for_G / nalt
                } else { 
					simPowerVec[kk] = 0 
					}
            } 

            // Select Gstar
            indStar = windStar = .;
            if (l_G > 0) { 
				maxindex(simPowerVec,1,indStar,windStar); if (rows(indStar) > 1) indStar = indStar[1]; if (missing(indStar)) indStar = l_G 
				} else { 
				indStar = . 
				}
            Gstar = missing(indStar) ? G_max : G_vec[indStar]
            if (missing(Gstar) || Gstar < 2) Gstar = G_max
            clusteringStar = missing(indStar) ? clusteringSet[., l_G] : clusteringSet[., indStar]
            GstarVec[iCov] = Gstar

            // --- 3e. Final CCE Estimation and Inference ---
            Coef_full = beta_s_2sls[iCov]
            resid_full = Y_s - X_s_all * beta_s_2sls
            se_final_vec = cluster_se_iv(X_s_all, Z_s_all, resid_full, clusteringStar)
            SE_final = se_final_vec[iCov]

            tVal=.; pStar=.; pValAdj=.; CI_lower=.; CI_upper=.;
            if (!missing(SE_final) && SE_final > 1e-9) {
                tVal = Coef_full/SE_final
                tVal_adj = tVal * sqrt((Gstar-1)/Gstar)
                pStar = 2*t(Gstar-1, -abs(tVal_adj))
                currentSigLevelAdj = (l_G > 0 && !missing(indStar) && indStar <= rows(sigLevelAdjVec)) ? sigLevelAdjVec[indStar] : sigLevel
                if(missing(currentSigLevelAdj)) currentSigLevelAdj = sigLevel

                valid_pValSim_star = select(pValSim[.,indStar], pValSim[.,indStar] :!= .)
                if(rows(valid_pValSim_star)>0) { 
					pValAdj = mean(valid_pValSim_star :<= pStar) 
					} else { 
					pValAdj = pStar 
					}
                if(missing(pValAdj)) pValAdj = pStar

                gap = -invt(Gstar-1, currentSigLevelAdj/2) * SE_final * sqrt(Gstar/(Gstar-1))
                CI_lower = Coef_full - gap
                CI_upper = Coef_full + gap
            }

            resultsMat[iCov,.] = (Coef_full, SE_final, tVal, pValAdj, CI_lower, CI_upper, Gstar)
        } 
        st_numscalar("n_obs", n)
		st_matrix("resultsMat",resultsMat)
        st_matrix("beta_s_2sls", beta_s_2sls)
		st_matrix("X_s_all_mata", X_s_all)
        st_matrix("Y_s_mata", Y_s)
		st_numscalar("p_s_mata", p_s)
    }

    // --- 4. Final Stata Output (identical to others) ---
    mata {
        n = st_numscalar("n_obs"); p_s = st_numscalar("p_s_mata")
        betaHat_final_2sls = st_matrix("beta_s_2sls")
        X_s_all = st_matrix("X_s_all_mata"); Y_s = st_matrix("Y_s_mata")
        Yhat_final = X_s_all * betaHat_final_2sls; resid_final  = Y_s - Yhat_final
        RSS_final = sum((resid_final:^2)); TSS_final = sum((Y_s :- mean(Y_s)):^2)
        R2_final = (TSS_final > 1e-9) ? (1 - RSS_final/TSS_final) : .
        RootMSE_final= ((n-p_s) > 0) ? sqrt(RSS_final/(n-p_s)) : .
        st_numscalar("R2", R2_final); st_numscalar("RootMSE", RootMSE_final)
    }

    mat colnames resultsMat = Coefficient Std_err t PValue CI_lower CI_upper Gstar
    local rnames `all_struc_regressors' _cons
    mat rownames resultsMat = `rnames'
    mata{
		st_local("R2_val", strofreal(R2_final))
		st_local("RootMSE_val", strofreal(RootMSE_final))
	}

    local obs_text "Number of obs  = "
	local obs_num = string(n_obs, "%9.0g")
    local r2_text  "R^2 (2SLS)     = "
	local r2_num = string(`R2_val', "%9.4f")
    local rmse_text "Root MSE (2SLS)= "
	local rmse_num = string(`RootMSE_val', "%9.4f")

    di _n as text "CCE method (IV) with learned clusters"
    di _col(65) as text "`obs_text'" _col(5) as result "`obs_num'"
    di _col(65) as text "`r2_text'" _col(5) as result "`r2_num'"
    di _col(65) as text "`rmse_text'" _col(5) as result "`rmse_num'"
    di as text "{hline 85}"
    di as text %15s abbrev("`depvar'",15) _col(17) " {c |} Coefficient  Std. err.      t    P>|t|     [95% conf. interval]   Clusters"
    di as text "{hline 85}"
    local list_of_struc_param_names : rownames resultsMat
	local p_s_count : word count `list_of_struc_param_names'
    forvalues i = 1/`p_s_count' {
        local name : word `i' of `list_of_struc_param_names'
        local coef = resultsMat[`i',1]
		local se = resultsMat[`i',2]
		local t = resultsMat[`i',3]
        local p = resultsMat[`i',4]
		local lci = resultsMat[`i',5]
		local uci = resultsMat[`i',6]
        local gval = resultsMat[`i',7]
        di as text %15s abbrev("`name'",15) ///
            _col(17) " {c |} " as res %8.6f `coef'  ///
            as res _col(31) %8.6f `se'  ///
            as res _col(41) %8.3f `t'    ///
            as res _col(51) %8.6f `p'    ///
            as res _col(61) %8.6f `lci' "   " %8.6f `uci' ///
            as res _col(81) %5.0f `gval'
    }
    di as text "{hline 85}"
    capture drop id0; capture drop group*; capture drop `__timedummy'
end


mata:
function cluster_se_iv(X, Z, e, group) {
    real scalar k, n, G, n_groups
    n = rows(e)
    k = cols(X)
    
    uniq_groups = uniqrows(group)
    n_groups = rows(uniq_groups)
    G = max(group)
    if (missing(G) || G==0) return(J(k,1,.))

    // --- 1. Bread: (X'PzX)^-1 ---
    real matrix ZtZ, ZtZ_inv, Pz, bread, bread_inv
    
    ZtZ = Z'*Z
    if(cond(ZtZ) > 1e12) return(J(k,1,.))
    ZtZ_inv = invsym(ZtZ)
    if (isscalar(ZtZ_inv) && missing(ZtZ_inv)) return(J(k,1,.))

    Pz = Z*ZtZ_inv*Z'
    bread_inv = X'*Pz*X
    if(cond(bread_inv) > 1e12) return(J(k,1,.))
    bread = invsym(bread_inv)

    // --- 2. Meat: Sum over clusters of squared scores ---

    real matrix meat_Z, S_g, Z_g, e_g
    meat_Z = J(cols(Z), cols(Z), 0)
    
    for (ii=1; ii<=G; ii++) {
        fii = selectindex(group:==ii)
        if (rows(fii)==0) continue
        
        Z_g = Z[fii,.]
        e_g = e[fii,.]
        
        // Calculate score for this cluster
        S_g = Z_g' * e_g
        
        // Add its contribution to the meat matrix
        meat_Z = meat_Z + S_g * S_g'
    }
    
    // --- 3. Sandwich and small sample adjustment ---
    real matrix meat_params, vcluster_mat
    
    // Build the full meat for the parameters
    meat_params = (X'*Z*ZtZ_inv) * meat_Z * (ZtZ_inv*Z'*X)
    
    // Form the sandwich
    vcluster_mat = bread * meat_params * bread
    
    // Apply small sample adjustment
    if (n_groups > 1) {
        vcluster_mat = vcluster_mat * (n_groups/(n_groups-1)) * ((n-1)/(n-k))
    }
    
    real colvector se
    se = sqrt(diagonal(vcluster_mat))
    return(se)
}

// Robust Fama-Macbeth for IV with multiple regressors
void FamaMacbethIV(Y_s_arg, X_s_all_arg, Z_s_all_arg, clustering_labels, b_coeffs_clusterwise_out) {
    real scalar G_fm, p_s_fm 
    G_fm = max(clustering_labels)
    if (missing(G_fm) || G_fm < 1) {
        G_fm = 0 
    }
    p_s_fm = cols(X_s_all_arg)
    
    b_coeffs_clusterwise_out = J(G_fm, p_s_fm, .) // Initialize with missings

    for (ii=1; ii<=G_fm; ii++ ){
        fii = selectindex(clustering_labels :== ii)
        if (rows(fii) < p_s_fm) { // Not enough obs in cluster
            continue
        }

        Y_c = Y_s_arg[fii,.]
        X_s_c = X_s_all_arg[fii,.]
        Z_s_c = Z_s_all_arg[fii,.]
        
        // --- Robust 2SLS calculation for the cluster ---
        real matrix ZZ, ZZ_inv, XZ, Xprime_PZ_X, Xprime_PZ_X_inv, beta_c
        
        // Check rank of instruments in cluster
        if(rank(Z_s_c) < cols(Z_s_c)) continue;
        ZZ = Z_s_c'*Z_s_c
        if(cond(ZZ) > 1e12) continue; // Check condition number for stability
        ZZ_inv = invsym(ZZ)
        if (isscalar(ZZ_inv) && missing(ZZ_inv)) continue
        
        XZ = X_s_c'*Z_s_c
        Xprime_PZ_X = XZ * ZZ_inv * XZ'
        
        // Check rank for structural regressors projection
        if(rank(Xprime_PZ_X) < cols(Xprime_PZ_X)) continue;
        if(cond(Xprime_PZ_X) > 1e12) continue;
        Xprime_PZ_X_inv = invsym(Xprime_PZ_X)
        if (isscalar(Xprime_PZ_X_inv) && missing(Xprime_PZ_X_inv)) continue
        
        beta_c = Xprime_PZ_X_inv * (XZ * ZZ_inv * (Z_s_c'*Y_c))
        
        if (!any(missing(beta_c))) {
            b_coeffs_clusterwise_out[ii,.] = beta_c'
        }
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