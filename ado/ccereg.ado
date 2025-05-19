program define ccereg, eclass
    syntax varlist(min=2) [if] [in], cluster(varlist) [timeperiod(varname)]

    tempname __htscalar
    if "`timeperiod'" == "" {
        tempvar __htedummy
        quietly gen double `__htedummy' = 0 `if' `in'
        local timeperiod `__htedummy'
        scalar `__htscalar' = 0
    }
    else {
        quietly summarize `timeperiod' `if' `in'
        scalar `__htscalar' = (r(min) < r(max))
    }

    tokenize `varlist'
    local yvar `1'
    macro shift
    local xvars `*'

    mata {
        uniformseed(12345)
        rseed(12345)

        Y         = st_data(., "`yvar'")
        X         = st_data(., tokens("`xvars'"))
        X         = X, J(rows(X), 1, 1)
        n         = length(Y)
        coord     = st_data(., "`cluster'")
        timePeriod= st_data(., "`timeperiod'")
        p         = cols(X)

        G_max = ceil(n^(1/3))
        st_numscalar("G_max", G_max)

        G_vec = range(2, G_max, 1)
        l_G   = length(G_vec)
        st_numscalar("l_G", l_G)
    }

    gen id0 = _n
    matrix dissim dis_mat  = `cluster', L2
    matrix dissim time_mat = `timeperiod', L2
    forvalues i = 2/`=G_max' {
        quietly clpam group`i', distmat(dis_mat) id(id0) medoids(`i') ga
    }

    mata {
        hasTime     = st_numscalar("`__htscalar'")
        data_medoids= st_data(., .)
        clusteringSet = data_medoids[., (cols(data_medoids)-G_max+2)..cols(data_medoids)]

        M      = I(n) - X * invsym(X'*X) * X'
        Qd=Rd=ex=.  // for qrdp
        qrdp(M, Qd, Rd, ex)
        useQML = ex[1..(cols(ex)-p-1)]
        resid  = Y - X*(invsym(X'*X)*X')*Y
        dis_mat = st_matrix("dis_mat")
        time_mat= st_matrix("time_mat")

        S = optimize_init()
        if (hasTime) {
            optimize_init_evaluator    (S, &QMLE_new())
            optimize_init_which        (S, "min")
            optimize_init_evaluatortype(S, "d0")
            sig2_0 = log(mean(resid:^2))
            d_vec  = select(colshape(dis_mat,1), colshape(dis_mat:>0,1))
            rho_0  = log(mm_median(d_vec)/(-ln(0.30)))
            t_vec  = select(colshape(time_mat,1), colshape(time_mat:>0,1))
            tau_0  = log((rows(t_vec)>0?mm_median(t_vec):1)/(-ln(0.30)))
            optimize_init_params        (S, (sig2_0, rho_0, tau_0))
            optimize_init_argument      (S, 1, M)
            optimize_init_argument      (S, 2, useQML)
            optimize_init_argument      (S, 3, dis_mat)
            optimize_init_argument      (S, 4, time_mat)
            optimize_init_argument      (S, 5, resid)
        }
        else {
            optimize_init_evaluator    (S, &QMLE_bin())
            optimize_init_which        (S, "min")
            optimize_init_evaluatortype(S, "d0")
            sig2_0 = log(mean(resid:^2))
            d_vec  = select(colshape(dis_mat,1), colshape(dis_mat:>0,1))
            rho_0  = log(mm_median(d_vec)/(-ln(0.30)))
            optimize_init_params        (S, (sig2_0, rho_0))
            optimize_init_argument      (S, 1, M)
            optimize_init_argument      (S, 2, useQML)
            optimize_init_argument      (S, 3, dis_mat)
            optimize_init_argument      (S, 4, resid)
        }
        alphaHat = optimize(S)

        if (hasTime) {
            SigmaHat = Sigma_func_DGP(alphaHat, dis_mat, time_mat)
        }
        else {
            SigmaHat = Sigma_func_DGP_bin(alphaHat, dis_mat)
            epsHat   = mean(resid:^2)*1e-2
            SigmaHat = SigmaHat + I(rows(SigmaHat))*epsHat
        }
    }

    mata {
        betaHat    = invsym(X'*X)*X'*Y
        resid_full = Y - X*betaHat
    }

    mata {
        sigLevel = .05
        Bboot    = 1000
        CSHat    = cholesky(SigmaHat)
        UbootMat = CSHat' * rnormal(n, Bboot, 0, 1)
        resultsMat= J(p,7,0)
    }

    mata {
        for (iCov=1; iCov<=p; iCov++) {
            Coef_full = betaHat[iCov]
            if (iCov==1)             betaNull = (0 \ betaHat[2..p])
            else if (iCov==p)        betaNull = (betaHat[1..p-1] \ 0)
            else                     betaNull = (betaHat[1..iCov-1] \ 0 \ betaHat[iCov+1..p])

            simPowerVec    = J(l_G,1,0)
            pValSim        = J(Bboot,l_G,0)
            sigLevelAdjVec = J(l_G,1,0)

            for (kk=1; kk<=l_G; kk++) {
                clustering = clusteringSet[.,kk]
                G          = G_vec[kk]
                abootVec   = J(Bboot,1,0)
                sbootVec   = J(Bboot,1,0)
                for (rr=1; rr<=Bboot; rr++) {
                    Uboot        = UbootMat[,rr]
                    Yboot        = X*betaNull + Uboot
                    bh_boot      = invsym(X'*X)*X'*Yboot
                    resid_boot   = Yboot - X*bh_boot
                    se_boot      = cluster_se(X, resid_boot, invsym(X'*X), clustering)

                    abootVec[rr] = bh_boot[iCov]
                    sbootVec[rr] = se_boot[iCov]
                    t_boot       = abootVec[rr]/sbootVec[rr]*sqrt((G-1)/G)
                    pValSim[rr,kk] = 2*t(G-1, -abs(t_boot))
                }
                sigLevelAdjVec[kk] = min((mm_median(pValSim[.,kk]), 0.05))


                alternatives = (range(-10,-1,1)\range(1,10,1)) / sqrt(n)
                nalt         = rows(alternatives)
                powerMat     = J(Bboot, nalt, 0)
                for (a=1; a<=nalt; a++) {
                    powerMat[.,a] = (abootVec :- alternatives[a]) :/ sbootVec
                }
                pvals_alt    = 2*t(G-1, -abs(powerMat))
                simPowerVec[kk] = mean(mean((pvals_alt :< sigLevelAdjVec[kk]))')
            }


            indStar     = .
			windStar = .
            maxindex(simPowerVec,1,indStar,windStar)
			if (rows(indStar) > 1) indStar = indStar[1]
            Gstar       = G_vec[indStar]
            clusteringStar = clusteringSet[.,indStar]

            se_final_vec= cluster_se(X, resid_full, invsym(X'*X), clusteringStar)
            SE_final    = se_final_vec[iCov]
            tVal        = Coef_full/SE_final
            tVal_adj    = tVal*sqrt((Gstar-1)/Gstar)
            pStar       = 2*t(Gstar-1, -abs(tVal_adj))
            pValAdj     = mean(pStar :>= pValSim[.,indStar])
            gap         = -invt(Gstar-1, sigLevelAdjVec[indStar]/2)*SE_final
            CI_lower    = Coef_full - gap
            CI_upper    = Coef_full + gap

            resultsMat[iCov,] = (Coef_full, SE_final, tVal, pValAdj, CI_lower, CI_upper, Gstar)
        }
        st_numscalar("n_obs", n)
        st_matrix("resultsMat", resultsMat)
    }

    mata {
        R2_final     = 1 - sum((resid_full:^2)) / sum((Y :- mean(Y)):^2)
        RootMSE_final= sqrt(sum((resid_full:^2)) / (n - cols(X)))
        st_numscalar("R2", R2_final)
        st_numscalar("RootMSE", RootMSE_final)
    }
    mat colnames resultsMat = Coefficient Std_err t PValue CI_lower CI_upper Gstar
    mat rownames resultsMat = `xvars' _cons
    mata: st_local("R2_val", strofreal(R2_final))
    mata: st_local("RootMSE_val", strofreal(RootMSE_final))

    local obs_text "Number of obs  = "
    local obs_num  = string(n_obs,"%9.0g")
    local r2_text  "R^2            = "
    local r2_num   = string(`R2_val',"%9.4f")
    local rmse_text"Root MSE       = "
    local rmse_num = string(`RootMSE_val',"%9.4f")

    di _n as text "CCE method with learned cluster"
    di _col(65) as text "`obs_text'"  _col(5) as result "`obs_num'"
    di _col(65) as text "`r2_text'"   _col(5) as result "`r2_num'"
    di _col(65) as text "`rmse_text'" _col(5) as result "`rmse_num'"
    di as text "{hline 85}"
    di as text %12s abbrev("`yvar'",12) _col(14) " {c |} Coefficient  Std. err.      t    P>|t|     [95% conf. interval]   Clusters"
    di as text "{hline 85}"
    forvalues i = 1/`=rowsof(resultsMat)' {
        local name : word `i' of `xvars' _cons
        local coef = resultsMat[`i',1]
        local se   = resultsMat[`i',2]
        local t    = resultsMat[`i',3]
        local p    = resultsMat[`i',4]
        local lci  = resultsMat[`i',5]
        local uci  = resultsMat[`i',6]
        local gval = resultsMat[`i',7]
        di as text %12s abbrev("`name'",12) ///
           _col(14) " {c |} " as res %8.6f `coef'  ///
           as res _col(28) %8.6f `se'   ///
           as res _col(38) %8.3f `t'    ///
           as res _col(48) %8.6f `p'    ///
           as res _col(58) %8.6f `lci' "   " %8.6f `uci' ///
           as res _col(78) %5.0f `gval'
    }
    di as text "{hline 85}"
end



// -----------------------------------------------------------------------------
// Mata Helper Functions
// -----------------------------------------------------------------------------

mata:
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
end
