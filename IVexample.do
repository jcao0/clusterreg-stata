*=========================================================================*
*  Simulated panel data for learned-cluster IV — SAR error specification  *
*=========================================================================*
clear all
set more off
set seed 12345

*----------------------------Parameters-----------------------------------*
local N_cs     100          // 
local T_time   3            // timeperiod
local rho_sar  0.30         // SAR parameter ρ  (|ρ| < 1/λ_max(W))
local d_cut    1.5          //
local rho_uv   0.60         // corr(U , V)  ⇒ endog
local sigma_u  1            //
local sigma_v  1            //
*-------------------------------------------------------------------------*

local total_obs = `N_cs'*`T_time'
set obs `total_obs'

gen long id_cs       = .
gen byte time_period = .

forvalues i = 1/`N_cs' {
    forvalues t = 1/`T_time' {
        local row = (`i'-1)*`T_time'+`t'
        quietly replace id_cs       = `i' in `row'
        quietly replace time_period = `t' in `row'
    }
}

gen coord1 = (ceil(id_cs/10) - .5 + runiform()-.5)
gen coord2 = (mod(id_cs-1,10) + .5 + runiform()-.5)
su coord1 , meanonly
replace coord1 = 10*(coord1-r(min))/(r(max)-r(min))
su coord2 , meanonly
replace coord2 = 10*(coord2-r(min))/(r(max)-r(min))

bysort id_cs: gen X_k = rnormal(1,2)  if _n==1
bysort id_cs: replace X_k = X_k[1]
bysort id_cs: gen Z_i = rnormal(0,1.5) if _n==1
bysort id_cs: replace Z_i = Z_i[1]

mata:
    Ns = strtoreal(st_local("N_cs"))
    Tspan = strtoreal(st_local("T_time"))
    rho = strtoreal(st_local("rho_sar"))
    du = strtoreal(st_local("sigma_u"))
    dv = strtoreal(st_local("sigma_v"))
    dcut = strtoreal(st_local("d_cut"))
    rhoUV = strtoreal(st_local("rho_uv"))

    st_view(X1=., .,"coord1")
    st_view(X2=., .,"coord2")
    st_view(T  =., .,"time_period")

    D = J(Ns, Ns, .)
    for(i=1; i<=Ns; i++) {
        for(j=1; j<=Ns; j++) {
            D[i,j] = sqrt((X1[i]-X1[j])^2 + (X2[i]-X2[j])^2)
        }
    }
    W_small = (D:<dcut) :!= I(Ns)
    rowsum  = rowsum(W_small)
    W_small = diag(1:/rowsum) * W_small

    I_T   = I(Tspan)
	W_big = J(Ns*Tspan, Ns*Tspan, 0)
	for (tt = 1; tt <= Tspan; tt++) {
		W_big[(tt-1)*Ns+1 .. tt*Ns, (tt-1)*Ns+1 .. tt*Ns] = W_small
	}

    Ntot = Ns*Tspan

    epsV  = dv * rnormal(Ntot,1,0,1)
    Verr  = invsym(I(Ntot) - rho*W_big) * epsV

    epsU  = du * rnormal(Ntot,1,0,1)
    Ustar = invsym(I(Ntot) - rho*W_big) * epsU
    Uerr  = rhoUV*Verr + sqrt(1-rhoUV^2)*Ustar

    (void) st_addvar("double", "V_error")
    (void) st_addvar("double", "U_error")
    st_store(., "V_error", Verr)
    st_store(., "U_error", Uerr)
end

gen X_e = 0.5 + 1.2*Z_i + 0.8*X_k + V_error
gen Y   = 1.0 + 2.0*X_e - 0.5*X_k + U_error

imivreg_test Y X_k ///
       , iv(X_e = Z_i) ///
       cluster(coord1 coord2) ///
       timeperiod(time_period)

crsivreg_test Y X_k ///
       , iv(X_e = Z_i) ///
       cluster(coord1 coord2) ///
       timeperiod(time_period)

cceivreg Y X_k ///
       , iv(X_e = Z_i) ///
       cluster(coord1 coord2) ///
       timeperiod(time_period)