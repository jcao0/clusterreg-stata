*========================================*
*  Example 1: Single Exogenous Variable  *
*========================================*
clear all
set more off
set seed 12345

*----------------------------Parameters-----------------------------------*
local N_cs     100
local T_time   3
local rho_sar  0.30
local d_cut    1.5
local rho_uv   0.60
local sigma_u  1
local sigma_v  1
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

clusterivreg Y X_k (X_e = Z_i), coord(coord1 coord2) time(time_period) type("CRS")



*===========================================*
*  Example 2: Multiple Exogenous Variables  *
*===========================================*
clear all
set more off
set seed 12345

*----------------------------Parameters-----------------------------------*
local N_cluster   5
local N_per_clust 50
local sigma_u     1.2
local sigma_v     1
local rho_uv      0.65
*-------------------------------------------------------------------------*

set obs `=`N_cluster'*`N_per_clust''
gen long id_cs      = _n
gen byte time_period = 2 

matrix C = ( 0,0 \ 10,0 \ 0,10 \ 10,10 \ 5,5 )

gen double coord1 = .
gen double coord2 = .
forvalues g = 1/`N_cluster' {
    local start = (`g'-1)*`N_per_clust' + 1
    local end   = `g'*`N_per_clust'
    replace coord1 = C[`g',1] + runiform(-1.8,1.8) in `start'/`end'
    replace coord2 = C[`g',2] + runiform(-1.8,1.8) in `start'/`end'
}

egen true_cluster = group(coord1 coord2), label

gen w1 = rnormal(0,1)
gen w2 = rnormal(0,1)

gen z1 = rnormal(0,1)
gen z2 = rnormal(0,1)

tempvar u_raw v1_raw v2_raw
gen `u_raw'  = rnormal(0,`sigma_u')
gen `v1_raw' = rnormal(0,`sigma_v')
gen `v2_raw' = rnormal(0,`sigma_v')

egen u_cluster_eff  = mean(`u_raw' ), by(true_cluster)
egen v1_cluster_eff = mean(`v1_raw'), by(true_cluster)
egen v2_cluster_eff = mean(`v2_raw'), by(true_cluster)

gen u  = sqrt(1-`rho_uv'^2)*`u_raw' + `rho_uv'*v1_cluster_eff
gen v1 = v1_cluster_eff
gen v2 = v2_cluster_eff

gen x1 =  1   + 1.5*z1 + 0.5*z2 + 0.8*w1 + v1
gen x2 = -0.5 + 0.7*z1 + 1.2*z2 + 0.4*w2 + v2
gen y = 2.0*x1 - 1.5*x2 + 1.0*w1 + 0.5*w2 + u


clusterivreg y w1 w2 (x1 x2 = z1 z2), coord(coord1 coord2)