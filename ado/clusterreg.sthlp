{smcl}
{* *! version May2025}{...}

{title:Title}

{p2colset 5 22 15 2}{...}
{p2col :{cmd:clusterreg} {hline 2}}Data-driven Inference for Dependent Data with Learned Clusters{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:clusterreg}
{depvar} {indepvars}
{ifin}
{it:{weight}}
{cmd:,}
{opt coor(varlist)}
[{opt time(varname)}]
[{opt type(string)}]
{p_end}

{marker options}{...}
{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt :{opt coor(varlist)}}Coordinate variables for spatial or feature-based clustering (e.g., latitude longitude){p_end}

{syntab:Optional}
{synopt :{opt time(varname)}}Optional time dimension variable. If provided, estimates spatio-temporal dependence; otherwise, estimates spatial/feature dependence only.{p_end}
{synopt :{opt type(string)}}Inference method: {it:IM} (Ibragimov-Müller, default), {it:CCE} (Bester-Conley-Hansen Cluster Covariance Estimator), or {it:CRS} (Canay-Romano-Shaikh randomization){p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}{cmd:clusterreg} implements data-driven cluster-robust inference for dependent data (spatial, temporal, feature-based, or network) using learned clusters. It provides three inference methods:{p_end}

{p 8 12}{bf:1. Data-Driven Clustering}: Uses the Partitioning Around Medoids (PAM) algorithm with Euclidean distance on the variables provided in {opt coor()} to automatically group observations into clusters. The number of clusters is chosen adaptively.{p_end}

{p 8 12}{bf:2. Optimal Cluster Selection}: Selects the number of clusters (`G*`) by simulating data under the null (using a QMLE-estimated spatial or spatio-temporal covariance structure) and choosing the `G` that maximizes the simulated power of the chosen test ({opt type()}) under local alternatives.{p_end}

{p 8 12}{bf:3. Robust Inference Methods}: Provides three approaches valid with data-driven clusters:{p_end}
{p 12 15}{it:IM} (Default): Ibragimov and Muller (2010) t-test using cluster-level statistics.{p_end}
{p 12 15}{it:CRS}: Canay, Romano, and Shaikh (2017) randomization test (permutation/sign-flipping).{p_end}
{p 12 15}{it:CCE}: Bester, Conley, and Hansen (2011) cluster-robust covariance estimator.{p_end}

{pstd}The method is particularly valuable when:{p_end}
{p 8 12}- Natural clustering boundaries are unknown or subjective.{p_end}
{p 8 12}- The dependence structure is complex (e.g., spatial and temporal).{p_end}
{p 8 12}- Traditional clustering approaches might yield misspecified groups.{p_end}

{marker theoretical}{...}
{title:Theoretical Foundations (based on Cao et al., 2024)}

{pstd}The procedure's validity relies on conditions including:{p_end}
{pmore}1. {bf:Metric Space Properties}: Assumptions on the space defined by {opt coor()} variables (e.g., Ahlfors regularity).{p_end}
{pmore}2. {bf:Mixing Condition}: Dependence between observations decays sufficiently fast with distance (in coordinate/time space).{p_end}
{pmore}3. {bf:Cluster Properties}: The clustering algorithm produces reasonably balanced clusters with vanishing inter-cluster dependence asymptotically.{p_end}

{pstd}Key theoretical results allow for:{p_end}
{pmore}- Asymptotic normality of cluster-level estimators (for IM).{p_end}
{pmore}- Asymptotic validity of the randomization procedure (for CRS).{p_end}
{pmore}- Valid Type I error control under the specified conditions, accounting for data-driven cluster selection.{p_end}

{marker examples}{...}
{title:Examples}

{pstd}Basic IM regression (spatial dependence only):{p_end}
{phang2}{cmd:. clusterreg crime income pop, coor(lat lon)}{p_end}

{pstd}Basic IM regression (spatio-temporal dependence):{p_end}
{phang2}{cmd:. clusterreg crime income pop, coor(lat lon) time(year)}{p_end}

{pstd}CRS method with weights (spatial only):{p_end}
{phang2}{cmd:. clusterreg pollution temp wind [aweight=area], coor(x y) type(CRS)}{p_end}

{pstd}Subset analysis with CCE (spatio-temporal):{p_end}
{phang2}{cmd:. clusterreg growth gdp trade if year>2000, coor(lon lat) time(quarter) type(CCE)}{p_end}

{pstd}Florida homicide analysis (Example 2 from README):{p_end}
{phang2}{cmd:. clear all}{p_end}
{phang2}{cmd:. set more off}{p_end}
{phang2}{cmd:. set seed 100}{p_end}
{phang2}{cmd:. webuse homicide_1960_1990}{p_end}
{phang2}{cmd:. keep if sname == "Florida"}{p_end}
{phang2}{cmd:/* Spatio-temporal using IM */}{p_end}
{phang2}{cmd:. clusterreg hrate divorce unemployment ln_income poverty, coord(_CX _CY) time(year) type("IM")}{p_end}
{phang2}{cmd:/* Spatial only using CRS */}{p_end}
{phang2}{cmd:. clusterreg hrate divorce unemployment ln_income poverty, coord(_CX _CY) type("CRS")}{p_end}



{marker saved_results}{...}
{title:Saved results}

{pstd}{cmd:clusterreg} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}Number of observations{p_end}
{synopt:{cmd:e(k)}}Optimal number of clusters selected (`G*`){p_end}
{* Add other scalars if relevant, e.g., e(rank), e(r2), e(rmse) - check the code output *}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(rmse)}}Root Mean Squared Error{p_end}
{synopt:{cmd:e(N_clust)}}Optimal number of clusters (same as e(k)){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:clusterreg}{p_end}
{synopt:{cmd:e(depvar)}}Name of dependent variable{p_end}
{synopt:{cmd:e(indeps)}}Names of independent variables{p_end}
{synopt:{cmd:e(clustvar)}}Name(s) of coordinate variables{p_end}
{synopt:{cmd:e(timevar)}}Name of time variable (if specified){p_end}
{synopt:{cmd:e(method)}}Inference method used (IM, CRS, or CCE){p_end}
{synopt:{cmd:e(predict)}}Program used for predict{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}Coefficient vector{p_end}
{synopt:{cmd:e(V)}}Variance-covariance matrix (robust to selected clustering){p_end}
{* Add other matrices if relevant - check the code output *}

{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}Marks estimation sample{p_end}
{p2colreset}{...}

{marker references}{...}
{title:References}

{pstd}{bf:Core methodology}:{p_end}
{pmore}Cao, J., Hansen, C., Kozbur, D., & Villacorta, L. (2024). Inference for dependent data with learned clusters. {it:Review of Economics and Statistics}, 1-45.{p_end}

{pstd}Underlying inference procedures adapted by the package:{p_end}

{pstd}{bf:IM method}:{p_end}
{pmore}Ibragimov, R., & Müller, U. K. (2010). t-Statistic based correlation and heterogeneity robust inference. {it:Journal of Business & Economic Statistics}, 28(4), 453-468.{p_end}

{pstd}{bf:CRS method}:{p_end}
{pmore}Canay, I. A., Romano, J. P., & Shaikh, A. M. (2017). Randomization tests under an approximate symmetry assumption. {it:Econometrica}, 85(3), 1013-1030.{p_end}

{pstd}{bf:CCE method}:{p_end}
{pmore}Bester, C. A., Conley, T. G., & Hansen, C. B. (2011). Inference with dependent data using cluster covariance estimators. {it:Journal of Econometrics}, 165(2), 137-151.{p_end}
