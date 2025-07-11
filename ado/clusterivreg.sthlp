{smcl}
{* *! version July 2025}{...}

{title:Title}

{p2colset 5 24 17 2}{...}
{p2col :{cmd:clusterivreg} {hline 2}}Data-driven IV Inference for Dependent Data with Learned Clusters{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2}
{cmd:clusterivreg} {depvar} [{it:exog_vars}] ({it:endog_vars} {cmd:=} {it:inst_vars}) {ifin} {it:{weight}} {cmd:,} {opt coord(varlist)} [{opt time(varname)}] [{opt type(string)}]
{p_end}

{pstd}
where {it:endog_vars} are the endogenous regressors and {it:inst_vars} are the external instrumental variables. This syntax is analogous to Stata's official {cmd:ivregress} command.
{p_end}

{marker options}{...}
{title:Options}

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt :{opt coord(varlist)}}Specifies the numeric variables representing the coordinates or features used for calculating distances and forming clusters (e.g., latitude, longitude). At least one variable must be provided.{p_end}

{syntab:Optional}
{synopt :{opt time(varname)}}Specifies a numeric variable indicating the time period for each observation. If provided, the estimation of the error covariance structure can incorporate temporal decay alongside spatial/feature decay.{p_end}
{synopt :{opt type(string)}}Specifies the inference method to use. The default is {cmd:"IM"}.{p_end}
{p 12 24 2}
{it:IM}: Ibragimov and Muller (2010) method.{p_end}
{p 12 24 2}
{it:CRS}: Canay, Romano, and Shaikh (2017) randomization method.{p_end}
{p 12 24 2}
{it:CCE}: Bester, Conley, and Hansen (2011) cluster covariance estimator method.{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:clusterivreg} is the instrumental variables (IV) counterpart to {cmd:clusterreg}. It implements data-driven cluster-robust inference for linear models with endogenous regressors and dependent data (spatial, temporal, or feature-based). It provides three inference methods, all adapted for the 2SLS/IV framework:
{p_end}

{p 8 12}{bf:1. Data-Driven Clustering}: Uses the Partitioning Around Medoids (PAM) algorithm with Euclidean distance on the variables provided in {opt cluster()} to automatically group observations. The number of clusters is chosen adaptively.{p_end}

{p 8 12}{bf:2. Optimal Cluster Selection}: Selects the number of clusters (`G*`) by simulating data under the null (using a QMLE-estimated spatial or spatio-temporal covariance structure for both the structural and first-stage errors) and choosing the `G` that maximizes the simulated power of the chosen test ({opt type()}) under local alternatives.{p_end}

{p 8 12}{bf:3. Robust IV Inference Methods}: Provides three approaches for IV estimation:{p_end}
{p 12 15}{it:IM} (Default): Ibragimov and Muller (2010) t-test using cluster-level 2SLS statistics.{p_end}
{p 12 15}{it:CRS}: Canay, Romano, and Shaikh (2017) randomization test applied to cluster-level 2SLS statistics.{p_end}
{p 12 15}{it:CCE}: Bester, Conley, and Hansen (2011) cluster-robust covariance estimator applied to the full-sample 2SLS estimates.{p_end}

{pstd}This command is essential when you suspect that one or more regressors are endogenous {it:and} that the errors exhibit complex dependence structures not easily captured by standard clustering methods.{p_end}

{marker theoretical}{...}
{title:Theoretical Foundations (based on Cao et al., 2024)}

{pstd}
The procedure's validity relies on conditions including:{p_end}
{pmore}1. {bf:Metric Space Properties}: Assumptions on the space defined by {opt coord()} variables (e.g., Ahlfors regularity).{p_end}
{pmore}2. {bf:Mixing Condition}: Dependence between observations decays sufficiently fast with distance (in coordinate/time space).{p_end}
{pmore}3. {bf:IV Validity}: Standard instrument relevance and exogeneity conditions must hold.{p_end}

{pstd}
Key theoretical results allow for valid Type I error control under the specified conditions, accounting for both endogeneity and data-driven cluster selection.
{p_end}

{marker examples}{...}
{title:Examples}

{pstd}A simple IV model with one endogenous and one instrument variable:{p_end}
{phang2}{cmd:. clusterivreg wage tenure (educ = father_educ), coord(lat lon)}{p_end}

{pstd}A model with exogenous controls and spatio-temporal dependence, using the CRS method:{p_end}
{phang2}{cmd:. clusterivreg log_price sqft (crime_rate = police_funding), coord(x_coord y_coord) time(quarter) type(CRS)}{p_end}

{pstd}A model with multiple endogenous variables and multiple instruments:{p_end}
{phang2}{cmd:. clusterivreg growth gdp (invest trade = tariff infrastructure), coord(region_x region_y) type(CCE)}{p_end}

{pstd}For a full demonstration, please see the provided {bf:IVexample.do} file.{p_end}

{marker saved_results}{...}
{title:Saved results}

{pstd}
{cmd:clusterivreg} saves the following in {cmd:e()}:
{p_end}

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}Number of observations{p_end}
{synopt:{cmd:e(k)}}Optimal number of clusters selected (`G*`){p_end}
{synopt:{cmd:e(r2_2sls)}}R-squared from the full-sample 2SLS regression{p_end}
{synopt:{cmd:e(rmse_2sls)}}Root Mean Squared Error from the full-sample 2SLS regression{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:clusterivreg}{p_end}
{synopt:{cmd:e(depvar)}}Name of dependent variable{p_end}
{synopt:{cmd:e(exog_vars)}}Names of exogenous regressors{p_end}
{synopt:{cmd:e(endog_vars)}}Names of endogenous regressors{p_end}
{synopt:{cmd:e(inst_vars)}}Names of instrumental variables{p_end}
{synopt:{cmd:e(coord)}}Name(s) of coordinate variables{p_end}
{synopt:{cmd:e(timevar)}}Name of time variable (if specified){p_end}
{synopt:{cmd:e(method)}}Inference method used (IM, CRS, or CCE){p_end}
{synopt:{cmd:e(predict)}}Program used for predict{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}Coefficient vector{p_end}
{synopt:{cmd:e(V)}}Variance-covariance matrix (robust to selected clustering){p_end}

{p2col 5 22 26 2: Functions}{p_end}
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
