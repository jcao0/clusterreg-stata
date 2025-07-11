--------------------------------------------------------------------------------
README for the clusterreg and clusterivreg Stata Package
Version: July2025
--------------------------------------------------------------------------------

Authors: Jianfei Cao, Chris Hansen, Damian Kozbur, Lucciano Villacorta, Romano Li, Jiacheng Liang

Correspondence: Jianfei Cao, Northeastern University

Associated Paper: Cao, J., Hansen, C., Kozbur, D., & Villacorta, L. (2024). Inference for dependent data with learned clusters. Review of Economics and Statistics, 1-45.


--------------------------------------------------------------------------------
1. Purpose
--------------------------------------------------------------------------------

This package provides Stata commands to perform inference in linear regression models with spatially or feature-dependent data for both standard OLS and Instrumental Variables settings. Standard errors and inference procedures are adjusted using data-driven ("learned") clusters based on coordinates or other features provided by the user. The number of clusters is selected adaptively to balance the size-power trade-off. 

The package implements three main approaches accessible via the primary `clusterreg` and `clusterivreg` commands:
* IM: The t-statistic based procedure of Ibragimov and Müller (2010).
* CRS: The randomization inference procedure of Canay, Romano, and Shaikh (2017).
* CCE: The full-sample Cluster Covariance Estimator based on Bester, Conley, and Hansen (2011).

--------------------------------------------------------------------------------
2. Installation and Dependencies
--------------------------------------------------------------------------------

To use the `clusterreg` package locally, you can copy all files in the "/ado" folder into your personal Stata ado directory.

**Files to Copy:**
* Dispatchers: `clusterreg.ado`, `clusterivreg.ado`
* OLS Methods: `imreg.ado`, `crsreg.ado`, `ccereg.ado`
* IV Methods: `imivreg.ado`, `crsivreg.ado`, `cceivreg.ado`
* Core Helper: `clpam.ado`
* Help Files: `clusterreg.sthlp`, `clusterivreg.sthlp` (recommended)

**Finding Your Personal Ado Directory:**

* **Windows:** Typically `c:\ado\personal`, but it might be different. In Stata, type `personal` to see the path.
* **Mac:** Typically `~/Documents/Stata/ado/personal` or `~/Library/Application Support/Stata/ado/personal`. In Stata, type `personal` to see the path.
* **Unix/Linux:** Typically `~/ado/personal`. In Stata, type `personal` to see the path.

For more details on personal ado directories, see the Stata FAQ:
https://www.stata.com/support/faqs/programming/personal-ado-directory/

After copying the files, you should be able to use the `clusterreg` and `clusterivreg` commands in Stata. You might need to type `discard` or restart Stata for it to recognize the new commands.

This package requires the `moremata` package for some Mata functions. You can install it or ensure it's up-to-date by typing the following command in Stata:

    ssc install moremata, replace

The package also relies heavily on embedded Mata functions for calculations. Ensure your Stata installation includes Mata.

--------------------------------------------------------------------------------
3. Usage and Syntax
--------------------------------------------------------------------------------
**Syntax**

`clusterreg` (for OLS)

    clusterreg depvar indepvars [if] [in] [weight], cluster(varlist) [time(varname) type(string)]

`clusterivreg` (for IV)



    clusterivreg depvar [exog_vars] (endog_vars = inst_vars) [if] [in] [weight], cluster(varlist) [time(varname) type(string)]

The IV syntax is standard, with endogenous variables and instruments specified within parentheses. It fully supports multiple endogenous variables and multiple instruments.

**Options**
* cluster(varlist): (Required) Specifies the numeric variables representing the coordinates or features used for forming clusters (e.g., latitude, longitude).

* time(varname): (Optional) Specifies a numeric variable for the time period. If provided, the model can estimate spatio-temporal dependence. If omitted, only spatial/feature dependence is modeled.

* type(string): (Optional) Specifies the inference method. Default is "IM".
	* "IM": Ibragimov and Müller (2010) method.
	* "CRS": Canay, Romano, and Shaikh (2017) method.
	* "CCE": Bester, Conley, and Hansen (2011) cluster covariance estimator.

**Method Overview**

Regardless of the command or type, the procedure follows these steps:

1. Clustering: A set of potential data partitions is created using the k-medoids algorithm on the variables specified in `cluster()`. The number of clusters, `G`, is varied from 2 up to a maximum of $\lceil n^{1/3}\rceil$.

2. Covariance Estimation: The procedure estimates the parameters of a spatial/spatio-temporal covariance function. This is done via Quasi-Maximum Likelihood Estimation on the residuals from the initial full-sample OLS or 2SLS estimation.

3. Optimal Cluster Selection: The package simulates data under the null hypothesis using the estimated covariance structure. It then calculates the statistical power for each potential number of clusters (`G`) and selects the `G*` that yields the highest simulated power while maintaining correct test size.

4. Final Inference: The final p-values and confidence intervals are reported based on the chosen inference method using the optimal number of clusters, `G*`.

--------------------------------------------------------------------------------
4. Examples
--------------------------------------------------------------------------------

The package includes two example do-files (`example.do` and `IVexample.do`) that replicate the analyses from the associated paper and demonstrate the commands in practice.