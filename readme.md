--------------------------------------------------------------------------------
README for the clusterreg Stata Package
Version: May2025
--------------------------------------------------------------------------------

Authors: Jianfei Cao, Chris Hansen, Damian Kozbur, Lucciano Villacorta, Romano Li, Jiacheng Liang

Correspondence: Jianfei Cao

Associated Paper: Cao, J., Hansen, C., Kozbur, D., & Villacorta, L. (2024). Inference for dependent data with learned clusters. Review of Economics and Statistics, 1-45.


--------------------------------------------------------------------------------
1. Purpose
--------------------------------------------------------------------------------

This package provides Stata commands to perform inference in linear regression models with spatially or feature-dependent data. Standard errors and inference procedures are adjusted using data-driven ("learned") clusters based on coordinates or other features provided by the user. The number of clusters is selected adaptively to balance the size-power trade-off. 

The package implements three main approaches accessible via the primary `clusterreg` command: the procedure of Ibragimov and Muller (2010) (IM), the procedure of Canay et al. (2017) (CRS), and inference based on the cluster covariance estimator as described in Bester et al. (2011) (CCE).

--------------------------------------------------------------------------------
2. Installation and Dependencies
--------------------------------------------------------------------------------

To use the `clusterreg` package locally, you can copy all files in the "/ado" folder into your personal Stata ado directory.

**Files to Copy:**
* `clusterreg.ado`
* `clusterreg.sthlp`
* `imreg.ado`
* `crsreg.ado`
* `ccereg.ado`
* `clpam.ado`

**Finding Your Personal Ado Directory:**

* **Windows:** Typically `c:\ado\personal`, but it might be different. In Stata, type `personal` to see the path.
* **Mac:** Typically `~/Documents/Stata/ado/personal` or `~/Library/Application Support/Stata/ado/personal`. In Stata, type `personal` to see the path.
* **Unix/Linux:** Typically `~/ado/personal`. In Stata, type `personal` to see the path.

For more details on personal ado directories, see the Stata FAQ:
https://www.stata.com/support/faqs/programming/personal-ado-directory/

After copying the files, you should be able to use the `clusterreg` command in Stata. You might need to type `discard` or restart Stata for it to recognize the new command.

This package requires the `moremata` package for some Mata functions (specifically `mm_median`). You can install it or ensure it's up-to-date by typing the following command in Stata:

ssc install moremata, replace

The package also relies heavily on embedded Mata functions for calculations. Ensure your Stata installation includes Mata.

<!-- **Important Note:** The individual method files (`imreg.ado`, `crsreg.ado`, `ccereg.ado`) rely on Mata helper functions defined at the end of the `clusterreg.ado` file and the `clpam.ado` program. If you intend to run these method files directly (without calling them through `clusterreg`), you must first ensure these Mata functions are defined in your Stata session (e.g., by running the `clusterreg.ado` file once, or by copying the Mata function block into each individual method file or a separate shared `.mata` file that is compiled) and also ensure `clpam.ado` is run or available. -->

--------------------------------------------------------------------------------
3. Usage
--------------------------------------------------------------------------------

The main command is `clusterreg`.

**Syntax:**

clusterreg depvar indepvars [if] [in] [weight] , coord(varlist) [ time(varname) type(string) ]

**Required Option:**

* `coord(varlist)`: Specifies the numeric variables representing the coordinates or features used for calculating distances and forming clusters (e.g., latitude, longitude, or other relevant characteristics). At least one variable must be provided.

**Optional Options:**

* `time(varname)`: Specifies a numeric variable indicating the time period for each observation. If provided, the estimation of the covariance structure (via QMLE) can incorporate temporal decay alongside spatial/feature decay. If omitted, only spatial/feature dependence is modeled in the covariance estimation step.
* `type(string)`: (default `"IM"`) Specifies the inference method to use. Options are:
    * `"IM"`: Ibragimov and Muller (2010) method (Default). [Runs `imreg.ado`]
    * `"CRS"`: Canay, Romano, and Shaikh (2017) randomization method. [Runs `crsreg.ado`]
    * `"CCE"`: Bester, Conley, and Hansen (2011) cluster covariance estimator method. [Runs `ccereg.ado`]
    
All methods (`IM`, `CRS`, `CCE`) first involve generating a set of potential clusterings of the data based on the provided `coord()` variables using k-medoids by the Partitioning Around Medoids algorithm (PAM). The number of clusters (`G`) ranges from 2 up to $ceil(n^{1/3})$.

The package then estimates the parameters of a spatial/spatio-temporal covariance function using Quasi-Maximum Likelihood Estimation (QMLE) on the OLS residuals.

Finally, for each potential number of clusters (`G`), the package simulates data under the null hypothesis (using the estimated covariance structure) and calculates the power of the chosen inference test (`IM`, `CRS`, or `CCE`) under specific local alternatives. Given correct size controlling, The number of clusters (`G*`) that yields the highest average simulated power is selected. The final inference (p-values, confidence intervals) is then performed using this chosen `G*` and the corresponding clustering.

--------------------------------------------------------------------------------
4. Examples
--------------------------------------------------------------------------------

Here are two examples:

**Example 1: Election Violence Data (Spatio-Temporal)**

clusterreg total_v2_agcho10 windspeed_06z windspeed_12z, coord(_cx _cy) time(first)

**Example 2: Homicide Data (Spatial Only)**

webuse homicide_1960_1990
keep if sname == "Florida"

clusterreg hrate divorce unemployment ln_income poverty, coord(_CX _CY) time(year) type("IM")

clusterreg hrate divorce unemployment ln_income poverty, coord(_CX _CY) type("CRS")

clusterreg hrate divorce unemployment ln_income poverty, coord(_CX _CY) type("CCE")