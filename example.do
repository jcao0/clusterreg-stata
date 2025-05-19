// Example 1: Election Violence Data (Spatio-Temporal)

clear all
set more off

// Import example data (replace with actual path/source)
import delimited ElectionViolenceTABLE2.csv

clusterreg total_v2_agcho10 windspeed_06z windspeed_12z, coord(_cx _cy) time(first) 



// Example 2: Homicide Data (Spatial Only)

// Load necessary components
clear all
set more off

// Use Stata's built-in webuse data
webuse homicide_1960_1990
keep if sname == "Florida"

// Run regression using IM method (default) with spatial coordinates (_CX, _CY)
clusterreg hrate divorce unemployment ln_income poverty, coord(_CX _CY) time(year) type("IM")

// Alternatively, run using the CRS method, with spatial coordinates only:
clusterreg hrate divorce unemployment ln_income poverty, coord(_CX _CY) type("CRS")

// Or the CCE method:
clusterreg hrate divorce unemployment ln_income poverty, coord(_CX _CY) type("CCE")
