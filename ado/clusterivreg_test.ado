*! clusterivreg
*! Dispatches to learned-cluster IV routines (IM, CRS, CCE).
program define clusterivreg_test, eclass
    version 17
    syntax anything [if] [in] [aweight fweight iw pw], ///
           coord(varlist) [time(varname) type(string)]

    local allargs `anything'
    
    local iv_spec ""
    if regexm("`allargs'", "\(([^)]+)\)") {
        local iv_spec = regexs(1)
        local main_vars = regexr("`allargs'", "\s*\(([^)]+)\)\s*", " ")
        local main_vars = trim("`main_vars'")
    }
    else {
        di as error "IV specification in parentheses, e.g., (endog = inst), is required."
        exit 198
    }
    
    // --- CCE-IV ---
    if ("`type'" == "CCE" | "`type'" == "cce") {
        di as text "--> Calling CCE-IV method"
        cceivreg_test `main_vars' `if' `in' `weight', iv(`iv_spec') cluster(`coord') timeperiod(`time')
    }
    
    // --- CRS-IV ---
    else if ("`type'" == "CRS" | "`type'" == "crs") {
        di as text "--> Calling CRS-IV method"
        crsivreg_test `main_vars' `if' `in' `weight', iv(`iv_spec') cluster(`coord') timeperiod(`time')
    }
    
    // --- IM-IV (Default) ---
    else {
        if ("`type'"=="" | "`type'"=="IM" | "`type'"=="im") {
            di as text "--> Calling IM-IV method (default)"
            imivreg_test `main_vars' `if' `in' `weight', iv(`iv_spec') cluster(`coord') timeperiod(`time')
        }
        else {
            di as error "type(`type') not recognized. Available types are IM, CRS, CCE."
            exit 198
        }
    }

    capture drop id0
    capture drop group*

end