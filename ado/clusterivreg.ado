*! clusterivreg.ado
*! Version 2.2: Final dispatcher for learned-cluster IV routines.
*! Dispatches to imivreg, crsivreg, and cceivreg.
program define clusterivreg, eclass
    version 17
    syntax anything [if] [in] [aweight fweight iw pw], ///
           coord(varlist) [TIME(varname) TYPE(string)]

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
        di as text "--> Calling CCE-IV method with learned clusters"
        // Call the finalized cceivreg program
        cceivreg `main_vars' `if' `in' `weight', iv(`iv_spec') coord(`coord') timeperiod(`time')
    }
    
    // --- CRS-IV ---
    else if ("`type'" == "CRS" | "`type'" == "crs") {
        di as text "--> Calling CRS-IV method with learned clusters"
        // Call the finalized crsivreg program
        crsivreg `main_vars' `if' `in' `weight', iv(`iv_spec') coord(`coord') timeperiod(`time')
    }
    
    // --- IM-IV (Default) ---
    else {
        if ("`type'"=="" | "`type'"=="IM" | "`type'"=="im") {
            di as text "--> Calling IM-IV method with learned clusters (default)"
            // Call the finalized imivreg program
            imivreg `main_vars' `if' `in' `weight', iv(`iv_spec') coord(`coord') timeperiod(`time')
        }
        else {
            di as error "type(`type') not recognized. Available types are IM, CRS, CCE."
            exit 198
        }
    }

    // Note: The cleanup of temp variables is now handled within each sub-program.
end