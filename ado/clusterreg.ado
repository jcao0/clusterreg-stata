program define clusterreg
    syntax anything [if] [in] [aweight fweight iweight pweight /], coord(varlist) [time(varname) type(string)]
    
    capture program drop clpam
    mata: mata clear
    matrix dissim dd = `coord', L2Squared
    
    if "`type'"=="CCE"{
        capture program drop cce
        if "`time'" != "" {
            ccereg `anything' `if' `in' `weight', coord(`coord') timeperiod(`time')
        }
        else {
            ccereg `anything' `if' `in' `weight', coord(`coord')
        }
    }
    else if "`type'"=="CRS"{
        capture program drop crs
        if "`time'" != "" {
            crsreg `anything' `if' `in' `weight', coord(`coord') timeperiod(`time')
        }
        else {
            crsreg `anything' `if' `in' `weight', coord(`coord')
        }
    }
    else{
        capture program drop im
        if "`time'" != "" {
            imreg `anything' `if' `in' `weight', coord(`coord') timeperiod(`time')
        }
        else {
            imreg `anything' `if' `in' `weight', coord(`coord')
        }
    }
    
    capture drop id0
    capture drop group*
    
end
