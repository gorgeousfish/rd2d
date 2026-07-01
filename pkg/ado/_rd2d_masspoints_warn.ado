*! version 1.1.0 24jun2026
*! _rd2d_masspoints_warn.ado
*! Shared mass-point warning utility for the rd2d package.
*! Scans a massinfo matrix (neval x 4) for mass_ratio >= 0.2 and displays
*! the standard mass-point warning messages.
*!
*! Usage (from ado-level Stata code):
*!   _rd2d_masspoints_warn, matrix(`massinfo') neval(`neval') masspoints(`masspoints')
*!
*! Interface contract:
*!   massinfo matrix layout - each row is one evaluation point:
*!     col 1 = M   (total unique values)
*!     col 2 = M0  (unique values on control/negative side)
*!     col 3 = M1  (unique values on treated/positive side)
*!     col 4 = mass_ratio  (1 - M/N; proportion of duplicated locations)
*!   When masspoints="off", rows carry M=N, M0=N0, M1=N1, mass_ratio=0.

program define _rd2d_masspoints_warn
    version 16
    syntax , MATrix(name) Neval(integer) Masspoints(string)

    // Only check/adjust paths need warnings; "off" is a no-op.
    if !inlist("`masspoints'", "check", "adjust") exit

    local has_mass = 0
    forvalues j = 1/`neval' {
        if (`matrix'[`j', 4] >= 0.2) local has_mass = 1
    }
    if (`has_mass') {
        di as txt "warning: mass points detected in the running variables."
        if ("`masspoints'" == "check") {
            di as txt "warning: try using option masspoints(adjust)."
        }
    }
end
