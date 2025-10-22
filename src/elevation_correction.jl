"""
    elevation_correction(elevation)

Calculates smooth polynomial approximations of atmospheric constituent correction factors 
as a function of elevation (based on Kearney's modification of the ALTFCT array originally 
from SOLAR.DAT). Input `elevation` is the elevation in meters and can include units.

# Description


The array `ELEVFCT(i, j)` represents the **ratio of the total amount of a given 
atmospheric constituent (index j) above the elevation of interest (index i) to that 
above sea level**. The constituent indices are:

- j = 1: Molecular
- j = 2: Aerosol
- j = 3: Ozone
- j = 4: Water vapor

For j = 1–3, values are derived from standard profiles. For water vapor (j = 4), no 
standard profile exists, so only `ELEVFCT(1, 4)` is defined as 1.00.

The elevation index i runs from 1 to 21, corresponding to elevations from sea level to 
20 km in 1 km steps.

This function implements fitted polynomials to reproduce this correction smoothly 
from `elevation` (in meters) using continuous approximation.

# Returns

A `NamedTuple` with the fields:
- `molecular_corr`
- `aerosol_corr`
- `ozone_corr`
- `water_vapor_corr`
"""
function elevation_correction(elevation)
    elev_km = ustrip(u"km", elevation + 1.0u"km")

    molecular_corr =0.00007277 * elev_km^3 +
               0.00507293 * elev_km^2 -
               0.12482149 * elev_km +
               1.11687469

    aerosol_corr =  8.35656e-7 * elev_km^6 -
               6.26384e-5 * elev_km^5 +
               1.86967e-3 * elev_km^4 -
               2.82585e-2 * elev_km^3 +
               2.26739e-1 * elev_km^2 -
               9.25268e-1 * elev_km +
               1.71321

    ozone_corr =    1.07573e-6 * elev_km^5 -
               5.14511e-5 * elev_km^4 +
               7.97960e-4 * elev_km^3 -
               4.90904e-3 * elev_km^2 +
               2.99258e-3 * elev_km +
               1.00238

    water_vapour_corr = 1.0

    return (; molecular_corr, aerosol_corr, ozone_corr, water_vapour_corr)
end