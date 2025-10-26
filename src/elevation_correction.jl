"""
    elevation_correction(elevation)

Calculates smooth polynomial approximations of atmospheric constituent correction factors 
(molecular, aerosol, ozone and water vapour) as a function of elevation, derived from data
in Table 3 of McCullough & Porter (1971). 

Water vapour is included for completeness but no standard profile is available so
the correction factor is 1.0

# Returns

A `NamedTuple` with the fields:
- `molecular`
- `aerosol`
- `ozone`
- `water`
"""
function elevation_correction(elevation)
    elev_km = ustrip(u"km", elevation + 1.0u"km")

    molecular =0.00007277 * elev_km^3 +
               0.00507293 * elev_km^2 -
               0.12482149 * elev_km +
               1.11687469

    aerosol =  8.35656e-7 * elev_km^6 -
               6.26384e-5 * elev_km^5 +
               1.86967e-3 * elev_km^4 -
               2.82585e-2 * elev_km^3 +
               2.26739e-1 * elev_km^2 -
               9.25268e-1 * elev_km +
               1.71321

    ozone =    1.07573e-6 * elev_km^5 -
               5.14511e-5 * elev_km^4 +
               7.97960e-4 * elev_km^3 -
               4.90904e-3 * elev_km^2 +
               2.99258e-3 * elev_km +
               1.00238

    water = 1.0

    return (; molecular, aerosol, ozone, water)
end