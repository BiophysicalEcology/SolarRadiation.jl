"""
    hour_angle(t::Quantity, longitude_correction::Quantity) -> Quantity

Compute the solar hour angle `h` in radians.

# Arguments
- `t`: Local solar hour (e.g., `14.0`)
- `longitude_correction`: Longitude correction in hours (e.g., `0.5`)

# Returns

- Hour angle `h` as a `Quantity` in radians
- Time at solar noon, `tsn` as a time in hours

# Reference
McCullough & Porter 1971, Eq. 6
"""
function hour_angle(t::Real, longitude_correction::Real=0)
    tsn = 12.0 + longitude_correction      # solar noon time
    h = (π / 12) * (t - tsn) * u"rad"      # convert hours to radians
    return h, tsn
end

abstract type AbstractSolarGeometryModel end

"""
    solar_geometry(d::Real, latitude::Quantity, h::Quantity; d0::Real = 80, ω::Real = 2π/365, ϵ::Real = 0.0167, se::Real = 0.39779)

Computes key solar geometry parameters based on McCullough & Porter (1971):

- `ζ`: Auxiliary solar longitude (radians)
- `δ`: Solar declination (radians)
- `z`: Solar zenith angle (radians)
- `ar²`: Square of Earth-to-Sun radius factor (unitless)

# Arguments
- `d`: Day of year (1–365)
- `latitude`: Latitude (with angle units, e.g. `u"°"` or `u"rad"`)
- `h`: Hour angle (radians)

- `d0`: Reference day (default: 80)
- `ω`: Angular frequency of Earth’s orbit (default: `2π/365`)
- `ϵ`: Orbital eccentricity (default: `0.0167`)
- `se`: Constant for solar declination amplitude (default: `0.39779`)

# Returns
Tuple: `(ζ, δ, z, ar²)` with angle quantities in radians and ar² unitless.

# Reference
McCullough & Porter (1971)
"""
@kwdef struct McCulloughPorterSolarGeometry <: AbstractSolarGeometryModel
    d0::Real = 80
    ω::Real = 2π / 365
    ϵ::Real = 0.0167238
    se::Real = 0.39784993 #0.39779
end

solar_geometry(::McCulloughPorterSolarGeometry, ::Missing, ; kwargs...) = missing
function solar_geometry(sm::McCulloughPorterSolarGeometry, latitude::Quantity; # =83.07305u"°",
    d::Real, # =1.0,
    h::Quantity, # =-2.87979u"rad",
)
    (; d0, ω, ϵ, se) = sm

    ζ = (ω * (d - d0)) + 2.0ϵ * (sin(ω * d) - sin(ω * d0))          # eq.5 McCullough & Porter (1971)
    δ = asin(se * sin(ζ))                                           # eq.4 McCullough & Porter (1971)
    cosZ = cos(latitude) * cos(δ) * cos(h) + sin(latitude) * sin(δ) # Eq.3 McCullough & Porter (1971)
    z = acos(cosZ)u"rad"                                          
    AR2 = 1.0 + (2.0ϵ) * cos(ω * d)                                 # eq.2 McCullough & Porter (1971)
    δ = δ * u"rad"
    ζ = ζ * u"rad"
    return(; ζ, δ, z, ar²)
end