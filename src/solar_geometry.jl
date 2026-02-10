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
    solar_geometry(latitude; day_of_year, hour_angle)

Computes key solar geometry parameters based on McCullough & Porter (1971).

# Arguments
- `latitude`: Latitude (with angle units, e.g. `u"°"` or `u"rad"`)
- `day_of_year`: Day of year (1–365)
- `hour_angle`: Hour angle (radians)

# Returns
NamedTuple with:
- `solar_longitude`: Auxiliary solar longitude (radians)
- `solar_declination`: Solar declination (radians)
- `zenith_angle`: Solar zenith angle (radians)
- `sun_distance_factor`: Square of Earth-to-Sun radius factor (unitless)

# Reference
McCullough & Porter (1971)
"""
@kwdef struct McCulloughPorterSolarGeometry <: AbstractSolarGeometryModel
    reference_day::Real = 80
    orbital_angular_frequency::Real = 2π / 365
    orbital_eccentricity::Real = 0.0167238
    declination_amplitude::Real = 0.39784993
end

solar_geometry(::McCulloughPorterSolarGeometry, ::Missing; kwargs...) = missing
function solar_geometry(sm::McCulloughPorterSolarGeometry, latitude::Quantity;
    day_of_year::Real,
    hour_angle::Quantity,
)
    (; reference_day, orbital_angular_frequency, orbital_eccentricity, declination_amplitude) = sm
    # Use short aliases for equations (standard notation)
    d0, ω, ϵ, se = reference_day, orbital_angular_frequency, orbital_eccentricity, declination_amplitude
    d, h = day_of_year, hour_angle

    ζ = (ω * (d - d0)) + 2.0ϵ * (sin(ω * d) - sin(ω * d0))          # eq.5 McCullough & Porter (1971)
    δ = asin(se * sin(ζ))                                           # eq.4 McCullough & Porter (1971)
    cosZ = cos(latitude) * cos(δ) * cos(h) + sin(latitude) * sin(δ) # Eq.3 McCullough & Porter (1971)
    z = acos(cosZ)
    ar² = 1.0 + (2.0ϵ) * cos(ω * d)                                 # eq.2 McCullough & Porter (1971)
    return (;
        solar_longitude = ζ,
        solar_declination = δ,
        zenith_angle = z,
        sun_distance_factor = ar²
    )
end