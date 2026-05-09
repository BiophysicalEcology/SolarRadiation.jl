"""
    hour_angle(t::Real, longitude_correction::Real=0) -> (h, tsn)

Compute the solar hour angle `h` in radians.

# Arguments
- `t`: Local solar hour (e.g., `14.0`)
- `longitude_correction`: Longitude correction in hours (e.g., `0.5`), default `0`

# Returns
- `h`: Hour angle as a `Quantity` in radians
- `tsn`: Time at solar noon in hours

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
    McCulloughPorterSolarGeometry

Solar geometry model based on McCullough & Porter (1971).

# Fields
- `reference_day`: Vernal equinox day of year (default: 80)
- `orbital_angular_frequency`: Earth's orbital angular frequency (default: 2π/365)
- `orbital_eccentricity`: Earth's orbital eccentricity (default: 0.0167238)
- `declination_amplitude`: Solar declination amplitude (default: 0.39784993)
"""
@kwdef struct McCulloughPorterSolarGeometry{RD,OAF,OE,DA} <: AbstractSolarGeometryModel
    reference_day::RD = 80
    orbital_angular_frequency::OAF = 2π / 365
    orbital_eccentricity::OE = 0.0167238
    declination_amplitude::DA = 0.39784993
end

"""
    solar_geometry(model::McCulloughPorterSolarGeometry, latitude; day_of_year, hour_angle)

Compute solar geometry parameters based on McCullough & Porter (1971).

# Arguments
- `model`: Solar geometry model with orbital parameters
- `latitude`: Observer latitude (with angle units, e.g. `u"°"` or `u"rad"`)
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

"""
    solar_azimuth_angle(hour_angle, latitude, declination)

Compute solar azimuth angle with quadrant correction.

# Arguments
- `hour_angle`: Hour angle (radians, negative=morning, positive=afternoon)
- `latitude`: Observer latitude (with angle units)
- `declination`: Solar declination (radians)

# Returns
Solar azimuth angle in degrees (0-360°, measured clockwise from north)
"""
function solar_azimuth_angle(hour_angle, latitude, declination)
    h, ϕ, δ = hour_angle, latitude, declination

    tan_azimuth = sin(h) / (cos(ϕ) * tan(δ) - sin(ϕ) * cos(h))
    azimuth = atan(tan_azimuth)

    # Correct for hemisphere/quadrant
    azimuth = if h == 0.0 # Special case: solar noon
        # Sun is south when observer is poleward of sub-solar point (denominator < 0
        # at h=0), north otherwise.
        (cos(ϕ) * tan(δ) - sin(ϕ)) < 0 ? 180.0u"°" : 0.0u"°"
    elseif h <= 0.0 # Morning - east of reference
        if azimuth <= 0.0u"°"
            -azimuth  # 1st Quadrant (0–90°)
        else
            180.0u"°" - azimuth  # 2nd Quadrant (90–180°)
        end
    else # Afternoon - west of reference
        if azimuth < 0.0u"°"
            180.0u"°" - azimuth  # 3rd Quadrant (180–270°)
        else
            360.0u"°" - azimuth  # 4th Quadrant (270–360°)
        end
    end

    return azimuth
end
