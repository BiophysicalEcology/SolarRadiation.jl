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
    tsn = 12.0 + longitude_correction # time of solar noon
    h = 15u"°" * (t - tsn) # eq. 6 in McCullough & Porter (1971), 15° per hour
    return uconvert(u"rad", h), tsn
end

abstract type AbstractSolarGeometryModel end

"""
    orbital_angular_frequency(days_in_year::Real=365)

Compute Earth's orbital angular frequency for a given year length.

Handles leap years (366 days) and non-standard calendars (e.g., 360-day).
"""
orbital_angular_frequency(days_in_year::Real=365) = 2π / days_in_year

"""
    McCulloughPorterSolarGeometry

Solar geometry model based on McCullough & Porter (1971).

# Fields
- `reference_day`: Vernal equinox day of year (default: 80)
- `orbital_eccentricity`: Earth's orbital eccentricity (default: 0.0167238)
- `declination_amplitude`: Solar declination amplitude (default: 0.39784993)

Note: `orbital_angular_frequency` is now computed dynamically based on the year length
to handle leap years and non-standard calendars. Use `orbital_angular_frequency(days_in_year)`.
"""
@kwdef struct McCulloughPorterSolarGeometry{RD,OE,DA} <: AbstractSolarGeometryModel
    reference_day::RD = 80
    orbital_eccentricity::OE = 0.0167238
    declination_amplitude::DA = 0.39784993
end

"""
    solar_geometry(model::McCulloughPorterSolarGeometry, latitude; day_of_year, hour_angle, days_in_year=365)

Compute solar geometry parameters based on McCullough & Porter (1971).

# Arguments
- `model`: Solar geometry model with orbital parameters
- `latitude`: Observer latitude (with angle units, e.g. `u"°"` or `u"rad"`)
- `day_of_year`: Day of year (1–365, 1–366 for leap years, or 1–360 for 360-day calendars)
- `hour_angle`: Hour angle (radians)
- `days_in_year`: Number of days in the year (default: 365). Use 366 for leap years,
   or other values for non-standard calendars (e.g., 360 for 360-day calendar).

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
    days_in_year::Real=365,
)
    (; reference_day, orbital_eccentricity, declination_amplitude) = sm

    d₀ = reference_day
    ϵ = orbital_eccentricity
    se = declination_amplitude
    ϕ = latitude
    d = day_of_year
    h = hour_angle
    ω = orbital_angular_frequency(days_in_year)

    ζ = ω * (d - d₀) + 2ϵ * (sin(ω * d) - sin(ω * d₀)) # eq. 5 in McCullough & Porter (1971)
    δ = asin(se * sin(ζ)) # eq. 4
    cosZ = cos(ϕ) * cos(δ) * cos(h) + sin(ϕ) * sin(δ) # eq. 3
    Z = acos(cosZ)
    ar² = 1 + 2ϵ * cos(ω * d) # eq. 2

    return (;
        solar_longitude = ζ,
        solar_declination = δ,
        zenith_angle = Z,
        sun_distance_factor = ar²,
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
    h = hour_angle
    ϕ = latitude
    δ = declination

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
