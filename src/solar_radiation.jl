"""
    twilight_irradiance(zenith_angle)

Compute twilight skylight irradiance for zenith angles between 88° and 107°.

Based on Rozenberg (1966) "Twilight" and Diem (1966) "Documenta Geigy Scientific Tables".

Returns `nothing` if zenith angle is outside twilight range.
"""
function twilight_irradiance(zenith_angle)
    z = zenith_angle
    if 88u"°" < z < 107u"°"
        log_illuminance = TWILIGHT_LOG_INTERCEPT - TWILIGHT_LOG_SLOPE * ustrip(u"°", z)
        return (10.0^log_illuminance) * LUX_TO_WATTS_PER_M2 * u"W/m^2"
    end
    return nothing
end

"""
    is_sun_up(time_from_noon, sunrise_hour_angle)

Check if the sun is above the horizon.

# Arguments
- `time_from_noon`: Hours from solar noon (negative=before, positive=after)
- `sunrise_hour_angle`: Hour angle at sunrise in hours

# Returns
`true` if sun is above horizon, `false` otherwise.
"""
function is_sun_up(time_from_noon, sunrise_hour_angle)
    ts = time_from_noon
    H₋ = sunrise_hour_angle

    if ts <= 0.0 && abs(ts) > H₋
        return false
    elseif ts > 0.0 && ts >= H₋
        return false
    end
    return true
end

"""
    slope_zenith_angle(zenith_angle, terrain, solar_azimuth) -> (; zenith_angle, cosine_zenith)

Calculate the effective zenith angle on the terrain surface.

On flat terrain (slope = 0°) returns the solar zenith angle unchanged.
On sloped terrain adjusts for slope and aspect using Eq. 3.15 of Sellers (1965).
"""
@inline function slope_zenith_angle(zenith_angle, terrain::AbstractTerrain, solar_azimuth)
    z = zenith_angle
    if terrain.slope > 0u"°"
        czsl = cos(z) * cos(terrain.slope) + sin(z) * sin(terrain.slope) * cos(solar_azimuth - terrain.aspect)
        zsl = acos(clamp(czsl, -1.0, 1.0))  # guard against floating-point overshoot past ±1
        zsl = min(uconvert(u"°", zsl), 90u"°")  # cap at 90° if sun is below slope horizon
    else
        czsl = cos(z)
        # `z` is a plain Float64 in radians (from `acos(cosZ)`); avoid the
        # `uconvert(u"°", z * u"rad")` round-trip which heap-allocates a
        # transient Quantity{rad} per call.
        zsl = (z * (180.0 / π)) * u"°"
    end
    return (; zenith_angle=zsl, cosine_zenith=czsl)
end

"""
    terrain_irradiance(global_horizontal, cosine_zenith, cosine_slope_zenith, zenith_angle, terrain)

Calculate terrain-adjusted global irradiance.

On flat terrain returns `global_horizontal` unchanged. On sloped terrain, scales
horizontal irradiance by the ratio of slope to horizontal cosine zenith when the
sun is above the horizon.
"""
function terrain_irradiance(global_horizontal, cosine_zenith, cosine_slope_zenith, zenith_angle, terrain::AbstractTerrain)
    terrain.slope > 0u"°" && zenith_angle < 90.0u"°" ?
        max(0.0u"W/m^2", (global_horizontal / cosine_zenith) * cosine_slope_zenith) : global_horizontal
end

"""
    refraction_correction(zenith_angle)

Apply atmospheric refraction correction to zenith angle.

Only applies for zenith angles > 88° (McCullough & Porter 1971).
"""
function refraction_correction(zenith_angle)
    z = zenith_angle
    if z < REFRACTION_ZENITH_THRESHOLD
        return z
    end
    # Note: McCullough & Porter (1971) give R = 16 + (Za - 88) * 10 arcminutes,
    # but this implementation uses a different coefficient (15 * 90/π ≈ 7.5 vs 10).
    # The formula may have been tuned empirically or accounts for Z vs Za difference.
    refraction = REFRACTION_BASE_ARCMIN + ((z - REFRACTION_ZENITH_THRESHOLD) * 15) / (π / 90)
    refraction = (refraction / 60) * (π / 180)
    return z - refraction
end

"""
    optical_air_mass(zenith_angle)

Calculate optical air mass using Rozenberg (1966) formula.

Reference: p.159 eq. III.3.17 in "Twilight" by Rozenberg (1966).
"""
function optical_air_mass(zenith_angle)
    z = zenith_angle
    return 1.0 / (cos(z) + AIR_MASS_A * exp(-AIR_MASS_B * cos(z)))
end

"""
    ozone_depth_lookup(latitude, day_of_year, year, ozone_column)

Look up ozone column depth from latitude/month table.

# Arguments
- `latitude`: Observer latitude (with angle units)
- `day_of_year`: Day of year (1-365)
- `year`: Year (for leap year handling)
- `ozone_column`: Lookup table (19×12 matrix)

# Returns
Ozone depth in cm.
"""
function ozone_depth_lookup(latitude, day_of_year, year, ozone_column)
    # Convert latitude to nearest 10-degree index
    lat_index = round(Int, (latitude + 100u"°") / 10u"°")
    lat_index = clamp(lat_index, 1, size(ozone_column, 1))
    mon = month(Date(year, 1, 1) + Day(day_of_year - 1))
    return ozone_column[lat_index, mon]
end

"""
    spectral_optical_depth(wavelength_index, atmospheric_pressure, mixing_ratio_height,
        rayleigh_optical_depth, ozone_optical_depth, aerosol_optical_depth, water_optical_depth,
        elevation_factors, ozone_depth, air_mass, precipitable_water)

Rayleigh and total optical depth at one wavelength (eqs. 13-14 McCullough & Porter 1971).

Returns `(; rayleigh, total)`.
"""
function spectral_optical_depth(
    wavelength_index, atmospheric_pressure, mixing_ratio_height,
    rayleigh_optical_depth, ozone_optical_depth, aerosol_optical_depth, water_optical_depth,
    elevation_factors, ozone_depth, air_mass, precipitable_water,
)
    (; molecular, aerosol, ozone, water) = elevation_factors

    n = wavelength_index
    P = atmospheric_pressure
    MR₀ = mixing_ratio_height
    τR = rayleigh_optical_depth
    τO = ozone_optical_depth
    τA = aerosol_optical_depth
    τW = water_optical_depth
    A₁ = molecular
    A₂ = aerosol
    A₃ = ozone
    A₄ = water
    X = ozone_depth
    m_Zₐ = air_mass
    w = precipitable_water

    λτR = (P / REFERENCE_PRESSURE) * τR[n] * A₁
    λτA = (REFERENCE_VISIBILITY / MR₀) * τA[n] * A₂
    λτO = (X / REFERENCE_OZONE_DEPTH_CM) * τO[n] * A₃
    λτW = τW[n] * sqrt(m_Zₐ * w * A₄) # eq. 13 in McCullough & Porter (1971)
    λτ = (λτR + λτA + λτO) * m_Zₐ + λτW # eq. 14
    λτ = min(λτ, MAX_OPTICAL_DEPTH) # avoids numerical issues at low sun angles
    return (; rayleigh=λτR, total=λτ)
end

"""
    direct_irradiance(solar_spectral_irradiance, sun_distance_factor, cosine_zenith,
        total_optical_depth, rayleigh_optical_depth, air_mass)

Direct and Rayleigh-only direct spectral irradiance at one wavelength, from the
solar spectral irradiance at that wavelength.

Returns `(; direct, rayleigh)` in W/m²/nm.
"""
function direct_irradiance(
    solar_spectral_irradiance, sun_distance_factor, cosine_zenith,
    total_optical_depth, rayleigh_optical_depth, air_mass,
)
    Sλ = solar_spectral_irradiance
    ar² = sun_distance_factor
    cosZ = cosine_zenith
    λτ = total_optical_depth
    λτR = rayleigh_optical_depth
    m_Zₐ = air_mass

    I₀λ = Sλ * ar² * cosZ # eq. 1 in McCullough & Porter (1971)
    T = λτ > 0.0 ? exp(-λτ) : 0.0
    Iλ = T < MIN_IRRADIANCE ? zero(I₀λ) : I₀λ * T # eq. 10
    Iλ = max(Iλ, MIN_IRRADIANCE * oneunit(I₀λ))
    Iᵣλ = I₀λ * exp(-λτR * m_Zₐ)

    return (; direct=Iλ, rayleigh=Iᵣλ)
end

"""
    horizon_angle_at_azimuth(solar_azimuth, horizon_angles)

Look up the horizon angle at a given solar azimuth.
"""
function horizon_angle_at_azimuth(solar_azimuth, horizon_angles)
    n = length(horizon_angles)
    step = 360u"°" / n
    best_idx = 1
    best_diff = abs(solar_azimuth - 0u"°")
    @inbounds for i in 2:n
        diff = abs(solar_azimuth - (i - 1) * step)
        if diff < best_diff
            best_diff = diff
            best_idx = i
        end
    end
    return @inbounds horizon_angles[best_idx]
end

"""
    trapezoidal_integrate!(direct_integral, rayleigh_integral, diffuse_integral, global_integral,
        direct_spectrum, rayleigh_spectrum, diffuse_spectrum, global_spectrum,
        wavelengths, wavelength_index)

One step of cumulative trapezoidal integration over wavelength, for each spectrum.
"""
function trapezoidal_integrate!(
    direct_integral, rayleigh_integral, diffuse_integral, global_integral,
    direct_spectrum, rayleigh_spectrum, diffuse_spectrum, global_spectrum,
    wavelengths, wavelength_index,
)
    ∫I = direct_integral
    ∫Iᵣ = rayleigh_integral
    ∫D = diffuse_integral
    ∫G = global_integral
    Iλ = direct_spectrum
    Iᵣλ = rayleigh_spectrum
    Dλ = diffuse_spectrum
    Gλ = global_spectrum
    λ = wavelengths
    n = wavelength_index

    if n == 1
        ∫D[1] = 0.0u"W/m^2"
        ∫Iᵣ[1] = 0.0u"W/m^2"
        ∫I[1] = 0.0u"W/m^2"
        ∫G[1] = 0.0u"W/m^2"
    else
        Δλ = λ[n] - λ[n-1]
        ∫I[n] = ∫I[n-1] + Δλ * Iλ[n-1] + 0.5Δλ * (Iλ[n] - Iλ[n-1])
        ∫Iᵣ[n] = ∫Iᵣ[n-1] + Δλ * Iᵣλ[n-1] + 0.5Δλ * (Iᵣλ[n] - Iᵣλ[n-1])
        ∫D[n] = ∫D[n-1] + Δλ * Dλ[n-1] + 0.5Δλ * (Dλ[n] - Dλ[n-1])
        ∫G[n] = ∫G[n-1] + Δλ * Gλ[n-1] + 0.5Δλ * (Gλ[n] - Gλ[n-1])
    end
end

"""
    allocate_output_arrays(nsteps, ndays, nmax)

Allocate all output arrays for solar radiation computation.
"""
function allocate_output_arrays(nsteps, ndays, nmax)
    return (;
        zenith_angle = fill(90.0u"°", nsteps),
        zenith_slope_angle = fill(90.0u"°", nsteps),
        # 90° is the sun-below-horizon sentinel; using a plain Vector{Quantity{°}}
        # avoids the Union{Missing,...} boxing that would heap-allocate every
        # azimuth write inside the per-step loop.
        azimuth_angle = fill(90.0u"°", nsteps),
        hour_angle_sunrise = fill(0.0, ndays),
        hour_solar_noon = fill(0.0, ndays),
        day_of_year = Vector{Int}(undef, nsteps),
        hour = Vector{Float64}(undef, nsteps),
        rayleigh_horizontal = fill(0.0u"W/m^2", nsteps),
        direct_horizontal = fill(0.0u"W/m^2", nsteps),
        diffuse_horizontal = fill(0.0u"W/m^2", nsteps),
        global_horizontal = fill(0.0u"W/m^2", nsteps),
        global_terrain = fill(0.0u"W/m^2", nsteps),
        rayleigh_spectra = fill(0.0u"W/nm/m^2", nsteps, nmax),
        direct_spectra = fill(0.0u"W/nm/m^2", nsteps, nmax),
        diffuse_spectra = fill(0.0u"W/nm/m^2", nsteps, nmax),
        global_spectra = fill(0.0u"W/nm/m^2", nsteps, nmax),
    )
end

"""
    sunrise_hour_angle(declination, latitude)

Calculate sunrise/sunset hour angles (eq.7 McCullough & Porter 1971).

Returns `(; cosine_hour_angle_sunset, hour_angle_sunset, hour_angle_sunrise)`:
- `cosine_hour_angle_sunset`: `-tan δ tan ϕ`, used to detect polar day and night
- `hour_angle_sunset`: hour angle at sunset (radians)
- `hour_angle_sunrise`: time of sunrise before solar noon (hours)
"""
function sunrise_hour_angle(declination, latitude)
    δ = declination
    ϕ = latitude

    cosH₊ = -tan(δ) * tan(ϕ) # eq. 7 in McCullough & Porter (1971)
    H₊ = abs(cosH₊) >= 1 ? float(π) : abs(acos(cosH₊))
    H₋ = 12.0 * H₊ / π

    return (;
        cosine_hour_angle_sunset = cosH₊,
        hour_angle_sunset = H₊,
        hour_angle_sunrise = H₋,
    )
end

"""
    compute_spectral_irradiance!(buffers, params, sun_below_horizon)

Compute spectral irradiance for all wavelengths at a single timestep.

Modifies `buffers` in place with computed spectral values.
"""
function compute_spectral_irradiance!(buffers, params::SpectralParams, sun_below_horizon)
    (; global_integral, rayleigh_integral, direct_integral, diffuse_integral,
       global_spectrum, rayleigh_spectrum, direct_spectrum, diffuse_spectrum) = buffers
    (; wavelength_count, atmospheric_pressure, mixing_ratio_height, rayleigh_optical_depth,
       ozone_optical_depth, aerosol_optical_depth, water_optical_depth, wavelengths,
       ozone_depth, precipitable_water, elevation_factors, diffuse_model, air_mass,
       solar_spectral_irradiance, sun_distance_factor, cosine_zenith) = params

    ∫G = global_integral
    ∫Iᵣ = rayleigh_integral
    ∫I = direct_integral
    ∫D = diffuse_integral
    Gλ = global_spectrum
    Iᵣλ = rayleigh_spectrum
    Iλ = direct_spectrum
    Dλ = diffuse_spectrum
    Sλ = solar_spectral_irradiance
    ar² = sun_distance_factor
    cosZ = cosine_zenith

    for n in 1:wavelength_count
        τ = spectral_optical_depth(n, atmospheric_pressure, mixing_ratio_height,
                                   rayleigh_optical_depth, ozone_optical_depth,
                                   aerosol_optical_depth, water_optical_depth,
                                   elevation_factors, ozone_depth, air_mass, precipitable_water)

        direct = direct_irradiance(Sλ[n], ar², cosZ, τ.total, τ.rayleigh, air_mass)
        Iλ[n] = direct.direct
        Iᵣλ[n] = direct.rayleigh

        if sun_below_horizon
            Iλ[n] = MIN_IRRADIANCE * u"W/m^2/nm"
            Iᵣλ[n] = MIN_IRRADIANCE * u"W/m^2/nm"
        end

        Dλ[n] = diffuse_irradiance(diffuse_model, n, τ.rayleigh, params, buffers)
        Gλ[n] = Dλ[n] + Iλ[n]

        trapezoidal_integrate!(∫I, ∫Iᵣ, ∫D, ∫G, Iλ, Iᵣλ, Dλ, Gλ, wavelengths, n)
    end
end

"""
    solar_radiation!(out, buffers, solar_model; kwargs...)

Mutating version of `solar_radiation` that writes results into pre-allocated buffers.

Use `allocate_output_arrays` and `allocate_buffers` to create the required buffers
for reuse across multiple calls.

See `solar_radiation` for argument documentation.
"""
function solar_radiation!(out, buffers, solar_model::AbstractSolarRadiation;
    solar_terrain::AbstractTerrain,
    days::Vector{<:Real}=[15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349],
    year::Real=1975,
    hours::AbstractVector{<:Real}=0:1:23,
    longitude_correction::Real=0.0,
    days_in_year::Real=365,
)
    (; solar_geometry_model, precipitable_water, diffuse_model, mixing_ratio_height,
       wavelength_count, wavelengths, ozone_column, rayleigh_optical_depth, ozone_optical_depth,
       aerosol_optical_depth, water_optical_depth, solar_spectral_irradiance) = solar_model
    (; elevation, albedo, atmospheric_pressure, latitude) = solar_terrain

    nmax = wavelength_count
    λ = wavelengths
    τR = rayleigh_optical_depth
    τO = ozone_optical_depth
    τA = aerosol_optical_depth
    τW = water_optical_depth
    Sλ = solar_spectral_irradiance
    w = precipitable_water
    MR₀ = mixing_ratio_height
    ϕ = latitude
    P = atmospheric_pressure
    A = albedo

    elevation_factors = elevation_correction(elevation)
    step = 1
    H₋ = 0.0
    tsn = 0.0

    for i in eachindex(days)
        for j in eachindex(hours)
            d = days[i]
            t = hours[j]
            h, tsn = hour_angle(t, longitude_correction)
            (; solar_declination, zenith_angle, sun_distance_factor) =
                solar_geometry(solar_geometry_model, ϕ; day_of_year=d, hour_angle=h, days_in_year)
            δ = solar_declination
            z = zenith_angle
            ar² = sun_distance_factor
            zsl = z

            # Twilight handling
            skylight = twilight_irradiance(z)
            if skylight !== nothing
                buffers.diffuse_integral[nmax] = skylight
                buffers.global_integral[nmax] = skylight
                out.global_horizontal[step] = skylight
                out.diffuse_horizontal[step] = skylight
            end

            # Sunrise/sunset calculation
            (; cosine_hour_angle_sunset, hour_angle_sunrise) = sunrise_hour_angle(δ, ϕ)
            H₋ = hour_angle_sunrise
            sun_up = is_sun_up(t - tsn, H₋)

            # 90° = sun-below-horizon sentinel (matches the Vector pre-fill);
            # avoids `Union{Missing, Quantity{°}}` typing of `solar_azimuth`
            # which would force boxing through the slope_zenith_angle call below.
            solar_azimuth = 90.0u"°"
            if sun_up || cosine_hour_angle_sunset == 1
                alt = (π / 2 - z)u"rad"
                solar_azimuth = solar_azimuth_angle(h, ϕ, δ)
                ahoriz = horizon_angle_at_azimuth(solar_azimuth, solar_terrain.horizon_angles)

                # Slope geometry
                slope_geometry = slope_zenith_angle(z, solar_terrain, solar_azimuth)
                zsl = slope_geometry.zenith_angle
                czsl = slope_geometry.cosine_zenith

                # Refraction correction
                z = refraction_correction(z)
                cosZ = cos(z)
                intcz = floor(Int, 100.0 * cosZ + 1.0)
                m_Zₐ = optical_air_mass(z)
                X = ozone_depth_lookup(ϕ, d, year, ozone_column)

                # Compute spectral irradiance
                params = SpectralParams(
                    nmax, P, MR₀, τR, τO, τA, τW, Sλ, λ, ar², cosZ, intcz, m_Zₐ,
                    X, w, elevation_factors, diffuse_model, A, z
                )
                compute_spectral_irradiance!(buffers, params, alt < ahoriz)

                # Store results
                (; global_integral, rayleigh_integral, direct_integral, diffuse_integral,
                   global_spectrum, rayleigh_spectrum, direct_spectrum, diffuse_spectrum) = buffers
                out.global_spectra[step, :] .= global_spectrum
                out.rayleigh_spectra[step, :] .= rayleigh_spectrum
                out.direct_spectra[step, :] .= direct_spectrum
                out.diffuse_spectra[step, :] .= diffuse_spectrum
                out.global_horizontal[step] = global_integral[nmax]
                out.global_terrain[step] = terrain_irradiance(global_integral[nmax], cosZ, czsl, z, solar_terrain)
                out.rayleigh_horizontal[step] = rayleigh_integral[nmax]
                out.direct_horizontal[step] = direct_integral[nmax]
                out.diffuse_horizontal[step] = diffuse_integral[nmax]
            end

            out.zenith_angle[step] = uconvert(u"°", z)
            out.zenith_slope_angle[step] = uconvert(u"°", zsl)
            out.azimuth_angle[step] = uconvert(u"°", solar_azimuth)
            out.day_of_year[step] = d
            out.hour[step] = t
            step += 1
        end
        out.hour_angle_sunrise[i] = H₋
        out.hour_solar_noon[i] = tsn
    end

    return out
end

"""
    solar_radiation(solar_model; solar_terrain, days, hours, ...)
    solar_radiation(solar_model; solar_terrain, dates, hours, ...)

Compute solar radiation for a given model and terrain configuration.

Allocates output buffers internally. For repeated calls with the same dimensions,
use `solar_radiation!` with pre-allocated buffers for better performance.

# Arguments
- `solar_model::AbstractSolarRadiation`: Solar radiation model parameters
- `solar_terrain::AbstractTerrain`: Terrain configuration (elevation, slope, etc.)

# Numeric interface (R-compatible)
- `days`: Vector of days of year (default: mid-month days)
- `hours`: Hours of day to compute (default: 0:23)
- `year`: Year for leap year handling (default: 1975)
- `longitude_correction`: Longitude correction in hours (default: 0.0)
- `days_in_year`: Number of days in the year (default: 365). Use 366 for leap years,
   or other values for non-standard calendars (e.g., 360 for 360-day calendar).

# DateTime interface
- `dates`: Vector of `AbstractDateTime` instances (one per day to simulate)
- `hours`: Hour offsets as `Period` values (e.g., `Hour(0):Hour(1):Hour(23)`)
- `timezone_offset`: Timezone offset from UTC in hours (default: 0.0)

The DateTime interface extracts calendar information (leap years, 360-day calendars)
from the date type automatically.

# Examples
```julia
# Numeric interface
solar_radiation(model; solar_terrain, days=1:10, hours=0:23, year=2024)

# DateTime interface - mid-month days for a year
using Dates
dates = [Date(2024, m, 15) for m in 1:12]
solar_radiation(model; solar_terrain, dates, hours=Hour(0):Hour(1):Hour(23))

# With CFTime (no dependency required)
using CFTime
dates = [DateTimeNoLeap(2024, m, 15) for m in 1:12]
solar_radiation(model; solar_terrain, dates, hours=Hour(0):Hour(1):Hour(23))
```

# Returns
NamedTuple with zenith/azimuth angles, integrated irradiances, and spectral data.
"""
function solar_radiation(solar_model::AbstractSolarRadiation;
    solar_terrain::AbstractTerrain,
    # Numeric interface
    days::Union{Nothing, AbstractVector{<:Real}}=nothing,
    year::Real=1975,
    hours::AbstractVector=0:1:23,
    longitude_correction::Real=0.0,
    days_in_year::Real=365,
    # DateTime interface (accepts Date, DateTime, or CFTime types)
    dates::Union{Nothing, AbstractVector{<:Dates.TimeType}}=nothing,
    timezone_offset::Real=0.0,
)
    # Determine which interface is being used
    if !isnothing(dates)
        # DateTime interface: extract from dates + hour offsets
        days_in_year = Dates.daysinyear(first(dates))
        year = Dates.year(first(dates))
        longitude_correction = timezone_offset
        # Convert dates to day-of-year, hours to numeric
        days_numeric = [Dates.dayofyear(d) for d in dates]
        hours_numeric = [Dates.value(Dates.Millisecond(h)) / 3_600_000 for h in hours]
    elseif !isnothing(days)
        # Numeric interface
        days_numeric = collect(days)
        hours_numeric = collect(hours)
    else
        # Default mid-month days
        days_numeric = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
        hours_numeric = collect(hours)
    end

    nmax = solar_model.wavelength_count
    ndays = length(days_numeric)
    nsteps = ndays * length(hours_numeric)

    out = allocate_output_arrays(nsteps, ndays, nmax)
    buffers = allocate_buffers(nmax, solar_model.diffuse_model)

    return solar_radiation!(out, buffers, solar_model;
        solar_terrain,
        days=days_numeric,
        year,
        hours=hours_numeric,
        longitude_correction,
        days_in_year)
end
