"""
    twilight_irradiance(zenith_angle)

Compute twilight skylight irradiance for zenith angles between 88° and 107°.

Based on Rozenberg (1966) "Twilight" and Diem (1966) "Documenta Geigy Scientific Tables".

Returns `nothing` if zenith angle is outside twilight range.
"""
function twilight_irradiance(z)
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
    ts, H₋ = time_from_noon, sunrise_hour_angle
    if ts <= 0.0 && abs(ts) > H₋
        return false
    elseif ts > 0.0 && ts >= H₋
        return false
    end
    return true
end

"""
    slope_zenith_angle(zenith, slope, solar_azimuth, aspect)

Calculate the effective zenith angle on a sloped surface.

Based on Eq. 3.15 in Sellers (1965) "Physical Climatology".

# Returns
NamedTuple with `(; zenith_angle, cosine_zenith)`
"""
function slope_zenith_angle(z, slope, solar_azimuth, aspect)
    if slope > 0u"°"
        czsl = cos(z) * cos(slope) + sin(z) * sin(slope) * cos(solar_azimuth - aspect)
        zsl = acos(czsl)
        zsl = min(uconvert(u"°", zsl), 90u"°")  # cap at 90° if sun is below slope horizon
    else
        czsl = cos(z)
        zsl = z
    end
    return (; zenith_angle=zsl, cosine_zenith=czsl)
end

"""
    refraction_correction(zenith_angle)

Apply atmospheric refraction correction to zenith angle.

Only applies for zenith angles > 88° (McCullough & Porter 1971).
"""
function refraction_correction(z)
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
function optical_air_mass(z)
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
    diffuse_irradiance(wavelength_index, scattered, scattered_uv, λτR, gamma_buffers,
                       cz, intcz, Sλ_n, ar², m_Zₐ, albedo, z, FD, FDQ, s̄)

Calculate diffuse (scattered) spectral irradiance for a single wavelength.

# Arguments
- `wavelength_index`: Index into wavelength array (1-111)
- `scattered`: If false, disable all scattering
- `scattered_uv`: If true, use full Chandrasekhar model; otherwise use lookup tables
- `λτR`: Rayleigh optical depth at this wavelength
- `gamma_buffers`: Pre-allocated buffers for Chandrasekhar calculation
- `cz`: Cosine of zenith angle
- `intcz`: Integer index for zenith (1-101)
- `Sλ_n`: Solar spectral irradiance at this wavelength
- `ar²`: Sun distance factor
- `m_Zₐ`: Optical air mass
- `albedo`: Surface albedo
- `z`: Zenith angle
- `FD`, `FDQ`: Dave & Furukawa lookup tables
- `s̄`: Single scattering albedo (vector or scalar depending on mode)

# Returns
Diffuse spectral irradiance in W/m²/nm
"""
function diffuse_irradiance(n, scattered, scattered_uv, λτR, gamma_buffers,
                            cz, intcz, Sλ_n, ar², m_Zₐ, albedo, z, FD, FDQ, s̄)
    if !scattered
        return 0.0u"W/m^2/nm"
    end

    if scattered_uv
        # Full Chandrasekhar scattering model
        if λτR >= MIN_RAYLEIGH_OPTICAL_DEPTH_CHANDRASEKHAR
            γᵣ, γₗ, s̄_scalar = scattered_radiation!(gamma_buffers, λτR)
            I₀_λ = cz * Sλ_n * ar² / 1000.0
            # eq. 15 McCullough & Porter 1971
            return (((float(γₗ[intcz]) + float(γᵣ[intcz])) / (2.0 * (1.0 - albedo * float(s̄_scalar))))
                    - exp(-float(λτR) * m_Zₐ)) * I₀_λ
        else
            return 0.0u"W/m^2/nm"
        end
    else
        # Dave & Furukawa lookup table method (UV wavelengths only)
        if n > 11
            return 0.0u"W/m^2/nm"
        end
        B = ustrip(u"°", z) / 5
        k = trunc(Int, B) + 1 + (B % 1 > 0.5)
        flux_down = FD[n, k]
        flux_down_div_Q = FDQ[n, k]
        Q = albedo / (1.0 - albedo * s̄[n])  # eq. 31 in Dave & Furukawa 1966
        Dλ = (Sλ_n / π) * (flux_down + flux_down_div_Q * Q) / 1000.0
        return Dλ * ar²
    end
end

"""
    spectral_optical_depth(n, P, MR₀, τR, τO, τA, τW, A1, A2, A3, A4, ozone_depth, m_Zₐ, cmH2O)

Calculate wavelength-specific optical depths and total optical depth.

Returns `(; rayleigh, total)` - Rayleigh optical depth and total optical depth.
"""
function spectral_optical_depth(n, P, MR₀, τR, τO, τA, τW, A1, A2, A3, A4, ozone_depth, m_Zₐ, cmH2O)
    λτR = (P / REFERENCE_PRESSURE) * τR[n] * A1
    λτA = (REFERENCE_VISIBILITY / MR₀) * τA[n] * A2
    λτO = (ozone_depth / REFERENCE_OZONE_DEPTH_CM) * τO[n] * A3
    λτW = τW[n] * sqrt(m_Zₐ * cmH2O * A4)  # eq. 13 McCullough & Porter
    λτ = ((float(λτR) + λτA + λτO) * m_Zₐ) + λτW  # eq. 14 McCullough & Porter
    λτ = min(λτ, MAX_OPTICAL_DEPTH)  # clamp to avoid numerical issues at low sun angles
    return (; rayleigh=λτR, total=λτ)
end

"""
    direct_irradiance(Sλ_n, ar², cz, λτ, λτR, m_Zₐ)

Calculate direct spectral irradiance and Rayleigh-only direct irradiance.

Returns `(; direct, rayleigh)` in W/m²/nm.
"""
function direct_irradiance(Sλ_n, ar², cz, λτ, λτR, m_Zₐ)
    part1 = Sλ_n * ar² * cz
    part2 = λτ > 0.0 ? exp(-λτ) : 0.0

    if part2 < MIN_IRRADIANCE
        Iλ = 0.0u"W/m^2/nm"
    else
        Iλ = ((ustrip(u"W/m^2/nm", part1) * part2) / 1000.0) * u"W/m^2/nm"
    end

    # Clamp to minimum value for numerical stability
    Iλ = max(Iλ, MIN_IRRADIANCE * u"W/m^2/nm")

    # Rayleigh-only direct irradiance
    Iᵣλ = (Sλ_n * ar² * cz) * exp(-float(λτR) * m_Zₐ) / 1000.0

    return (; direct=Iλ, rayleigh=Iᵣλ)
end

"""
    horizon_angle_at_azimuth(solar_azimuth, horizon_angles)

Look up the horizon angle at a given solar azimuth.
"""
function horizon_angle_at_azimuth(solar_azimuth, horizon_angles)
    azi = range(0u"°", stop=360u"°" - 360u"°" / length(horizon_angles), length=length(horizon_angles))
    return horizon_angles[argmin(abs.(solar_azimuth .- azi))]
end

"""
    trapezoidal_integrate!(∫vals, vals, λ, n)

Perform one step of trapezoidal integration for spectral values.
"""
function trapezoidal_integrate!(∫I, ∫Iᵣ, ∫D, ∫G, Iλ, Iᵣλ, Dλ, Gλ, λ, n)
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
        azimuth_angle = fill!(Vector{Union{Missing,typeof(0.0u"°")}}(undef, nsteps), 90.0u"°"),
        hour_angle_sunrise = fill(0.0, ndays),
        hour_solar_noon = fill(0.0, ndays),
        day_of_year = Vector{Int}(undef, nsteps),
        hour = Vector{Real}(undef, nsteps),
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
    allocate_buffers(nmax)

Allocate working buffers for solar radiation computation.

Returns a NamedTuple containing spectral integration arrays and
Chandrasekhar scattering buffers for reuse across multiple calls to `solar_radiation!`.
"""
function allocate_buffers(nmax)
    ∫G = fill(0.0u"W/m^2", nmax)
    ∫Iᵣ = fill(0.0u"W/m^2", nmax)
    ∫I = fill(0.0u"W/m^2", nmax)
    ∫D = fill(0.0u"W/m^2", nmax)
    Gλ = ∫G * u"1/nm"
    Iᵣλ = ∫G * u"1/nm"
    Iλ = ∫G * u"1/nm"
    Dλ = ∫G * u"1/nm"
    gamma = allocate_scattered_radiation()
    return (; ∫G, ∫Iᵣ, ∫I, ∫D, Gλ, Iᵣλ, Iλ, Dλ, gamma)
end

"""
    sunrise_hour_angle(declination, latitude)

Calculate sunrise/sunset hour angles (eq.7 McCullough & Porter 1971).

Returns `(; tanδ_tanϕ, H₊, H₋)`:
- `tanδ_tanϕ`: Product of tangents (used for polar day/night detection)
- `H₊`: Hour angle at sunset (radians)
- `H₋`: Hour angle at sunrise (hours)
"""
function sunrise_hour_angle(δ, ϕ)
    tanδ_tanϕ = -tan(δ) * tan(ϕ)
    H₊ = abs(tanδ_tanϕ) >= 1 ? π : abs(acos(tanδ_tanϕ))
    H₋ = 12.0 * H₊ / π
    return (; tanδ_tanϕ, H₊, H₋)
end

"""
    compute_spectral_irradiance!(buffers, params, sun_below_horizon)

Compute spectral irradiance for all wavelengths at a single timestep.

Modifies `buffers` in place with computed spectral values.
"""
function compute_spectral_irradiance!(buffers, params, sun_below_horizon)
    (; ∫G, ∫Iᵣ, ∫I, ∫D, Gλ, Iᵣλ, Iλ, Dλ) = buffers
    (; nmax, P, MR₀, τR, τO, τA, τW, Sλ, λ, ar², cz, intcz, m_Zₐ,
       ozone_depth, cmH2O, elevation_factors, scattered, scattered_uv,
       gamma_buffers, A, z, FD, FDQ, s̄) = params
    (; molecular, aerosol, ozone, water) = elevation_factors

    for n in 1:nmax
        τ = spectral_optical_depth(n, P, MR₀, τR, τO, τA, τW,
                                   molecular, aerosol, ozone, water,
                                   ozone_depth, m_Zₐ, cmH2O)

        direct = direct_irradiance(Sλ[n], ar², cz, τ.total, τ.rayleigh, m_Zₐ)
        Iλ[n], Iᵣλ[n] = direct.direct, direct.rayleigh

        if sun_below_horizon
            Iλ[n] = MIN_IRRADIANCE * u"W/m^2/nm"
            Iᵣλ[n] = MIN_IRRADIANCE * u"W/m^2/nm"
        end

        Dλ[n] = diffuse_irradiance(n, scattered, scattered_uv, τ.rayleigh, gamma_buffers,
                                   cz, intcz, Sλ[n], ar², m_Zₐ, A, z, FD, FDQ, s̄)
        Gλ[n] = Dλ[n] + Iλ[n]

        trapezoidal_integrate!(∫I, ∫Iᵣ, ∫D, ∫G, Iλ, Iᵣλ, Dλ, Gλ, λ, n)
    end
end

"""
    solar_radiation!(out, buffers, solar_model; kwargs...)

Mutating version of `solar_radiation` that writes results into pre-allocated buffers.

Use `allocate_output_arrays` and `allocate_buffers` to create the required buffers
for reuse across multiple calls.

See `solar_radiation` for argument documentation.
"""
function solar_radiation!(out, buffers, solar_model::SolarProblem;
    solar_terrain::SolarTerrain,
    days::Vector{<:Real}=[15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349],
    year::Real=1975,
    hours::AbstractVector{<:Real}=0:1:23,
    longitude_correction::Real=0.0,
)
    # Unpack model parameters with short aliases for equations
    (; solar_geometry_model, precipitable_water, scattered_uv, scattered, mixing_ratio_height,
       wavelength_count, wavelengths, ozone_column, rayleigh_optical_depth, ozone_optical_depth,
       aerosol_optical_depth, water_optical_depth, solar_spectral_irradiance,
       diffuse_sky_irradiance, diffuse_ground_reflected, single_scattering_albedo) = solar_model
    nmax, λ = wavelength_count, wavelengths
    τR, τO, τA, τW = rayleigh_optical_depth, ozone_optical_depth, aerosol_optical_depth, water_optical_depth
    Sλ, FD, FDQ, s̄ = solar_spectral_irradiance, diffuse_sky_irradiance, diffuse_ground_reflected, single_scattering_albedo
    cmH2O, MR₀ = precipitable_water, mixing_ratio_height

    # Unpack terrain
    (; horizon_angles, elevation, slope, aspect, albedo, atmospheric_pressure, latitude) = solar_terrain
    ϕ, P, A = latitude, atmospheric_pressure, albedo

    ndays, ntimes = length(days), length(hours)
    elevation_factors = elevation_correction(elevation)

    (; ∫G, ∫Iᵣ, ∫I, ∫D, Gλ, Iᵣλ, Iλ, Dλ, gamma) = buffers
    gamma_buffers = gamma

    step = 1
    H₋, tsn = 0.0, 0.0

    for i in 1:ndays
        for j in 1:ntimes
            d, t = days[i], hours[j]
            h, tsn = hour_angle(t, longitude_correction)
            solar_geom = solar_geometry(solar_geometry_model, ϕ; day_of_year=d, hour_angle=h)
            δ, z, ar² = solar_geom.solar_declination, solar_geom.zenith_angle, solar_geom.sun_distance_factor
            zsl = z

            # Twilight handling
            skylight = twilight_irradiance(z)
            if skylight !== nothing
                ∫D[nmax], ∫G[nmax] = skylight, skylight
                out.global_horizontal[step], out.diffuse_horizontal[step] = skylight, skylight
            end

            # Sunrise/sunset calculation
            sunrise = sunrise_hour_angle(δ, ϕ)
            H₋ = sunrise.H₋
            sun_up = is_sun_up(t - tsn, H₋)

            solar_azimuth = missing
            if sun_up || sunrise.tanδ_tanϕ == 1
                alt = (π / 2 - z)u"rad"
                solar_azimuth = solar_azimuth_angle(h, ϕ, δ)
                ahoriz = horizon_angle_at_azimuth(solar_azimuth, horizon_angles)

                # Slope geometry
                slope_geom = slope_zenith_angle(z, slope, solar_azimuth, aspect)
                zsl, czsl = slope_geom.zenith_angle, slope_geom.cosine_zenith

                # Refraction correction
                z = refraction_correction(z)
                cz = cos(z)
                intcz = floor(Int, 100.0 * cz + 1.0)
                m_Zₐ = optical_air_mass(z)
                ozone_depth = ozone_depth_lookup(ϕ, d, year, ozone_column)

                # Compute spectral irradiance
                params = (; nmax, P, MR₀, τR, τO, τA, τW, Sλ, λ, ar², cz, intcz, m_Zₐ,
                           ozone_depth, cmH2O, elevation_factors, scattered, scattered_uv,
                           gamma_buffers, A, z, FD, FDQ, s̄)
                compute_spectral_irradiance!(buffers, params, alt < ahoriz)

                # Store results
                out.global_spectra[step, :] .= Gλ
                out.rayleigh_spectra[step, :] .= Iᵣλ
                out.direct_spectra[step, :] .= Iλ
                out.diffuse_spectra[step, :] .= Dλ
                out.global_horizontal[step] = ∫G[nmax]
                out.global_terrain[step] = (slope > 0.0 && z < 90.0u"°") ? max(0.0u"W/m^2", (out.global_horizontal[step] / cz) * czsl) : out.global_horizontal[step]
                out.rayleigh_horizontal[step] = ∫Iᵣ[nmax]
                out.direct_horizontal[step] = ∫I[nmax]
                out.diffuse_horizontal[step] = ∫D[nmax]
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
    solar_radiation(solar_model; solar_terrain, days, year, hours, longitude_correction)

Compute solar radiation for a given model and terrain configuration.

Allocates output buffers internally. For repeated calls with the same dimensions,
use `solar_radiation!` with pre-allocated buffers for better performance.

# Arguments
- `solar_model::SolarProblem`: Solar radiation model parameters
- `solar_terrain::SolarTerrain`: Terrain configuration (elevation, slope, etc.)
- `days`: Vector of days of year (default: mid-month days)
- `year`: Year for leap year handling (default: 1975)
- `hours`: Hours of day to compute (default: 0:23)
- `longitude_correction`: Longitude correction in hours (default: 0.0)

# Returns
NamedTuple with zenith/azimuth angles, integrated irradiances, and spectral data.
"""
function solar_radiation(solar_model::SolarProblem;
    solar_terrain::SolarTerrain,
    days::Vector{<:Real}=[15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349],
    year::Real=1975,
    hours::AbstractVector{<:Real}=0:1:23,
    longitude_correction::Real=0.0,
)
    nmax = solar_model.wavelength_count
    ndays, ntimes = length(days), length(hours)
    nsteps = ndays * ntimes

    out = allocate_output_arrays(nsteps, ndays, nmax)
    buffers = allocate_buffers(nmax)

    return solar_radiation!(out, buffers, solar_model;
        solar_terrain, days, year, hours, longitude_correction)
end
