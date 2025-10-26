function solar_radiation(model::SolarProblem, args...; kwargs...)
    any(ismissing, args) && return missing
    return solrad_core(model, args...; kwargs...)
end
function solar_radiation(solar_model::SolarProblem;
    days::Vector{<:Real}=[15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349],
    year::Real=1975, # to deal with leap years in obtaining month from day of year, default to non leap year
    latitude::Number,
    solar_terrain::SolarTerrain,
    longitude_correction::Real=0.0, # longitude correction, hours
    hours::AbstractVector{<:Real}=0:1:23,
)
    (; solar_geometry_model, cmH2O, scattered_uv, scattered, MR₀, nmax, λ, ozone_column, τR, τO, τA, τW, Sλ, FD, FDQ, s̄) = solar_model
    (; elevation, horizon_angles, slope, aspect, P_atmos, albedo) = solar_terrain
    ϕ = latitude
    ndays = length(days)    # number of days
    ntimes = length(hours)  # number of times
    nsteps = ndays * ntimes # total time steps

    # arrays to hold every time step's radiation between 300 and 320 nm in 2 nm steps
    λG = fill(0.0u"mW/nm/cm^2", nsteps, nmax) # wavelength-specific global radiation
    λIᵣ = fill(0.0u"mW/nm/cm^2", nsteps, nmax)# wavelength-specific direct Rayleigh radiation
    λI = fill(0.0u"mW/nm/cm^2", nsteps, nmax) # wavelength-specific direct radiation
    λD = fill(0.0u"mW/nm/cm^2", nsteps, nmax) # wavelength-specific scattered radiation
    G = fill(0.0u"mW/cm^2", nsteps)           # total global radiation
    Iᵣ = fill(0.0u"mW/cm^2", nsteps)          # total direct Rayleigh radiation
    I = fill(0.0u"mW/cm^2", nsteps)           # total direct radiation
    D = fill(0.0u"mW/cm^2", nsteps)           # total scattered radiation
    gamma_buffers = allocate_scattered_radiation()

    # arrays to hold zenith and azimuth angles each step
    zenith_angle = fill(90.0u"°", nsteps)       # zenith angles
    zenith_slope_angle = fill(90.0u"°", nsteps) # slope zenith angles
    azimuth_angle = Vector{Union{Missing,typeof(0.0u"°")}}(undef, nsteps)
    fill!(azimuth_angle, 90.0u"°")   
    hour_angle_sunrise = fill(0.0, ndays)       # hour angles
    day_of_year = Vector{Int}(undef, nsteps)    # day of year
    hour = Vector{Real}(undef, nsteps)          # time
    step = 1
    H₋ = 0.0 # initialise sunrise hour angle
    for i in 1:ndays
        # arrays to hold radiation for a given hour between 300 and 320 nm in 2 nm steps
        ∫G = fill(0.0u"mW/cm^2", nmax)   # integrated global radiation component (direct + scattered)
        ∫Iᵣ = fill(0.0u"mW/cm^2", nmax)  # integrated direct Rayleigh radiation component
        ∫I = fill(0.0u"mW/cm^2", nmax)   # integrated direct radiation component
        ∫D = fill(0.0u"mW/cm^2", nmax)   # integrated scattered radiation component
        Aλ = fill(0.0u"nm", nmax)
        Gλ = ∫G * u"1/nm"               # wavelength-specific global radiation component (direct + scattered)
        Iᵣλ = ∫G * u"1/nm"              # wavelength-specific direct Rayleigh radiation component
        Iλ = ∫G * u"1/nm"               # wavelength-specific direct radiation component
        Dλ = ∫G * u"1/nm"               # wavelength-specific scattered radiation component
        A = albedo#[i]
        for j in 1:ntimes
            d = days[i]
            t = hours[j]
            h, tsn = hour_angle(t, longitude_correction) # hour angle (radians)
            (; ζ, δ, z, ar²) = solar_geometry(solar_geometry_model, ϕ; d, h) # compute ecliptic, declination, zenith angle and (a/r)^2
            Z = uconvert(u"°", z)
            Zsl = Z
      
            # Compute twilight skylight irradiance (Rozenberg 1966; Diem 1966) # TODO add refs to doc, Rozenberg = Twilight. Plenum Press.
            if 88u"°" < Z < 107u"°"
                # Log10 of illuminance (lux) as a linear function of solar zenith angle
                log_illuminance = 41.34615384 - 0.423076923 * ustrip(u"°", Z) # p. 18,19 Rozenberg 1966
                
                # Convert lux → W/m² (via 1.46×10⁻³ kW/lumen)
                # p. 239 Documenta Geigy Scientific Tables. 1966. 6th ed. K. Diem, ed.
                skylight = (10.0^log_illuminance) * 1.46e-3u"mW/cm^2"

                # Assign twilight irradiance values
                ∫D[nmax] = skylight
                ∫G[nmax] = ∫D[nmax]
                G[step] = ∫G[nmax]
                D[step] = ∫D[nmax]
            end

            # testing cos(h) to see if it exceeds +1 or -1
            tanδ_tanϕ = -tan(δ) * tan(ϕ) # from eq.7 McCullough & Porter 1971
            # Hour angle at sunrise/sunset (radians)
            # For polar day/night (|TDTL| ≥ 1), clamp to π
            H₊ = abs(tanδ_tanϕ) >= 1 ? π : abs(acos(tanδ_tanϕ))

            # check if sunrise
            H₋ = 12.0 * H₊ / π # hour angle at sunrise
            ts = t - tsn

            sun_up = true

            if ts <= 0.0 && abs(ts) > H₋
                sun_up = false
            elseif ts > 0.0 && ts >= H₋
                sun_up = false
            end

            if sun_up || tanδ_tanϕ == 1 # sun is up, proceed
                alt = (π / 2 - z)u"rad"
                altdeg = uconvert(u"°", alt).val
                # tan_azimuth corresponds to tangent of azimuth
                tan_azimuth = sin(h) / (cos(ϕ) * tan(δ) - sin(ϕ) * cos(h))
                # sun azimuth in radians
                solar_azimuth = atan(tan_azimuth) * sign(latitude)
                # azimuth in degrees
                solar_azimuth_deg = uconvert(u"°", solar_azimuth) # TODO is this needed?
                # correcting for hemisphere/quadrant
                if h <= 0.0
                    # Morning - east of reference
                    if solar_azimuth_deg <= 0.0u"°"
                        # 1st Quadrant (0–90°)
                        solar_azimuth_deg = -1.0 * solar_azimuth_deg
                    else
                        # 2nd Quadrant (90–180°)
                        solar_azimuth_deg = 180.0u"°" - solar_azimuth_deg
                    end
                else
                    # Afternoon - west of reference
                    if solar_azimuth_deg < 0.0u"°"
                        # 3rd Quadrant (180–270°)
                        solar_azimuth_deg = 180.0u"°" - solar_azimuth_deg
                    else
                        # 4th Quadrant (270–360°)
                        solar_azimuth_deg = 360.0u"°" - solar_azimuth_deg
                    end
                end
                # Special case: hour angle = 0
                if h == 0.0
                    solar_azimuth_deg = 180.0u"°"
                end

                cz = cos(z)
                intcz = floor(Int, 100.0 * cz + 1.0)
                Z = uconvert(u"°", z)  # zenith angle in degrees

                # horizon angle - check this works when starting at 0 rather than e.g. 15 deg
                azi = range(0u"°", stop=360u"°" - 360u"°" / length(horizon_angles), length=length(horizon_angles))
                ahoriz = horizon_angles[argmin(abs.(solar_azimuth_deg .- azi))]

                # slope zenith angle calculation (Eq. 3.15 in Sellers 1965. Physical Climatology. U. Chicago Press)
                if slope > 0u"°"
                    czsl = cos(z) * cos(slope) + sin(z) * sin(slope) * cos(solar_azimuth_deg - aspect)
                    zsl = acos(czsl)
                    Zsl = min(uconvert(u"°", zsl), 90u"°") # cap at 90 degrees if sun is below slope horizon
                    intczsl = floor(Int, 100.0 * czsl + 1.0)
                else
                    czsl = cz
                    zsl = z
                    Zsl = Z
                    intczsl = intcz
                end

                # refraction correction check
                if z < 1.5358896
                    # skip refraction correction
                else
                    refraction = 16.0 + ((z - 1.53589) * 15) / (π / 90)
                    refraction = (refraction / 60) * (π / 180)
                    z -= refraction
                end

                # optical air mass (Rozenberg 1966 formula p.159 in book 'Twilight') ---
                m_Zₐ = 1.0 / (cos(z) + (0.025 * exp(-11.0 * cos(z))))
                cz = cos(z)
                intcz = floor(Int, 100.0 * cz + 1.0)
                Z = uconvert(u"°", z)  # zenith angle in degrees

                # atmospheric ozone lookup
                # convert latitude in degrees to nearest 10-degree index
                llat = round(Int, (ϕ + 100u"°") / 10u"°")
                # clamp llat index to valid range
                mon = month(Date(year, 1, 1) + Day(d - 1)) # month from day of year
                llat = clamp(llat, 1, size(ozone_column, 1))
                ozone_depth = ozone_column[llat, mon]  # ozone thickness (cm) from lookup table

                (; molecular, aerosol, ozone, water) = elevation_correction(elevation)
                A1 = molecular
                A2 = aerosol
                A3 = ozone
                A4 = water
                P = P_atmos

                for n in 1:nmax
                    λτR = (P / 101300u"Pa") * τR[n] * A1 # TODO 101300u"Pa" add as constant
                    λτA = (25.0u"km" / MR₀) * τA[n] * A2 # TODO add 25.0u"km" as constant
                    λτO = (ozone_depth / 0.34) * τO[n] * A3 # TODO add 0.34 as constant with units
                    λτW = τW[n] * sqrt(m_Zₐ * cmH2O * A4) # eq. 13 McCullough & Porter
                    λτ = ((float(λτR) + λτA + λτO) * m_Zₐ) + λτW # eq. 14 McCullough & Porter

                    if λτ > 80.0 # making sure that at low sun angles air mass doesn't make λτ too large
                        λτ = 80.0
                    end

                    part1 = Sλ[n] * ar² * cz
                    part2 = λτ > 0.0 ? exp(-λτ) : 0.0
                    if part2 < 1.0e-24
                        Iλ[n] = 0.0u"mW / cm^2 / nm"
                    else
                        # TODO: ustrip to what
                        Iλ[n] = ((ustrip(part1) * part2) / 1000.0) * u"mW / cm^2 / nm"
                    end

                    # so the integrator doesn't get confused at very low sun angles
                    if Iλ[n] < 1.0e-25u"mW / cm^2 / nm"
                        Iλ[n] = 1.0e-25u"mW / cm^2 / nm"
                    end

                    Iᵣλ[n] = (Sλ[n] * ar² * cz) * exp(-float(λτR) * m_Zₐ) / 1000.0 # TODO fix units

                    if altdeg < ahoriz
                        Iλ[n] = 1.0e-25u"mW / cm^2 / nm"
                        Iᵣλ[n] = 1.0e-25u"mW / cm^2 / nm"
                    end

                    # Sky (Dλ) and Global Radiation (Gλ)
                    if scattered == false
                        Dλ[n] = 0.0u"mW / cm^2 / nm"
                    elseif scattered_uv
                        if λτR >= 0.03
                            γᵣ, γₗ, s̄ = scattered_radiation!(gamma_buffers, λτR)
                            I₀_λ = cz * Sλ[n] * ar² / 1000.0 # TODO fix units
                            # eq. 15 McCullough & Porter 1971
                            Dλ[n] = (
                                         ((float(γₗ[intcz]) + float(γᵣ[intcz])) / (2.0 * (1.0 - A * float(s̄))))
                                         -
                                         exp(-float(λτR) * m_Zₐ)
                                     ) * I₀_λ
                        else
                            Dλ[n] = 0.0u"mW / cm^2 / nm"
                        end
                    else
                        if n > 11
                            Dλ[n] = 0.0u"mW / cm^2 / nm"
                        else
                            # The option scattered_uv = false has caused the program to enter this section which
                            # computes scattered radiation (Dλ) for 290 nm to 360 nm using a theory
                            # of radiation scattered from a Rayleigh (molecular) atmosphere with
                            # ozone absorption. The functions needed for the computation are stored
                            # as FD(n,k) and FDQ(n,k) where n is the wavelength index and k is
                            # (zenith angle + 5)/5 rounded off to the nearest integer value.
                            # The arrays FD and FDQ are for sea level (P = 1013 mb).
                            B = ustrip(u"°", Z) / 5
                            k = trunc(Int, B) + 1 + (B % 1 > 0.5)
                            flux_down = FD[n, k]
                            flux_down_div_Q = FDQ[n, k]
                            Q = (A / (1.0 - (A * s̄[n]))) # eq. 31 in Dave & Furukawa 1966
                            Dλ[n] = (Sλ[n] / π) * (flux_down + flux_down_div_Q * Q) / 1000.0
                            Dλ[n] *= ar²
                        end
                    end

                    Gλ[n] = Dλ[n] + Iλ[n]

                    if n == 1
                        ∫D[1] = 0.0u"mW / cm^2"
                        ∫Iᵣ[1] = 0.0u"mW / cm^2"
                        ∫I[1] = 0.0u"mW / cm^2"
                        ∫G[1] = 0.0u"mW / cm^2"
                    else
                        Aλ[n] = λ[n]
                        Aλ[n-1] = λ[n-1]

                        Δλ = Aλ[n] - Aλ[n-1]

                        ∫I[n] = ∫I[n-1] + (Δλ * Iλ[n-1]) + (0.5 * Δλ * (Iλ[n] - Iλ[n-1]))
                        ∫Iᵣ[n] = ∫Iᵣ[n-1] + (Δλ * Iᵣλ[n-1]) + (0.5 * Δλ * (Iᵣλ[n] - Iᵣλ[n-1]))
                        ∫D[n] = ∫D[n-1] + (Δλ * Dλ[n-1]) + (0.5 * Δλ * (Dλ[n] - Dλ[n-1]))
                        ∫G[n] = ∫G[n-1] + (Δλ * Gλ[n-1]) + (0.5 * Δλ * (Gλ[n] - Gλ[n-1]))
                    end
                end
                λG[step, :] .= Gλ
                λIᵣ[step, :] .= Iᵣλ
                λI[step, :] .= Iλ
                λD[step, :] .= Dλ
                G[step] = ∫G[nmax]
                Iᵣ[step] = ∫Iᵣ[nmax]
                I[step] = ∫I[nmax]
                D[step] = ∫D[nmax]
            else # sunrise, sunset or long day
                solar_azimuth_deg = missing
            end
            # Store into row `step`
            zenith_angle[step] = Z
            zenith_slope_angle[step] = Zsl
            azimuth_angle[step] = solar_azimuth_deg
            day_of_year[step] = d
            hour[step] = t
            step += 1
        end
    end

    return (
        zenith_angle,
        zenith_slope_angle,
        azimuth_angle,
        hour_angle_sunrise,
        day_of_year,
        hour,
        # TODO remove all this allocation from broadcasts
        # why is this conversion needed, what is the 10 about
        rayleigh_total = Iᵣ .* (10u"W/m^2" / 1u"mW/cm^2"),
        direct_total = I .* (10u"W/m^2" / 1u"mW/cm^2"),
        diffuse_total = D .* (10u"W/m^2" / 1u"mW/cm^2"),
        global_total = G .* (10u"W/m^2" / 1u"mW/cm^2"),
        wavelength = λ,
        rayleigh_spectra = λIᵣ .* (10u"W/m^2" / 1u"mW/cm^2"),
        direct_spectra = λI .* (10u"W/m^2" / 1u"mW/cm^2"),
        diffuse_spectra = λD .* (10u"W/m^2" / 1u"mW/cm^2"),
        global_spectra = λG .* (10u"W/m^2" / 1u"mW/cm^2"),
    )
end