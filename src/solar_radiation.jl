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
    (; solar_geometry_model, cmH2O, iuv, scattered, amr, nmax, Iλ, OZ, τR, τO, τA, τW, Sλ, FD, FDQ, s̄) = solar_model
    (; elevation, horizon_angles, slope, aspect, P_atmos, albedo) = solar_terrain

    ndays = length(days)    # number of days
    ntimes = length(hours)  # number of times
    nsteps = ndays * ntimes # total time steps

    # arrays to hold every time step's radiation between 300 and 320 nm in 2 nm steps
    GRλs = fill(0.0u"mW/nm/cm^2", nsteps, nmax) # wavelength-specific global radiation
    DRRλs = fill(0.0u"mW/nm/cm^2", nsteps, nmax)# wavelength-specific direct Rayleigh radiation
    DRλs = fill(0.0u"mW/nm/cm^2", nsteps, nmax) # wavelength-specific direct radiation
    SRλs = fill(0.0u"mW/nm/cm^2", nsteps, nmax) # wavelength-specific scattered radiation
    GRs = fill(0.0u"mW/cm^2", nsteps)           # total global radiation
    DRRs = fill(0.0u"mW/cm^2", nsteps)          # total direct Rayleigh radiation
    DRs = fill(0.0u"mW/cm^2", nsteps)           # total direct radiation
    SRs = fill(0.0u"mW/cm^2", nsteps)           # total scattered radiation
    gamma_buffers = allocate_scattered_uv()

    # arrays to hold zenith and azimuth angles each step
    Zs = fill(90.0u"°", nsteps)                 # zenith angles
    ZSLs = fill(90.0u"°", nsteps)               # slope zenith angles
    AZIs = Vector{Union{Missing,typeof(0.0u"°")}}(undef, nsteps)
    fill!(AZIs, 90.0u"°")   
    HHs = fill(0.0, ndays)                      # hour angles
    DOYs = Vector{Int}(undef, nsteps)           # day of year
    times = Vector{Real}(undef, nsteps)         # time
    step = 1
    HH = 0.0 # initialise sunrise hour angle
    for i in 1:ndays
        # arrays to hold radiation for a given hour between 300 and 320 nm in 2 nm steps
        GRINT = fill(0.0u"mW/cm^2", nmax)   # integrated global radiation component (direct + scattered)
        DRRINT = fill(0.0u"mW/cm^2", nmax)  # integrated direct Rayleigh radiation component
        DRINT = fill(0.0u"mW/cm^2", nmax)   # integrated direct radiation component
        SRINT = fill(0.0u"mW/cm^2", nmax)   # integrated scattered radiation component
        AIλ = fill(0.0u"nm", nmax)
        GRλ = GRINT * u"1/nm"               # wavelength-specific global radiation component (direct + scattered)
        DRRλ = GRINT * u"1/nm"              # wavelength-specific direct Rayleigh radiation component
        DRλ = GRINT * u"1/nm"               # wavelength-specific direct radiation component
        SRλ = GRINT * u"1/nm"               # wavelength-specific scattered radiation component
        alb = albedo#[i]
        for j in 1:ntimes
            d = days[i]
            t = hours[j]
            h, tsn = hour_angle(t, longitude_correction) # hour angle (radians)
            (; ζ, δ, z, AR2) = solar_geometry(solar_geometry_model, latitude; d, h) # compute ecliptic, declination, zenith angle and (a/r)^2
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
                SRINT[nmax] = skylight
                GRINT[nmax] = skylight
                GRs[step] = skylight
                SRs[step] = skylight
            end

            # testing cos(h) to see if it exceeds +1 or -1
            TDTL = -tan(δ) * tan(latitude) # from eq.7 McCullough & Porter 1971
            # Hour angle at sunrise/sunset (radians)
            # For polar day/night (|TDTL| ≥ 1), clamp to π
            H = abs(TDTL) >= 1 ? π : abs(acos(TDTL))
            # check if sunrise
            HH = 12.0 * H / π
            ts = t - tsn

            sun_up = true

            if ts <= 0.0 && abs(ts) > HH
                sun_up = false
            elseif ts > 0.0 && ts >= HH
                sun_up = false
            end

            if sun_up || TDTL == 1 # sun is up, proceed
                alt = (π / 2 - z)u"rad"
                altdeg = uconvert(u"°", alt).val
                # tazsun corresponds to tangent of azimuth
                tazsun = sin(h) / (cos(latitude) * tan(δ) - sin(latitude) * cos(h))
                # sun azimuth in radians
                azsun = atan(tazsun) * sign(latitude)
                # azimuth in degrees
                dazsun = uconvert(u"°", azsun)
                # correcting for hemisphere/quadrant
                dazsun = if h == 0.0 # Special case: hour angle = 0
                    180u"°"
                elseif h <= 0.0      # Morning - east of reference
                    dazsun <= 0u"°" ? -dazsun : 180u"°" - dazsun
                else                 # Afternoon - west of reference
                    dazsun < 0u"°" ? 180u"°" - dazsun : 360u"°" - dazsun
                end

                cz = cos(z)
                intcz = floor(Int, 100.0 * cz + 1.0)
                Z = uconvert(u"°", z)  # zenith angle in degrees

                # horizon angle - check this works when starting at 0 rather than e.g. 15 deg
                azi = range(0u"°", stop=360u"°" - 360u"°" / length(horizon_angles), length=length(horizon_angles))
                ahoriz = horizon_angles[argmin(abs.(dazsun .- azi))]

                # slope zenith angle calculation (Eq. 3.15 in Sellers 1965. Physical Climatology. U. Chicago Press)
                if slope > 0u"°"
                    czsl = cos(z) * cos(slope) + sin(z) * sin(slope) * cos(dazsun - aspect)
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
                    refr = 16.0 + ((z - 1.53589) * 15) / (π / 90)
                    refr = (refr / 60) * (π / 180)
                    z -= refr
                end

                # optical air mass (Rozenberg 1966 formula p.159 in book 'Twilight') ---
                airms = 1.0 / (cos(z) + (0.025 * exp(-11.0 * cos(z))))
                cz = cos(z)
                intcz = floor(Int, 100.0 * cz + 1.0)
                Z = uconvert(u"°", z)  # zenith angle in degrees

                # atmospheric ozone lookup
                # convert latitude in degrees to nearest 10-degree index
                llat = round(Int, (latitude + 100u"°") / 10u"°")
                # clamp llat index to valid range
                mon = month(Date(year, 1, 1) + Day(d - 1)) # month from day of year
                llat = clamp(llat, 1, size(OZ, 1))
                ozone = OZ[llat, mon]  # ozone thickness (cm) from lookup table

                (; molecular_corr, aerosol_corr, ozone_corr, water_vapour_corr) = elevation_correction(elevation)

                P = P_atmos

                for N in 1:nmax
                    τλ1 = (P / 101300u"Pa") * τR[N] * molecular_corr
                    τλ2 = (25.0u"km" / amr) * τA[N] * aerosol_corr
                    τλ3 = (ozone / 0.34) * τO[N] * ozone_corr
                    τλ4 = τW[N] * sqrt(airms * cmH2O * water_vapour_corr)
                    τλ = ((float(τλ1) + τλ2 + τλ3) * airms) + τλ4

                    if τλ > 80.0 # making sure that at low sun angles air mass doesn't make τλ too large
                        τλ = 80.0
                    end

                    part1 = Sλ[N] * AR2 * cz
                    part2 = τλ > 0.0 ? exp(-τλ) : 0.0
                    if part2 < 1.0e-24
                        DRλ[N] = 0.0u"mW / cm^2 / nm"
                    else
                        # TODO: ustrip to what
                        DRλ[N] = ((ustrip(part1) * part2) / 1000.0) * u"mW / cm^2 / nm"
                    end

                    # so the integrator doesn't get confused at very low sun angles
                    if DRλ[N] < 1.0e-25u"mW / cm^2 / nm"
                        DRλ[N] = 1.0e-25u"mW / cm^2 / nm"
                    end

                    DRRλ[N] = (Sλ[N] * AR2 * cz) * exp(-float(τλ1) * airms) / 1000.0

                    if altdeg < ahoriz
                        DRλ[N] = 1.0e-25u"mW / cm^2 / nm"
                        DRRλ[N] = 1.0e-25u"mW / cm^2 / nm"
                    end

                    # Sky (SRλ) and Global Radiation (GRλ)
                    if scattered == false
                        SRλ[N] = 0.0u"mW / cm^2 / nm"
                    elseif iuv
                        if τλ1 >= 0.03
                            GAMR, GAML, SBAR = scattered_uv!(gamma_buffers, τλ1)
                            SRλ[N] = (
                                         ((float(GAML[intcz]) + float(GAMR[intcz])) / (2.0 * (1.0 - alb * float(SBAR))))
                                         -
                                         exp(-float(τλ1) * airms)
                                     ) * cz * Sλ[N] * AR2 / 1000.0
                        else
                            SRλ[N] = 0.0u"mW / cm^2 / nm"
                        end
                    else
                        if N > 11
                            SRλ[N] = 0.0u"mW / cm^2 / nm"
                        else
                            # The option iuv = false has caused the program to enter this section which
                            # computes scattered radiation (SRλ) for 290 nm to 360 nm using a theory
                            # of radiation scattered from a Rayleigh (molecular) atmosphere with
                            # ozone absorption. The functions needed for the computation are stored
                            # as FD(N,I) and FDQ(N,I) where N is the wavelength index and I is
                            # (zenith angle + 5)/5 rounded off to the nearest integer value.
                            # The arrays FD and FDQ are for sea level (P = 1013 mb).
                            # TODO: ustrip to what
                            B = ustrip(u"°", Z) / 5
                            I = trunc(Int, B) + 1 + (B % 1 > 0.5)
                            FDAV = FD[N, I]
                            FDQDAV = FDQ[N, I]
                            SRλ[N] = (Sλ[N] / π) * (FDAV + FDQDAV * (alb / (1.0 - (alb * s̄[N])))) / 1000.0
                            SRλ[N] *= AR2
                        end
                    end

                    GRλ[N] = SRλ[N] + DRλ[N]

                    if N == 1
                        SRINT[1] = 0.0u"mW / cm^2"
                        DRRINT[1] = 0.0u"mW / cm^2"
                        DRINT[1] = 0.0u"mW / cm^2"
                        GRINT[1] = 0.0u"mW / cm^2"
                    else
                        AIλ[N] = Iλ[N]
                        AIλ[N-1] = Iλ[N-1]

                        Δλ = AIλ[N] - AIλ[N-1]

                        DRINT[N] = DRINT[N-1] + (Δλ * DRλ[N-1]) + (0.5 * Δλ * (DRλ[N] - DRλ[N-1]))
                        DRRINT[N] = DRRINT[N-1] + (Δλ * DRRλ[N-1]) + (0.5 * Δλ * (DRRλ[N] - DRRλ[N-1]))
                        SRINT[N] = SRINT[N-1] + (Δλ * SRλ[N-1]) + (0.5 * Δλ * (SRλ[N] - SRλ[N-1]))
                        GRINT[N] = GRINT[N-1] + (Δλ * GRλ[N-1]) + (0.5 * Δλ * (GRλ[N] - GRλ[N-1]))
                    end
                end
                GRλs[step, :] .= GRλ
                DRRλs[step, :] .= DRRλ
                DRλs[step, :] .= DRλ
                SRλs[step, :] .= SRλ
                GRs[step] = GRINT[nmax]
                DRRs[step] = DRRINT[nmax]
                DRs[step] = DRINT[nmax]
                SRs[step] = SRINT[nmax]
            else # sunrise, sunset or long day
                dazsun = missing
            end
            # Store into row `step`
            Zs[step] = Z
            ZSLs[step] = Zsl
            AZIs[step] = dazsun
            DOYs[step] = d
            times[step] = t
            step += 1
        end
        HHs[i] = HH     # save today's sunrise hour angle
    end

    return (
        zenith_angle = Zs,
        zenith_slope_angle = ZSLs,
        azimuth_angle = AZIs,
        hour_angle_sunrise = HHs,
        day_of_year = DOYs,
        hour = times,
        # TODO remove all this allocation from broadcasts
        # why is this conversion needed, what is the 10 about
        rayleigh_total = DRRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        direct_total = DRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        diffuse_total = SRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        global_total = GRs .* (10u"W/m^2" / 1u"mW/cm^2"),
        wavelength = Iλ,
        rayleigh_spectra = DRRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
        direct_spectra = DRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
        diffuse_spectra = SRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
        global_spectra = GRλs .* (10u"W/m^2" / 1u"mW/cm^2"),
    )
end