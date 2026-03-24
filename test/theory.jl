using SolarRadiation, Test, Unitful
using SolarRadiation: hour_angle, solar_geometry, solar_azimuth_angle,
    McCulloughPorterSolarGeometry, optical_air_mass, refraction_correction,
    twilight_irradiance, elevation_correction, slope_zenith_angle,
    sunrise_hour_angle, SolarProblem, SolarTerrain

@testset "Solar Geometry" begin
    model = McCulloughPorterSolarGeometry()

    @testset "Hour angle" begin
        # At solar noon (t=12), hour angle should be 0
        h, tsn = hour_angle(12.0, 0.0)
        @test abs(h) < 1e-10u"rad"
        @test tsn == 12.0

        # Hour angle is linear in time: h = (π/12)(t - 12)
        h6, _ = hour_angle(6.0, 0.0)
        h18, _ = hour_angle(18.0, 0.0)
        @test abs(h6 - (-π/2) * u"rad") < 1e-10u"rad"
        @test abs(h18 - (π/2) * u"rad") < 1e-10u"rad"

        # Longitude correction shifts solar noon
        h, tsn = hour_angle(12.0, 1.0)  # 1 hour correction
        @test tsn == 13.0
        @test abs(h - (-π/12) * u"rad") < 1e-10u"rad"
    end

    @testset "Declination at equinoxes and solstices" begin
        # Vernal equinox (day ~80): δ ≈ 0°
        h, _ = hour_angle(12.0, 0.0)
        geom = solar_geometry(model, 0.0u"°"; day_of_year=80, hour_angle=h)
        @test abs(ustrip(u"°", geom.solar_declination)) < 1.0  # within 1° of 0

        # Summer solstice (day ~172): δ ≈ +23.44°
        geom = solar_geometry(model, 0.0u"°"; day_of_year=172, hour_angle=h)
        @test ustrip(u"°", geom.solar_declination) > 22.0
        @test ustrip(u"°", geom.solar_declination) < 24.5

        # Winter solstice (day ~355): δ ≈ -23.44°
        geom = solar_geometry(model, 0.0u"°"; day_of_year=355, hour_angle=h)
        @test ustrip(u"°", geom.solar_declination) < -22.0
        @test ustrip(u"°", geom.solar_declination) > -24.5
    end

    @testset "Declination bounds" begin
        h, _ = hour_angle(12.0, 0.0)
        for d in 1:365
            geom = solar_geometry(model, 0.0u"°"; day_of_year=d, hour_angle=h)
            δ_deg = abs(ustrip(u"°", geom.solar_declination))
            @test δ_deg ≤ 24.0  # declination bounded by ~23.5°
        end
    end

    @testset "Zenith angle at equator on equinox at noon" begin
        h, _ = hour_angle(12.0, 0.0)
        geom = solar_geometry(model, 0.0u"°"; day_of_year=80, hour_angle=h)
        @test ustrip(u"°", geom.zenith_angle) < 2.0  # nearly overhead
    end

    @testset "Zenith angle symmetric about solar noon" begin
        latitude = 45.0u"°"
        day = 172  # summer solstice

        # Test symmetry: Z(noon - Δt) ≈ Z(noon + Δt)
        for Δt in [1.0, 2.0, 3.0, 4.0]
            h_before, _ = hour_angle(12.0 - Δt, 0.0)
            h_after, _ = hour_angle(12.0 + Δt, 0.0)
            geom_before = solar_geometry(model, latitude; day_of_year=day, hour_angle=h_before)
            geom_after = solar_geometry(model, latitude; day_of_year=day, hour_angle=h_after)
            @test abs(geom_before.zenith_angle - geom_after.zenith_angle) < 0.01u"°"
        end
    end

    @testset "Zenith angle minimum at solar noon" begin
        latitude = 45.0u"°"
        day = 172

        h_noon, _ = hour_angle(12.0, 0.0)
        geom_noon = solar_geometry(model, latitude; day_of_year=day, hour_angle=h_noon)

        # Zenith should be larger at other times
        for t in [6.0, 8.0, 10.0, 14.0, 16.0, 18.0]
            h, _ = hour_angle(t, 0.0)
            geom = solar_geometry(model, latitude; day_of_year=day, hour_angle=h)
            @test geom.zenith_angle ≥ geom_noon.zenith_angle
        end
    end

    @testset "Sun-earth distance factor" begin
        h, _ = hour_angle(12.0, 0.0)

        # ar² should be close to 1.0 with small variations
        for d in 1:365
            geom = solar_geometry(model, 0.0u"°"; day_of_year=d, hour_angle=h)
            @test 0.96 < geom.sun_distance_factor < 1.04
        end

        # Perihelion (~day 3): ar² should be maximum (~1.034)
        geom_peri = solar_geometry(model, 0.0u"°"; day_of_year=3, hour_angle=h)
        # Aphelion (~day 186): ar² should be minimum (~0.967)
        geom_aph = solar_geometry(model, 0.0u"°"; day_of_year=186, hour_angle=h)
        @test geom_peri.sun_distance_factor > geom_aph.sun_distance_factor
    end

    @testset "Azimuth angle at solar noon" begin
        h, _ = hour_angle(12.0, 0.0)
        geom = solar_geometry(model, 45.0u"°"; day_of_year=172, hour_angle=h)
        δ = geom.solar_declination

        # At solar noon in NH, azimuth should be 180° (south)
        az = solar_azimuth_angle(h, 45.0u"°", δ)
        @test abs(az - 180.0u"°") < 0.1u"°"
    end

    @testset "Azimuth morning vs afternoon" begin
        geom = solar_geometry(model, 45.0u"°"; day_of_year=172, hour_angle=0.0u"rad")
        δ = geom.solar_declination

        # Morning: azimuth in [0°, 180°]
        h_morning, _ = hour_angle(8.0, 0.0)
        az_morning = solar_azimuth_angle(h_morning, 45.0u"°", δ)
        @test 0.0u"°" ≤ az_morning ≤ 180.0u"°"

        # Afternoon: azimuth in [180°, 360°]
        h_afternoon, _ = hour_angle(16.0, 0.0)
        az_afternoon = solar_azimuth_angle(h_afternoon, 45.0u"°", δ)
        @test 180.0u"°" ≤ az_afternoon ≤ 360.0u"°"
    end

    @testset "Sunrise hour angle" begin
        # At equinox, sunrise hour angle should give ~12h daylight
        h, _ = hour_angle(12.0, 0.0)
        geom = solar_geometry(model, 45.0u"°"; day_of_year=80, hour_angle=h)
        sunrise = sunrise_hour_angle(geom.solar_declination, 45.0u"°")
        daylight_hours = 2 * sunrise.H₋
        @test 11.5 < daylight_hours < 12.5

        # At equator on equinox, should be exactly 12h
        geom_eq = solar_geometry(model, 0.0u"°"; day_of_year=80, hour_angle=h)
        sunrise_eq = sunrise_hour_angle(geom_eq.solar_declination, 0.0u"°")
        @test sunrise_eq.H₋ ≈ 6.0 atol=0.5  # 6 hours from noon = 12h day
    end
end

@testset "Optical Air Mass" begin
    @testset "Reference values from Rozenberg (1966) Table 10" begin
        # Rozenberg formula should match tabulated values within ~5%
        # (tabulated values are from observations, formula is empirical fit)
        reference = [
            (0, 1.00),
            (30, 1.15),   # table shows 1.14-1.15
            (45, 1.41),   # table shows 1.41-1.45
            (60, 2.00),   # table shows 1.99-2.18
            (70, 2.90),   # table shows 2.90-3.30
            (75, 3.82),   # table shows 3.81-4.32
            (80, 5.60),   # table shows 5.56-6.03
        ]
        for (z_deg, m_ref) in reference
            z_rad = z_deg * π / 180
            m = optical_air_mass(z_rad)
            @test m ≈ m_ref rtol=0.05  # within 5%
        end
    end

    @testset "Air mass at zenith equals 1" begin
        # The Rozenberg formula gives m = 1/(1 + 0.025*exp(-11)) ≈ 0.9999959
        @test optical_air_mass(0.0) ≈ 1.0 rtol=1e-4
    end

    @testset "Air mass monotonically increases with zenith" begin
        z_values = 0:5:85
        m_prev = 0.0
        for z_deg in z_values
            z_rad = z_deg * π / 180
            m = optical_air_mass(z_rad)
            @test m > m_prev
            m_prev = m
        end
    end

    @testset "Air mass bounded" begin
        for z_deg in 0:88  # up to 88° (before extreme horizon)
            z_rad = z_deg * π / 180
            m = optical_air_mass(z_rad)
            @test 0.99 ≤ m ≤ 50.0  # 0.99 allows for Rozenberg formula at z=0
        end
    end
end

@testset "Refraction Correction" begin
    @testset "No correction below 88°" begin
        for z_deg in 0:5:85
            z_rad = z_deg * π / 180
            @test refraction_correction(z_rad) == z_rad
        end
    end

    @testset "Correction applied above 88°" begin
        z_rad = 89 * π / 180
        z_corrected = refraction_correction(z_rad)
        @test z_corrected < z_rad  # refraction reduces apparent zenith
    end
end

@testset "Twilight Irradiance" begin
    @testset "Returns nothing outside twilight range" begin
        @test twilight_irradiance(80u"°") === nothing
        @test twilight_irradiance(110u"°") === nothing
    end

    @testset "Returns value in twilight range" begin
        E = twilight_irradiance(95u"°")
        @test E !== nothing
        @test E > 0.0u"W/m^2"
    end

    @testset "Monotonically decreasing with zenith" begin
        E_90 = twilight_irradiance(90u"°")
        E_95 = twilight_irradiance(95u"°")
        E_100 = twilight_irradiance(100u"°")
        @test E_90 > E_95 > E_100
    end
end

@testset "Elevation Correction" begin
    @testset "Sea level correction factors" begin
        factors = elevation_correction(0.0u"m")
        @test factors.molecular ≈ 1.0 atol=0.02
        @test factors.aerosol ≈ 1.0 atol=0.02  # may not be exactly 1.0
        @test factors.ozone ≈ 1.0 atol=0.02
        @test factors.water ≈ 1.0 atol=0.02
    end

    @testset "Correction factors decrease with elevation" begin
        factors_0 = elevation_correction(0.0u"m")
        factors_1000 = elevation_correction(1000.0u"m")
        factors_3000 = elevation_correction(3000.0u"m")

        # Molecular (Rayleigh) decreases with elevation
        @test factors_1000.molecular < factors_0.molecular
        @test factors_3000.molecular < factors_1000.molecular

        # All factors should be positive
        @test factors_3000.molecular > 0
        @test factors_3000.aerosol > 0
        @test factors_3000.ozone > 0
        @test factors_3000.water > 0
    end
end

@testset "Slope Zenith Angle" begin
    @testset "Flat surface: Zsl = Z" begin
        z = 45.0u"°"
        slope = 0.0u"°"
        azimuth = 180.0u"°"
        aspect = 180.0u"°"

        result = slope_zenith_angle(z, slope, azimuth, aspect)
        @test abs(result.zenith_angle - z) < 0.01u"°"
        @test abs(result.cosine_zenith - cos(z)) < 1e-6
    end

    @testset "Slope zenith bounded to [0°, 90°]" begin
        for z_deg in 30:10:80
            for slope_deg in 10:10:40
                for aspect_offset in [0, 90, 180, 270]
                    z = z_deg * 1.0u"°"
                    slope = slope_deg * 1.0u"°"
                    azimuth = 180.0u"°"
                    aspect = (azimuth + aspect_offset * 1.0u"°")

                    result = slope_zenith_angle(z, slope, azimuth, aspect)
                    @test 0.0u"°" ≤ result.zenith_angle ≤ 90.0u"°"
                end
            end
        end
    end
end

@testset "Full Solar Radiation" begin
    terrain = SolarTerrain(;
        elevation=0.0u"m",
        horizon_angles=zeros(24) .* u"°",
        slope=0.0u"°",
        aspect=180.0u"°",
        albedo=0.2,
        atmospheric_pressure=101300.0u"Pa",
        latitude=45.0u"°",
        longitude=0.0u"°",
    )
    problem = SolarProblem()

    @testset "Irradiance positivity" begin
        result = solar_radiation(problem;
            solar_terrain=terrain,
            days=[172],
            hours=0:23
        )

        @test all(result.global_horizontal .≥ 0.0u"W/m^2")
        @test all(result.direct_horizontal .≥ 0.0u"W/m^2")
        @test all(result.diffuse_horizontal .≥ 0.0u"W/m^2")
        @test all(result.rayleigh_horizontal .≥ 0.0u"W/m^2")
    end

    @testset "Global = Direct + Diffuse" begin
        result = solar_radiation(problem;
            solar_terrain=terrain,
            days=[172],
            hours=6:18  # daylight hours
        )

        for i in eachindex(result.global_horizontal)
            G = result.global_horizontal[i]
            I = result.direct_horizontal[i]
            D = result.diffuse_horizontal[i]
            if G > 1.0u"W/m^2"  # skip near-zero values
                @test G ≈ I + D rtol=0.01
            end
        end
    end

    @testset "No NaN or Inf values" begin
        result = solar_radiation(problem;
            solar_terrain=terrain,
            days=[80, 172, 265, 355],
            hours=0:23
        )

        @test !any(isnan.(ustrip.(result.global_horizontal)))
        @test !any(isinf.(ustrip.(result.global_horizontal)))
        @test !any(isnan.(ustrip.(result.zenith_angle)))
        @test !any(isinf.(ustrip.(result.zenith_angle)))
    end

    @testset "Irradiance upper bound" begin
        result = solar_radiation(problem;
            solar_terrain=terrain,
            days=[172],
            hours=0:23
        )

        # Solar constant is ~1361 W/m², with distance factor max ~1.04
        # so max possible is ~1415 W/m²
        solar_max = 1420.0u"W/m^2"
        @test all(result.global_horizontal .≤ solar_max)
    end

    @testset "Daily zenith symmetry" begin
        result = solar_radiation(problem;
            solar_terrain=terrain,
            days=[172],
            hours=0:23
        )

        # Find solar noon (minimum zenith)
        noon_idx = argmin(result.zenith_angle)

        # Check symmetry around noon
        for offset in 1:5
            if noon_idx - offset ≥ 1 && noon_idx + offset ≤ length(result.zenith_angle)
                z_before = result.zenith_angle[noon_idx - offset]
                z_after = result.zenith_angle[noon_idx + offset]
                @test abs(z_before - z_after) < 0.5u"°"
            end
        end
    end

    @testset "Flat terrain: global_terrain = global_horizontal" begin
        result = solar_radiation(problem;
            solar_terrain=terrain,
            days=[172],
            hours=6:18
        )

        for i in eachindex(result.global_horizontal)
            if result.global_horizontal[i] > 1.0u"W/m^2"
                @test result.global_terrain[i] ≈ result.global_horizontal[i] rtol=0.001
            end
        end
    end
end

@testset "Spectral Properties" begin
    terrain = SolarTerrain(;
        elevation=0.0u"m",
        horizon_angles=zeros(24) .* u"°",
        slope=0.0u"°",
        aspect=180.0u"°",
        albedo=0.2,
        atmospheric_pressure=101300.0u"Pa",
        latitude=45.0u"°",
        longitude=0.0u"°",
    )
    problem = SolarProblem()

    result = solar_radiation(problem;
        solar_terrain=terrain,
        days=[172],
        hours=[12]  # solar noon
    )

    @testset "Spectral arrays have correct dimensions" begin
        nmax = problem.wavelength_count
        @test size(result.global_spectra, 2) == nmax
        @test size(result.direct_spectra, 2) == nmax
        @test size(result.diffuse_spectra, 2) == nmax
        @test size(result.rayleigh_spectra, 2) == nmax
    end

    @testset "Spectral irradiance positivity" begin
        @test all(result.global_spectra .≥ 0.0u"W/nm/m^2")
        @test all(result.direct_spectra .≥ 0.0u"W/nm/m^2")
    end
end
