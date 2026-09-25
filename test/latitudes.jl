using SolarRadiation
using Unitful
using Test

terrain(latitude) = SolarTerrain(;
    horizon_angles=fill(0.0u"°", 24), elevation=0.0u"m", slope=0.0u"°", aspect=0.0u"°",
    albedo=0.2, atmospheric_pressure=101325.0u"Pa", latitude, longitude=0.0u"°",
)
hours = collect(0.0:1.0:23.0)
run(latitude, day; diffuse_model=DaveFurukawaScattering()) = solar_radiation(
    SolarProblem(; diffuse_model); solar_terrain=terrain(latitude), days=[day], hours,
)
global_w(out) = ustrip.(u"W/m^2", out.global_horizontal)

june_solstice, december_solstice, march_equinox, september_equinox = 172.0, 355.0, 80.0, 264.0

@testset "latitudes" begin
    @testset "sunrise hour angle" begin
        polar_night = SolarRadiation.sunrise_hour_angle(deg2rad(-23.44), 80.0u"°")
        @test polar_night.H₊ == 0.0
        @test polar_night.H₋ == 0.0
        polar_day = SolarRadiation.sunrise_hour_angle(deg2rad(23.44), 80.0u"°")
        @test polar_day.H₊ == π
        @test polar_day.H₋ == 12.0
        @test SolarRadiation.sunrise_hour_angle(0.0, 45.0u"°").H₋ ≈ 6.0
        @test !SolarRadiation.is_sun_up(0.0, 0.0)
    end

    @testset "finite and non-negative at all latitudes, all seasons, all diffuse models" begin
        for diffuse_model in (NoScattering(), DaveFurukawaScattering(), ChandrasekharScattering())
            for latitude in (90.0, 80.0, 70.0, 23.44, 0.0, -23.44, -70.0, -80.0, -90.0),
                day in (june_solstice, december_solstice, march_equinox, september_equinox)
                g = global_w(run(latitude * u"°", day; diffuse_model))
                @test all(isfinite, g)
                @test all(>=(0), g)
            end
        end
    end

    @testset "polar night and midnight sun" begin
        @test maximum(global_w(run(90.0u"°", december_solstice))) < 1e-3
        @test run(90.0u"°", december_solstice).hour_angle_sunrise[1] == 0.0
        @test maximum(global_w(run(-90.0u"°", june_solstice))) < 1e-3
        midnight_sun = run(90.0u"°", june_solstice)
        @test minimum(global_w(midnight_sun)) > 100
        @test midnight_sun.hour_angle_sunrise[1] == 12.0
        @test minimum(global_w(run(70.0u"°", june_solstice))) > 10
    end

    @testset "equator" begin
        for day in (march_equinox, september_equinox)
            @test run(0.0u"°", day).hour_angle_sunrise[1] ≈ 6.0 atol = 0.01
        end
        g = global_w(run(0.0u"°", march_equinox))
        @test argmax(g) == 13  # local solar noon at 12:00
        @test g[1] < 1e-3 && g[24] < 1e-3
    end

    @testset "hemispheres mirror each other half a year apart" begin
        for latitude in (45.0, 60.0)
            north = run(latitude * u"°", june_solstice).hour_angle_sunrise[1]
            south = run(-latitude * u"°", december_solstice).hour_angle_sunrise[1]
            @test north ≈ south atol = 0.01
        end
    end
end
