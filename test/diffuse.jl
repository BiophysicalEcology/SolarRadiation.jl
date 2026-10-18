using SolarRadiation
using Unitful
using Test

terrain = SolarTerrain(;
    horizon_angles=fill(0.0u"°", 24), elevation=300.0u"m", slope=0.0u"°", aspect=0.0u"°",
    albedo=0.2, atmospheric_pressure=97000.0u"Pa", latitude=-37.0u"°", longitude=145.0u"°",
)
days = [15.0, 196.0]
hours = collect(0.0:1.0:23.0)
run(model) = solar_radiation(model; solar_terrain=terrain, days, hours)

models = (NoScattering(), DaveFurukawaScattering(), ChandrasekharScattering())

@testset "diffuse models" begin
    @testset "no scattering" begin
        out = run(SolarProblem(; diffuse_model=NoScattering()))
        @test all(iszero, out.diffuse_spectra)
    end

    @testset "Dave-Furukawa is UV only" begin
        out = run(SolarProblem(; diffuse_model=DaveFurukawaScattering()))
        @test any(!iszero, out.diffuse_spectra[:, 1:11])
        @test all(iszero, out.diffuse_spectra[:, 12:end])
    end

    @testset "Dave-Furukawa agrees with Chandrasekhar where ozone absorbs little" begin
        wavelengths = ustrip.(u"nm", SolarProblem().wavelengths)
        uv = findall(w -> 330 <= w <= 360, wavelengths)
        noon = 13 # first day, 12:00
        dave_furukawa = run(SolarProblem(; diffuse_model=DaveFurukawaScattering()))
        chandrasekhar = run(SolarProblem(; diffuse_model=ChandrasekharScattering()))
        @test ustrip.(u"W/m^2/nm", dave_furukawa.diffuse_spectra[noon, uv]) ≈
            ustrip.(u"W/m^2/nm", chandrasekhar.diffuse_spectra[noon, uv]) rtol = 0.05
    end

    @testset "model owns its data" begin
        default = run(SolarProblem(; diffuse_model=DaveFurukawaScattering()))
        spherical = fill(0.5, length(SolarRadiation.DEFAULT_SPHERICAL_ALBEDO))
        custom = run(SolarProblem(; diffuse_model=DaveFurukawaScattering(; spherical_albedo=spherical)))
        @test custom.diffuse_spectra != default.diffuse_spectra
    end

    @testset "removed keywords are rejected" begin
        @test_throws MethodError SolarProblem(; scattered_uv=true)
        @test_throws MethodError SolarProblem(; scattered=false)
        @test_throws MethodError SolarProblem(; diffuse_sky_irradiance=nothing)
    end

    @testset "inference and allocation" begin
        for model in models
            solar_model = SolarProblem(; diffuse_model=model)
            nmax = solar_model.wavelength_count
            out = allocate_output_arrays(length(days) * length(hours), length(days), nmax)
            buffers = allocate_buffers(nmax, model)
            kw = (; solar_terrain=terrain, days, hours)
            @inferred solar_radiation!(out, buffers, solar_model; kw...)
            allocs = @allocated solar_radiation!(out, buffers, solar_model; kw...)
            @test allocs < 1024
        end
    end
end
