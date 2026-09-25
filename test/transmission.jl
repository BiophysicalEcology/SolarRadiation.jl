using SolarRadiation
using SolarRadiation: direct_irradiance, DEFAULT_WAVELENGTHS, DEFAULT_WATER_OPTICAL_DEPTH, DEFAULT_MIXED_GAS_ABSORPTION,
    MAX_OPTICAL_DEPTH
using Unitful
using Test

@testset "direct irradiance" begin
    Sλ = 1500.0u"W/m^2/nm" # tabulated spectrum is stored ×1000
    # zero optical depth transmits the full beam
    I = direct_irradiance(Sλ, 1.0, 0.5, 0.0, 0.0, 1.0)
    @test I.direct ≈ 0.75u"W/m^2/nm"
    @test I.rayleigh ≈ 0.75u"W/m^2/nm"
    @test direct_irradiance(Sλ, 1.0, 0.5, 1.0, 0.0, 1.0).direct ≈ 0.75u"W/m^2/nm" * exp(-1.0)
end

@testset "strong water bands are opaque" begin
    terrain = SolarTerrain(;
        horizon_angles=fill(0.0u"°", 24), elevation=0.0u"m", slope=0.0u"°", aspect=0.0u"°",
        albedo=0.2, atmospheric_pressure=101325.0u"Pa", latitude=45.0u"°", longitude=0.0u"°",
    )
    radiation(τW) = solar_radiation(SolarProblem(; water_optical_depth=τW); solar_terrain=terrain, days=[172], hours=[8.0, 12.0])

    λ = ustrip.(u"nm", DEFAULT_WAVELENGTHS)
    opaque = findall(==(MAX_OPTICAL_DEPTH), DEFAULT_WATER_OPTICAL_DEPTH)
    @test λ[opaque] == [1400, 1420, 1900, 2600, 2700, 2800]

    out = radiation(DEFAULT_WATER_OPTICAL_DEPTH)
    @test all(<(1e-20u"W/m^2/nm"), out.direct_spectra[:, opaque])

    # integrated irradiance does not depend on the value chosen for opaque
    τW = Float64.(DEFAULT_WATER_OPTICAL_DEPTH)
    τW[opaque] .= MAX_OPTICAL_DEPTH / 2
    out_half = radiation(τW)
    @test out_half.direct_horizontal ≈ out.direct_horizontal rtol = 1e-12
    @test out_half.global_horizontal ≈ out.global_horizontal rtol = 1e-12

    # transparent band centres would add several percent to the direct beam
    τW[opaque] .= 0
    @test all(radiation(τW).direct_horizontal .> 1.03 .* out.direct_horizontal)
end

@testset "mixed gases do not scale with water vapour" begin
    terrain = SolarTerrain(;
        horizon_angles=fill(0.0u"°", 24), elevation=0.0u"m", slope=0.0u"°", aspect=0.0u"°",
        albedo=0.2, atmospheric_pressure=101325.0u"Pa", latitude=45.0u"°", longitude=0.0u"°",
    )
    radiation(w) = solar_radiation(SolarProblem(; precipitable_water=w); solar_terrain=terrain, days=[172], hours=[12.0])
    aM = DEFAULT_MIXED_GAS_ABSORPTION
    gas = findall(>(0), aM)
    @test ustrip.(u"nm", DEFAULT_WAVELENGTHS[gas]) == [1260, 1280, 2000, 2020, 2050]
    dry, humid = radiation(0.5u"cm"), radiation(3.0u"cm")
    @test dry.direct_spectra[:, gas] ≈ humid.direct_spectra[:, gas]
    no_gas = solar_radiation(SolarProblem(; precipitable_water=0.5u"cm", mixed_gas_absorption=zeros(length(aM)));
        solar_terrain=terrain, days=[172], hours=[12.0])
    @test all(dry.direct_spectra[:, gas] .< no_gas.direct_spectra[:, gas])
end
