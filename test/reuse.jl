using SolarRadiation
using Unitful
using Test

site(latitude) = SolarTerrain(;
    horizon_angles=fill(0.0u"°", 24), elevation=0.0u"m", slope=0.0u"°", aspect=0.0u"°",
    albedo=0.2, atmospheric_pressure=101325.0u"Pa", latitude, longitude=0.0u"°",
)

@testset "reused output arrays match fresh ones" begin
    model = SolarProblem()
    days = [172.0]
    hours = collect(0.0:6.0:18.0)
    nmax = model.wavelength_count
    out = allocate_output_arrays(length(hours), length(days), nmax)
    buffers = allocate_buffers(nmax, model.diffuse_model)

    solar_radiation!(out, buffers, model; solar_terrain=site(87.5u"°"), days, hours) # sun up all day
    solar_radiation!(out, buffers, model; solar_terrain=site(-87.5u"°"), days, hours) # sun down all day
    fresh = solar_radiation(model; solar_terrain=site(-87.5u"°"), days, hours)

    for name in (:rayleigh_horizontal, :direct_horizontal, :diffuse_horizontal, :global_horizontal, :global_terrain,
                 :rayleigh_spectra, :direct_spectra, :diffuse_spectra, :global_spectra)
        @test getproperty(out, name) == getproperty(fresh, name)
    end
    @test all(iszero, out.global_horizontal)
end
