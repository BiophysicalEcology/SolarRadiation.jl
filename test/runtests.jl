using Aqua, SolarRadiation, DataFrames, CSV, Test, SafeTestsets, Unitful

@testset "Aqua.jl quality assurance" begin
    Aqua.test_all(SolarRadiation)
end

@safetestset "Test against NicheMapR outputs" begin include("solar_radiation.jl") end
@safetestset "Diffuse models" begin include("diffuse.jl") end
@safetestset "Latitudes" begin include("latitudes.jl") end
@safetestset "Reused output arrays" begin include("reuse.jl") end
@safetestset "Transmission" begin include("transmission.jl") end
