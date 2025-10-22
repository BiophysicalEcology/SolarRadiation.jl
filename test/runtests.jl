using Aqua, SolarRadiation, DataFrames, CSV, Test, SafeTestsets, Unitful

@testset "Aqua.jl quality assurance" begin
    Aqua.test_all(SolarRadiation)
end

@safetestset "Test against NicheMapR outputs" begin include("solar_radiation.jl") end