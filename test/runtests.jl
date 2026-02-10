using Aqua, SolarRadiation, DataFrames, CSV, Test, SafeTestsets, Unitful

@testset "Aqua.jl quality assurance" begin
    Aqua.test_all(SolarRadiation)
end

@safetestset "Theory-based tests" begin include("theory.jl") end

@safetestset "Test against NicheMapR outputs" begin include("solar_radiation.jl") end