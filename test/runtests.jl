using Aqua, DataFrames, CSV, Test, SafeTestsets, Unitful, UnitfulMoles

@testset "Aqua.jl quality assurance" begin
    Aqua.test_all(SolarRadiation)
end

@safetestset "Test against NicheMapR outputs" begin include("solar_radiation.jl") end