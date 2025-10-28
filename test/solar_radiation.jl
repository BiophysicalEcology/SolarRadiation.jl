using SolarRadiation
using Unitful
using CSV, DataFrames
using Test

testdir = realpath(joinpath(dirname(pathof(SolarRadiation)), "../test"))

# NicheMapR simulation results
direct_spectra_nmr = Matrix(DataFrame(CSV.File("$testdir/data/drlam.csv"))[:, 4:114])
rayleigh_spectra_nmr = Matrix(DataFrame(CSV.File("$testdir/data/drrlam.csv"))[:, 4:114])
diffuse_spectra_nmr = Matrix(DataFrame(CSV.File("$testdir/data/srlam.csv"))[:, 4:114])
global_nmr = DataFrame(CSV.File("$testdir/data/global.csv"))[:, 2]
zenith_nmr = DataFrame(CSV.File("$testdir/data/zenith.csv"))[:, 2]
direct_spectra_nmr_units = direct_spectra_nmr*u"W/m^2/nm"
rayleigh_spectra_nmr_units = rayleigh_spectra_nmr*u"W/m^2/nm"
diffuse_spectra_nmr_units = diffuse_spectra_nmr*u"W/m^2/nm"

# NicheMapR simulation parameters
solarinput_vec = DataFrame(CSV.File("$testdir/data/input.csv"))[:, 2]
names = [
    :lon, :lat, :EC, :elev, :slope, :aspect, :lamb, :IUV, :REFL, :P_atmos
]
solarinput = (; zip(names, solarinput_vec)...)

slope = (solarinput[:slope])*1.0u"°"
aspect = (solarinput[:aspect])*1.0u"°"
elevation = (solarinput[:elev])*1.0u"m"
horizon_angles = (DataFrame(CSV.File("$testdir/data/horizon.csv"))[:, 2])*1.0u"°"
albedo = solarinput[:REFL]*1.0
P_atmos = solarinput[:P_atmos]*1.0u"Pa"
scattered_uv = Bool(Int(solarinput[:IUV]))
latitude = solarinput[:lat]*1.0u"°" # latitude
longitude =  solarinput[:lon]*1.0u"°" # longitude
#τA_nmr = (DataFrame(CSV.File("$testdir/data/TAI.csv"))[:, 2]*1.0)

hours = collect(0.0:1:23.0)
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]*1.0

solar_model = SolarProblem(; scattered_uv = scattered_uv)


solar_terrain = SolarTerrain(;
    horizon_angles,
)


@time solar_radiation_out = solar_radiation(solar_model, latitude, elevation, 
    slope, aspect, albedo, P_atmos;
    days,               # days of year
    hours,              # hours of day
    #latitude,           # latitude (degrees)
    solar_terrain,
);

zenith_angle = solar_radiation_out.zenith_angle
zenith_slope_angle = solar_radiation_out.zenith_slope_angle
zenith_angle[zenith_angle.>90u"°"] .= 90u"°"
azimuth_angle = solar_radiation_out.azimuth_angle

hour_angle_sunrise = solar_radiation_out.hour_angle_sunrise
global_horizontal = solar_radiation_out.global_horizontal
global_terrain = solar_radiation_out.global_terrain
direct_horizontal = solar_radiation_out.direct_horizontal
diffuse_horizontal = solar_radiation_out.diffuse_horizontal
rayleigh_horizontal = solar_radiation_out.rayleigh_horizontal
direct_spectra = solar_radiation_out.direct_spectra
diffuse_spectra = solar_radiation_out.diffuse_spectra
rayleigh_spectra = solar_radiation_out.rayleigh_spectra
λ = solar_radiation_out.wavelength

# Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
day_of_year = repeat(days, inner=length(hours))

# diffuse spectra test needs to be 1e-2 to pass with scattered_uv=true
@testset "solar radiation comparisons" begin
    @test ustrip.(u"°", zenith_angle) ≈ zenith_nmr rtol=1e-8
    @test ustrip.(u"W/m^2", global_horizontal) ≈ global_nmr rtol=1e-4
    @test direct_spectra ≈ direct_spectra_nmr_units rtol=1e-5
    @test diffuse_spectra ≈ diffuse_spectra_nmr_units rtol=1e-3
    @test rayleigh_spectra ≈ rayleigh_spectra_nmr_units rtol=1e-7
end
