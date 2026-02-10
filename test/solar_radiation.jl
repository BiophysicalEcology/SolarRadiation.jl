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
names = [:lon, :lat, :EC, :elev, :slope, :aspect, :lamb, :IUV, :REFL, :P_atmos]
solarinput = (; zip(names, solarinput_vec)...)

# Conversions to unitful units
slope = (solarinput[:slope])*1.0u"°"
aspect = (solarinput[:aspect])*1.0u"°"
elevation = (solarinput[:elev])*1.0u"m"
horizon_angles = (DataFrame(CSV.File("$testdir/data/horizon.csv"))[:, 2])*1.0u"°"
albedo = solarinput[:REFL]*1.0
atmospheric_pressure = solarinput[:P_atmos]*1.0u"Pa"
scattered_uv = Bool(Int(solarinput[:IUV]))
latitude = solarinput[:lat]*1.0u"°"
longitude = solarinput[:lon]*1.0u"°"
#τA_nmr = (DataFrame(CSV.File("$testdir/data/TAI.csv"))[:, 2]*1.0)

hours = collect(0.0:1:23.0)
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]*1.0

solar_model = SolarProblem(; scattered_uv = scattered_uv)

solar_terrain = SolarTerrain(;
    horizon_angles, elevation, slope, aspect, albedo, atmospheric_pressure, latitude, longitude
)

@time solar_radiation_out = solar_radiation(solar_model; 
    solar_terrain,
    days,               # days of year
    hours,              # hours of day
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
λ = solar_model.wavelengths

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

@testset "dates interface" begin
    using Dates

    # Create dates matching the numeric days (mid-month days for 1975)
    year = 1975
    # The numeric days are day-of-year values, convert to dates
    dates_vec = [Date(year, 1, 1) + Day(Int(d) - 1) for d in days]

    # Run with dates interface
    solar_out_dt = solar_radiation(solar_model;
        solar_terrain,
        dates=dates_vec,
        hours=Hour(0):Hour(1):Hour(23),
    )

    # Compare zenith angles (apply same capping as line 52)
    zenith_dt = copy(solar_out_dt.zenith_angle)
    zenith_dt[zenith_dt .> 90u"°"] .= 90u"°"
    @test zenith_dt ≈ zenith_angle
    @test solar_out_dt.global_horizontal ≈ global_horizontal
    @test solar_out_dt.direct_spectra ≈ direct_spectra
end

@testset "dates with Date type" begin
    using Dates

    # Simple test with Date objects
    dates_vec = [Date(2024, 1, 15), Date(2024, 2, 15)]
    out = solar_radiation(solar_model;
        solar_terrain,
        dates=dates_vec,
        hours=Hour(0):Hour(6):Hour(18),  # 4 hours per day
    )

    @test length(out.zenith_angle) == 2 * 4  # 2 days × 4 hours
    @test out.day_of_year[1] == 15  # Jan 15
    @test out.day_of_year[5] == 46  # Feb 15
end

@testset "dates with DateTime type" begin
    using Dates

    # Test with DateTime objects
    dates_vec = [DateTime(2024, 6, 1), DateTime(2024, 6, 2), DateTime(2024, 6, 3)]
    out = solar_radiation(solar_model;
        solar_terrain,
        dates=dates_vec,
        hours=Hour(0):Hour(1):Hour(23),
    )

    @test length(out.zenith_angle) == 3 * 24  # 3 days × 24 hours
end

@testset "leap year via dates" begin
    using Dates

    # Leap year - 2024 has 366 days
    dates_leap = [Date(2024, 12, 31)]
    out_leap = solar_radiation(solar_model;
        solar_terrain,
        dates=dates_leap,
        hours=[Hour(12)],  # just noon
    )
    @test out_leap.day_of_year[1] == 366

    # Non-leap year - 2023 has 365 days
    dates_normal = [Date(2023, 12, 31)]
    out_normal = solar_radiation(solar_model;
        solar_terrain,
        dates=dates_normal,
        hours=[Hour(12)],
    )
    @test out_normal.day_of_year[1] == 365
end

@testset "hour offsets" begin
    using Dates

    dates_vec = [Date(2024, 6, 15)]

    # Test different hour offset patterns
    out1 = solar_radiation(solar_model;
        solar_terrain,
        dates=dates_vec,
        hours=Hour(0):Hour(1):Hour(23),  # hourly
    )
    @test length(out1.zenith_angle) == 24

    out2 = solar_radiation(solar_model;
        solar_terrain,
        dates=dates_vec,
        hours=Hour(0):Hour(3):Hour(21),  # every 3 hours
    )
    @test length(out2.zenith_angle) == 8

    out3 = solar_radiation(solar_model;
        solar_terrain,
        dates=dates_vec,
        hours=[Hour(6), Hour(12), Hour(18)],  # specific hours
    )
    @test length(out3.zenith_angle) == 3
end

@testset "orbital_angular_frequency" begin
    # Standard year
    @test orbital_angular_frequency(365) ≈ 2π / 365

    # Leap year
    @test orbital_angular_frequency(366) ≈ 2π / 366

    # 360-day calendar
    @test orbital_angular_frequency(360) ≈ 2π / 360

    # Default
    @test orbital_angular_frequency() ≈ 2π / 365
end
