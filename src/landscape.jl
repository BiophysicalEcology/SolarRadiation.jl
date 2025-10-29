abstract type AbstractTerrain end

@kwdef struct SolarTerrain <: AbstractTerrain
    #elevation
    horizon_angles
    #slope
    #aspect
    #albedo
    #P_atmos
 end

 abstract type AbstractSolarRadiation end

"""
    SolarRadiation

# TODO who wrote this model what is it called

# Keyword Arguments

- `cmH2O::Real=1`: Precipitable water in cm for atmospheric column (e.g. 0.1: dry, 1.0: moist, 2.0: humid).
- `ϵ::Real=0.0167238`: Orbital eccentricity of Earth.
- `ω::Real=2π/365`: Mean angular orbital velocity of Earth (radians/day).
- `se::Real=0.39779`: Precomputed solar elevation constant.
- `d0::Real=80`: Reference day for declination calculations.
- `scattered_uv::Bool=false`: If `true`, uses the full scattered_uv model for diffuse radiation (expensive).
- `scattered::Bool=true`: If `true`, disables scattered light computations (faster).
- `MR₀::Quantity=25.0u"km"`: Mixing ratio height of the atmosphere.
- `nmax::Integer=111`: Maximum number of wavelength intervals.
- `λ::Vector{Quantity}`: Vector of wavelength bins (e.g. in `nm`).
- `ozone_column::Matrix{Float64}`: Ozone column depth table indexed by latitude band and month (size 19×12).
- `τR`, `τO`, `τA`, `τW`: Vectors of optical depths per wavelength for Rayleigh scattering, ozone, aerosols, and water vapor.
- `Sλ::Vector{Quantity}`: Solar spectral irradiance per wavelength bin (e.g. in `mW * cm^-2 * nm^-1`).
- `FD`, `FDQ`: Radiation scattered from the direct solar beam and reflected radiation
    rescattered downward as a function of wavelength, from tables in Dave & Furukawa (1966).
- `s̄`: a function of τR linked to molecular scattering in the UV range (< 360 nm)
"""
@kwdef struct SolarProblem <: AbstractSolarRadiation
    solar_geometry_model = McCulloughPorterSolarGeometry()
    cmH2O = 1.0 # precipitable cm H2O in air column 0.1 = very dry; 1 = moist air conditions; 2 = humid tropical conditions (note this is for the whole atmospheric profile not just near the ground)
    scattered_uv = false # if `true` uses the full scattered_uv model for diffuse radiation (expensive)
    scattered = true # if `false` disables scattered light computations (faster)
    MR₀ = 25.0u"km" # mixing ratio height of the atmosphere
    nmax = 111 # Maximum number of wavelength intervals
    # TODO better field names
    λ = DEFAULT_λ # Vector of wavelength bins (e.g. in `nm`)
    ozone_column = DEFAULT_OZONE_COLUMN # ozone column depth table indexed by latitude band and month (size 19×12)
    τR = DEFAULT_τR # vector of optical depths per wavelength for Rayleigh scattering
    τO = DEFAULT_τO # vector of optical depths per wavelength for ozone
    τA = DEFAULT_τA # vector of optical depths per wavelength for aerosols
    τW = DEFAULT_τW # vector of optical depths per wavelength for water vapor
    Sλ = DEFAULT_Sλ # solar spectral irradiance per wavelength bin (e.g. in `mW * cm^-2 * nm^-1`)
    FD = DEFAULT_FD # interpolated function of radiation scattered from the direct solar beam
    FDQ = DEFAULT_FDQ # interpolated function of radiation scattered from ground-reflected radiation
    s̄ = DEFAULT_s̄ # a function of τR linked to molecular scattering in the UV range (< 360 nm)
end