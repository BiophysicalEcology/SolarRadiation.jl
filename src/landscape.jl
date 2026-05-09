abstract type AbstractTerrain end

@kwdef struct SolarTerrain{E,HA,S,As,Al,AP,La,Lo} <: AbstractTerrain
    elevation::E
    horizon_angles::HA
    slope::S
    aspect::As
    albedo::Al
    atmospheric_pressure::AP
    latitude::La
    longitude::Lo
 end

 abstract type AbstractSolarRadiation end

"""
    SolarProblem

Solar radiation model parameters.

# Keyword Arguments

- `precipitable_water::Real=1`: Precipitable water in cm for atmospheric column (e.g. 0.1: dry, 1.0: moist, 2.0: humid).
- `scattered_uv::Bool=false`: If `true`, uses the full scattered_uv model for diffuse radiation (expensive).
- `scattered::Bool=true`: If `false`, disables scattered light computations (faster).
- `mixing_ratio_height::Quantity=25.0u"km"`: Mixing ratio height of the atmosphere.
- `wavelength_count::Integer=111`: Maximum number of wavelength intervals.
- `wavelengths::Vector{Quantity}`: Vector of wavelength bins (e.g. in `nm`).
- `ozone_column::Matrix{Float64}`: Ozone column depth table indexed by latitude band and month (size 19×12).
- `rayleigh_optical_depth`, `ozone_optical_depth`, `aerosol_optical_depth`, `water_optical_depth`: Vectors of optical depths per wavelength.
- `solar_spectral_irradiance::Vector{Quantity}`: Solar spectral irradiance per wavelength bin.
- `diffuse_sky_irradiance`, `diffuse_ground_reflected`: Radiation scattered from the direct solar beam and reflected radiation
    rescattered downward as a function of wavelength, from tables in Dave & Furukawa (1966).
- `single_scattering_albedo`: Molecular scattering function in the UV range (< 360 nm).
"""
@kwdef struct SolarProblem{SGM,PW,SUV,SC,MRH,WC,WL,OC,ROD,OOD,AOD,WOD,SSI,DSI,DGR,SSA} <: AbstractSolarRadiation
    solar_geometry_model::SGM = McCulloughPorterSolarGeometry()
    precipitable_water::PW = 1.0 # precipitable cm H2O in air column 0.1 = very dry; 1 = moist air conditions; 2 = humid tropical conditions (note this is for the whole atmospheric profile not just near the ground)
    scattered_uv::SUV = false # if `true` uses the full scattered_uv model for diffuse radiation (expensive)
    scattered::SC = true # if `false` disables scattered light computations (faster)
    mixing_ratio_height::MRH = 25.0u"km" # mixing ratio height of the atmosphere
    wavelength_count::WC = 111 # Maximum number of wavelength intervals
    wavelengths::WL = DEFAULT_WAVELENGTHS # Vector of wavelength bins (e.g. in `nm`)
    ozone_column::OC = DEFAULT_OZONE_COLUMN # ozone column depth table indexed by latitude band and month (size 19×12)
    rayleigh_optical_depth::ROD = DEFAULT_RAYLEIGH_OPTICAL_DEPTH # vector of optical depths per wavelength for Rayleigh scattering
    ozone_optical_depth::OOD = DEFAULT_OZONE_OPTICAL_DEPTH # vector of optical depths per wavelength for ozone
    aerosol_optical_depth::AOD = DEFAULT_AEROSOL_OPTICAL_DEPTH # vector of optical depths per wavelength for aerosols
    water_optical_depth::WOD = DEFAULT_WATER_OPTICAL_DEPTH # vector of optical depths per wavelength for water vapor
    solar_spectral_irradiance::SSI = DEFAULT_SOLAR_SPECTRAL_IRRADIANCE # solar spectral irradiance per wavelength bin (e.g. in `mW * cm^-2 * nm^-1`)
    diffuse_sky_irradiance::DSI = DEFAULT_DIFFUSE_SKY_IRRADIANCE # interpolated function of radiation scattered from the direct solar beam
    diffuse_ground_reflected::DGR = DEFAULT_DIFFUSE_GROUND_REFLECTED # interpolated function of radiation scattered from ground-reflected radiation
    single_scattering_albedo::SSA = DEFAULT_SINGLE_SCATTERING_ALBEDO # a function of τR linked to molecular scattering in the UV range (< 360 nm)
end
