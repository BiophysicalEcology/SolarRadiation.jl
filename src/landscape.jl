abstract type AbstractTerrain end

"""
    SolarTerrain

Terrain configuration for solar radiation computation.
"""
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
    SpectralParams

Per-timestep parameters for spectral irradiance computation.

Bundles model constants, terrain properties, and per-step solar geometry
into a single object passed to `compute_spectral_irradiance!` and `diffuse_irradiance`.
"""
struct SpectralParams{
    DM<:AbstractDiffuseModel,
    Tp,Tmr,Ttr,Tto,Tta,Ttw,
    Tsl,Tl,Tar,Tcz,Tmz,
    Tod,Tch,Tef,Ta,Tz
}
    wavelength_count::Int
    atmospheric_pressure::Tp
    mixing_ratio_height::Tmr
    rayleigh_optical_depth::Ttr
    ozone_optical_depth::Tto
    aerosol_optical_depth::Tta
    water_optical_depth::Ttw
    solar_spectral_irradiance::Tsl
    wavelengths::Tl
    sun_distance_factor::Tar
    cosine_zenith::Tcz
    cosine_zenith_index::Int
    air_mass::Tmz
    ozone_depth::Tod
    precipitable_water::Tch
    elevation_factors::Tef
    diffuse_model::DM
    albedo::Ta
    zenith_angle::Tz
end

"""
    SolarProblem

Solar radiation model parameters.

# Keyword Arguments

- `precipitable_water::Real=1`: Precipitable water in cm for atmospheric column (e.g. 0.1: dry, 1.0: moist, 2.0: humid).
- `diffuse_model::AbstractDiffuseModel=DaveFurukawaScattering()`: diffuse radiation model,
  one of [`DaveFurukawaScattering`](@ref), [`ChandrasekharScattering`](@ref), [`NoScattering`](@ref).
- `mixing_ratio_height::Quantity=25.0u"km"`: Mixing ratio height of the atmosphere.
- `wavelength_count::Integer=111`: Maximum number of wavelength intervals.
- `wavelengths::Vector{Quantity}`: Vector of wavelength bins (e.g. in `nm`).
- `ozone_column::Matrix{Float64}`: Ozone column depth table indexed by latitude band and month (size 19×12).
- `rayleigh_optical_depth`, `ozone_optical_depth`, `aerosol_optical_depth`, `water_optical_depth`: Vectors of optical depths per wavelength.
- `solar_spectral_irradiance::Vector{Quantity}`: Solar spectral irradiance per wavelength bin.
"""
@kwdef struct SolarProblem{SGM,DM<:AbstractDiffuseModel,PW,MRH,WC,WL,OC,ROD,OOD,AOD,WOD,SSI} <: AbstractSolarRadiation
    solar_geometry_model::SGM = McCulloughPorterSolarGeometry()
    diffuse_model::DM = DaveFurukawaScattering()
    precipitable_water::PW = 1.0 # precipitable cm H2O in air column 0.1 = very dry; 1 = moist air conditions; 2 = humid tropical conditions (note this is for the whole atmospheric profile not just near the ground)
    mixing_ratio_height::MRH = 25.0u"km" # mixing ratio height of the atmosphere
    wavelength_count::WC = 111 # Maximum number of wavelength intervals
    wavelengths::WL = DEFAULT_WAVELENGTHS # Vector of wavelength bins (e.g. in `nm`)
    ozone_column::OC = DEFAULT_OZONE_COLUMN # ozone column depth table indexed by latitude band and month (size 19×12)
    rayleigh_optical_depth::ROD = DEFAULT_RAYLEIGH_OPTICAL_DEPTH # vector of optical depths per wavelength for Rayleigh scattering
    ozone_optical_depth::OOD = DEFAULT_OZONE_OPTICAL_DEPTH # vector of optical depths per wavelength for ozone
    aerosol_optical_depth::AOD = DEFAULT_AEROSOL_OPTICAL_DEPTH # vector of optical depths per wavelength for aerosols
    water_optical_depth::WOD = DEFAULT_WATER_OPTICAL_DEPTH # vector of optical depths per wavelength for water vapor
    solar_spectral_irradiance::SSI = DEFAULT_SOLAR_SPECTRAL_IRRADIANCE # solar spectral irradiance per wavelength bin (e.g. in `mW * cm^-2 * nm^-1`)
end
