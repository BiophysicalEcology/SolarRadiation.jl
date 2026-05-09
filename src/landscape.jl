abstract type AbstractDiffuseModel end

"""
    NoScattering <: AbstractDiffuseModel

Disables all diffuse (scattered) radiation. Equivalent to `scattered=false` in the original
McCullough & Porter (1971) model.
"""
struct NoScattering <: AbstractDiffuseModel end

"""
    DaveFurukawaScattering <: AbstractDiffuseModel

Diffuse UV irradiance from pre-tabulated lookup tables (Dave & Furukawa 1966).
Only applies to the first 11 wavelength intervals (UV range, < ~380 nm).
This is the default diffuse model.
"""
struct DaveFurukawaScattering <: AbstractDiffuseModel end

"""
    ChandrasekharScattering <: AbstractDiffuseModel

Full diffuse irradiance using Chandrasekhar's iterative X and Y functions.
Applies to all wavelengths. More accurate but significantly more expensive than
`DaveFurukawaScattering`. Equivalent to `scattered_uv=true` in the original model.
"""
struct ChandrasekharScattering <: AbstractDiffuseModel end

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
    Tod,Tch,Tef,Ta,Tz,
    Tfd,Tfq,Ts
}
    nmax::Int
    P::Tp
    MR₀::Tmr
    τR::Ttr
    τO::Tto
    τA::Tta
    τW::Ttw
    Sλ::Tsl
    λ::Tl
    ar²::Tar
    cz::Tcz
    intcz::Int
    m_Zₐ::Tmz
    ozone_depth::Tod
    cmH2O::Tch
    elevation_factors::Tef
    diffuse_model::DM
    A::Ta
    z::Tz
    FD::Tfd
    FDQ::Tfq
    s̄::Ts
end

"""
    SolarProblem

Solar radiation model parameters.

# Keyword Arguments

- `precipitable_water::Real=1`: Precipitable water in cm for atmospheric column (e.g. 0.1: dry, 1.0: moist, 2.0: humid).
- `diffuse_model::AbstractDiffuseModel=DaveFurukawaScattering()`: Diffuse radiation algorithm. Options:
  - `DaveFurukawaScattering()` (default): lookup-table method for UV wavelengths only (Dave & Furukawa 1966).
  - `ChandrasekharScattering()`: full iterative X/Y function method, all wavelengths (expensive).
  - `NoScattering()`: disables all diffuse radiation.
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
@kwdef struct SolarProblem{SGM,DM<:AbstractDiffuseModel,PW,SUV,SC,MRH,WC,WL,OC,ROD,OOD,AOD,WOD,SSI,DSI,DGR,SSA} <: AbstractSolarRadiation
    solar_geometry_model::SGM = McCulloughPorterSolarGeometry()
    diffuse_model::DM = DaveFurukawaScattering()
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
