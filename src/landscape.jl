"""
    AbstractTerrain

Supertype of terrain descriptions passed to [`solar_radiation`](@ref) as `solar_terrain`.
"""
abstract type AbstractTerrain end

"""
    SolarTerrain(; elevation, horizon_angles, slope, aspect, albedo, atmospheric_pressure, latitude, longitude)

Site and terrain configuration for [`solar_radiation`](@ref). All angles and lengths are
`Unitful` quantities.

# Keywords
- `elevation`: height above sea level, sets the elevation correction of the optical depths.
- `horizon_angles`: horizon elevation angles, at equal azimuth steps clockwise from north
  (24 values give steps of 15°). Direct radiation is zero when the sun is below the horizon angle.
- `slope`: slope of the surface, from horizontal.
- `aspect`: azimuth the slope faces, clockwise from north.
- `albedo`: reflectance of the ground, from 0 to 1.
- `atmospheric_pressure`: pressure at the site, scales the Rayleigh optical depth. It can be calculated from the elevation
  with `atmospheric_pressure` of [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl).
- `latitude`: latitude, positive north.
- `longitude`: stored for reference, not used in the calculation. Solar time is set by
  the `longitude_correction` or `timezone_offset` keywords of `solar_radiation`.
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

"""
    AbstractSolarRadiation

Supertype of solar radiation models, such as [`SolarProblem`](@ref).
"""
abstract type AbstractSolarRadiation end

"""
    SpectralParams

Per-timestep parameters for spectral irradiance computation.

Bundles model constants, terrain properties, and per-step solar geometry
into a single object passed to `compute_spectral_irradiance!` and `diffuse_irradiance`.
"""
struct SpectralParams{
    DM<:AbstractDiffuseModel,
    Tp,Tmr,Ttr,Tto,Tta,Ttw,Tmg,
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
    mixed_gas_absorption::Tmg
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
    SolarProblem(; kw...)

Clear-sky solar radiation model for [`solar_radiation`](@ref), after McCullough & Porter (1971).

# Keywords
- `solar_geometry_model::AbstractSolarGeometryModel=McCulloughPorterSolarGeometry()`: orbit and declination model.
- `diffuse_model::AbstractDiffuseModel=DaveFurukawaScattering()`: one of [`DaveFurukawaScattering`](@ref),
  [`ChandrasekharScattering`](@ref), [`NoScattering`](@ref).
- `precipitable_water=1.0u"cm"`: precipitable water of the whole atmospheric column, a length
  (0.1 cm very dry, 1 cm moist, 2 cm humid tropical).
- `mixing_ratio_height=25.0u"km"`: sea-level meteorological range, the visibility at 0.55 μm.
- `wavelength_count=111`: number of wavelength intervals used.
- `wavelengths`: wavelength of each interval, in `nm`.
- `ozone_column`: total ozone in cm by latitude band (19 bands of 10°, from 90°S) and month, a 19×12 matrix.
- `rayleigh_optical_depth`, `ozone_optical_depth`, `aerosol_optical_depth`, `water_optical_depth`:
  sea-level vertical optical depths at each wavelength. The ozone depth applies at the reference
  column of 0.34 cm, and the water depth at 1 mm of precipitable water.
- `mixed_gas_absorption`: absorption coefficients of the uniformly mixed gases (O₂, CO₂) at each wavelength,
  after Bird & Riordan (1986).
- `solar_spectral_irradiance`: extraterrestrial solar spectrum at each wavelength. The values are
  stored ten times the tabulated ones, with a nominal unit of W m⁻² nm⁻¹, and the calculation divides
  by 1000 to give W m⁻² nm⁻¹.

The default tables are described in the manual, see `SolarRadiation.DEFAULT_WAVELENGTHS` and the other `DEFAULT_*` constants.
"""
@kwdef struct SolarProblem{SGM,DM<:AbstractDiffuseModel,PW,MRH,WC,WL,OC,ROD,OOD,AOD,WOD,MGA,SSI} <: AbstractSolarRadiation
    solar_geometry_model::SGM = McCulloughPorterSolarGeometry()
    diffuse_model::DM = DaveFurukawaScattering()
    precipitable_water::PW = 1.0u"cm" # precipitable water in air column 0.1 cm = very dry; 1 cm = moist air conditions; 2 cm = humid tropical conditions (note this is for the whole atmospheric profile not just near the ground)
    mixing_ratio_height::MRH = 25.0u"km" # mixing ratio height of the atmosphere
    wavelength_count::WC = 111 # Maximum number of wavelength intervals
    wavelengths::WL = DEFAULT_WAVELENGTHS # Vector of wavelength bins (e.g. in `nm`)
    ozone_column::OC = DEFAULT_OZONE_COLUMN # ozone column depth table indexed by latitude band and month (size 19×12)
    rayleigh_optical_depth::ROD = DEFAULT_RAYLEIGH_OPTICAL_DEPTH # vector of optical depths per wavelength for Rayleigh scattering
    ozone_optical_depth::OOD = DEFAULT_OZONE_OPTICAL_DEPTH # vector of optical depths per wavelength for ozone
    aerosol_optical_depth::AOD = DEFAULT_AEROSOL_OPTICAL_DEPTH # vector of optical depths per wavelength for aerosols
    water_optical_depth::WOD = DEFAULT_WATER_OPTICAL_DEPTH # vector of optical depths per wavelength for water vapor
    mixed_gas_absorption::MGA = DEFAULT_MIXED_GAS_ABSORPTION # vector of absorption coefficients per wavelength for O₂ and CO₂
    solar_spectral_irradiance::SSI = DEFAULT_SOLAR_SPECTRAL_IRRADIANCE # solar spectral irradiance per wavelength bin (e.g. in `mW * cm^-2 * nm^-1`)
end
