# Introduction

Solar radiation has a dominating effect on the microclimates available to organisms. Irradiance is the instantaneous amount
of solar radiation reaching the ground per unit area. It depends on

1. the radiation reaching the top of the atmosphere (extra-terrestrial radiation), which is set by the orbit of the earth 
   and is a function of the day of the year, the time of day and the latitude ([Solar geometry](solar_geometry.md));
2. the effect of the atmosphere as the radiation travels through it to the ground (terrestrial radiation), including 
   attenuation by Rayleigh scattering, aerosols, ozone, water vapour ([Atmosphere](atmosphere.md)) and associated scattering 
   of the direct solar beam into a diffuse component ([Diffuse models](diffuse_models.md));
3. the terrain: elevation, slope, aspect and hillshade ([Terrain](terrain.md)).

SolarRadiation.jl calculates these processes for cloudless skies. It follows McCullough and Porter (1971), with skylight 
before sunrise and after sunset from Rozenberg (1966) and including the effects of sloping surfaces and hills. It follows the 
solar radiation calculation of the microclimate model in [NicheMapR](https://github.com/mrke/NicheMapR) (Kearney and 
Porter 2017).

Cloud is not part of this package at present. Adjustments of the clear-sky irradiance for cloud cover, and the longwave 
radiation, are part of the microclimate model in [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl).

## Calculation

The total (global) irradiance on a horizontal surface is the sum of the direct beam ``I`` and the diffuse skylight ``D``,
which are calculated at each of 111 wavelengths ``\lambda`` from 290 to 4000 nm and then integrated:

```math
G = \int_0^\infty (I_\lambda + D_\lambda)\, d\lambda = I + D
```

[`solar_radiation`](@ref) does this for every day and time requested. At each time it

1. finds the position of the sun ([`solar_geometry`](@ref), [`hour_angle`](@ref)) and whether it is above the horizon,
2. for a sun above the horizon, corrects the zenith angle for refraction and finds the air mass,
3. calculates the optical depth of the atmosphere at every wavelength, and from it the direct irradiance,
4. calculates the diffuse irradiance with the chosen diffuse model,
5. integrates over the wavelengths with the trapezoidal rule, and
6. adjusts the global irradiance to the slope of the site.

Twilight, when the sun is at most 17° below the horizon, has a separate empirical calculation.

## Names and symbols

The names of the inputs and outputs of functions in the package are descriptive. Inside the functions the
symbols of the papers are used, so the code can be read alongside the equations. These symbols are used in the manual:

| Symbol | Name | Meaning |
| :----- | :--- | :------ |
| ``\phi`` | `latitude` | latitude of the site |
| ``\delta`` | `solar_declination` | latitude where the sun is overhead at noon |
| ``h`` | `hour_angle` | angle of the sun from solar noon, 15° per hour |
| ``H_-`` | `hour_angle_sunrise` | time of sunrise before solar noon, in hours |
| ``Z`` | `zenith_angle` | angle of the sun from the vertical |
| ``\cos Z`` | `cosine_zenith` | |
| ``(a/r)^2`` | `sun_distance_factor` | correction for the distance of the earth from the sun |
| ``m(Z_a)`` | `air_mass` | relative optical air mass, for the apparent zenith angle ``Z_a`` |
| ``S_\lambda`` | `solar_spectral_irradiance` | extraterrestrial solar spectrum |
| ``{}_\lambda\tau_R,\ {}_\lambda\tau_A,\ {}_\lambda\tau_O,\ {}_\lambda\tau_W`` | `rayleigh_optical_depth`, `aerosol_optical_depth`, `ozone_optical_depth`, `water_optical_depth` | vertical optical depths at wavelength ``\lambda`` |
| ``P`` | `atmospheric_pressure` | pressure at the site |
| ``MR_0`` | `mixing_ratio_height` | sea-level visibility |
| ``w`` | `precipitable_water` | precipitable water of the air column, in cm |
| ``A`` | `albedo` | reflectance of the ground |
| ``I_\lambda,\ D_\lambda,\ G_\lambda`` | `direct_spectrum`, `diffuse_spectrum`, `global_spectrum` | spectral irradiances |

## Interface

The choice of diffuse model is made by passing an instance of a type, as in the other packages of the group:

```julia
solar_radiation(SolarProblem(; diffuse_model = ChandrasekharScattering()); solar_terrain, days, hours)
```

The types, functions and their arguments are in the [API](../api.md).

## Related packages

Two packages of the [Climate Modeling Alliance](https://clima.caltech.edu) cover parts of the same calculation for
climate models:

- [Insolation.jl](https://github.com/CliMA/Insolation.jl) gives the position of the sun and the irradiance at the top of
  the atmosphere, from the date and time in UTC and the longitude. It includes the equation of time and, for past
  climates, orbital parameters from Laskar et al. (2004), but not the atmosphere or terrain. Its zenith and azimuth
  angles agree with [`solar_geometry`](@ref) to about 0.15°. SolarRadiation.jl uses local solar time, without the
  equation of time.
- [RRTMGP.jl](https://github.com/CliMA/RRTMGP.jl) calculates longwave and shortwave radiative transfer through the
  layers of an atmospheric column, with correlated-k gas absorption, clouds and aerosols (Pincus et al. 2019). It needs
  profiles of temperature, pressure and gas concentrations, takes the sun position from the host model, and gives
  broadband fluxes and heating rates. SolarRadiation.jl needs only column totals, and gives spectra and the effects of
  terrain. The water vapour optical depths were checked against RRTMGP.jl.

## References

McCullough EC, Porter WP (1971), Kearney MR, Porter WP (2017) and the other references are listed in the
[References](references.md).
