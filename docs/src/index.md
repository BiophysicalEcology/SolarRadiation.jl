```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "SolarRadiation.jl"
  text: "Clear-sky solar radiation"
  tagline: "spectral direct, diffuse and global solar irradiance anywhere on earth, for biophysical ecology, with units."
  actions:
    - theme: brand
      text: Get Started
      link: /get_started
    - theme: alt
      text: View on Github
      link: https://github.com/BiophysicalEcology/SolarRadiation.jl
    - theme: alt
      text: API Reference
      link: /api

features:
  - title: ☀️ Sun position
    details: Solar <a class="highlight-link">declination, zenith and azimuth angles</a>, day length and the sun distance factor for any latitude, day and time, including polar day and night.
    link: /manual/solar_geometry
  - title: 🌈 Spectral irradiance
    details: <a class="highlight-link">Direct and diffuse irradiance</a> at 111 wavelengths from 290 to 4000 nm, after Rayleigh scattering, aerosols, ozone and water vapour, integrated to global irradiance.
    link: /manual/atmosphere
  - title: 🌤️ Diffuse models
    details: Choose no scattering, the tabulated UV scattering of <a class="highlight-link">Dave and Furukawa</a> or the full <a class="highlight-link">Chandrasekhar</a> solution, by passing a type.
    link: /manual/diffuse_models
  - title: ⛰️ Terrain
    details: Elevation, pressure, albedo, <a class="highlight-link">slope, aspect and horizon angles</a> are all part of the site description.
    link: /manual/terrain
  - title: 🌍 Poles to equator
    details: Midnight sun, polar night, the equatorial double peak and the tropics, from the same functions.
    link: /manual/latitude_examples
  - title: 🗺️ Maps
    details: Aerosol profiles from the <a class="highlight-link">Global Aerosol Data Set</a> and rasters from <a class="highlight-link">RasterDataSources.jl</a> give global maps of clear-sky radiation.
    link: /tutorials/mapping
  - title: 📏 Units and speed
    details: Every input and output uses <a class="highlight-link">Unitful.jl</a>. Repeated calls reuse buffers and are type-stable and nearly non-allocating.
    link: /manual/performance
---
```

## How to install SolarRadiation.jl?

SolarRadiation.jl can be installed from the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/BiophysicalEcology/SolarRadiation.jl")
```

## Model

The calculations follow McCullough and Porter (1971), as used in the microclimate model of
[NicheMapR](https://github.com/mrke/NicheMapR). See the [Introduction](manual/introduction.md) for the model and
how the package is organised.
