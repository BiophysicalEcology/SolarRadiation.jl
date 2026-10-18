# SolarRadiation

[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://biophysicalecology.github.io/SolarRadiation.jl/stable/)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://biophysicalecology.github.io/SolarRadiation.jl/dev/)
[![Docs Build](https://github.com/BiophysicalEcology/SolarRadiation.jl/actions/workflows/Documenter.yml/badge.svg)](https://github.com/BiophysicalEcology/SolarRadiation.jl/actions/workflows/Documenter.yml)
[![CI](https://github.com/BiophysicalEcology/SolarRadiation.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/BiophysicalEcology/SolarRadiation.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/BiophysicalEcology/SolarRadiation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BiophysicalEcology/SolarRadiation.jl/tree/main)

Clear-sky spectral solar radiation for biophysical ecology, after McCullough and Porter (1971): direct, diffuse and global
irradiance at 111 wavelengths from 290 to 4000 nm, for any latitude, elevation, slope, aspect and horizon, with units.

```julia
using SolarRadiation, Unitful

terrain = SolarTerrain(;
    elevation = 0.0u"m", horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°",
    albedo = 0.2, atmospheric_pressure = 101325.0u"Pa", latitude = -37.8u"°", longitude = 145.0u"°",
)
out = solar_radiation(SolarProblem(); solar_terrain = terrain, days = [15, 196], hours = 0:23)
out.global_horizontal
```

See the [documentation](https://biophysicalecology.github.io/SolarRadiation.jl/dev/) for the model, the diffuse
radiation models, terrain, and examples from the poles to the equator.
