# Get started

SolarRadiation.jl calculates clear-sky solar irradiance at the surface of the earth. Inputs and outputs are
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities.

```julia
using Pkg
Pkg.add(url = "https://github.com/BiophysicalEcology/SolarRadiation.jl")
```

## The site

A site is described by a [`SolarTerrain`](@ref). Below is an example of a flat site at sea level near 
Melbourne, Australia:

```@example get_started
using SolarRadiation, Unitful

terrain = SolarTerrain(;
    elevation = 0.0u"m",
    horizon_angles = fill(0.0u"°", 24), # no hills; 24 directions, 15° apart from north
    slope = 0.0u"°",
    aspect = 0.0u"°",
    albedo = 0.2,
    atmospheric_pressure = 101325.0u"Pa",
    latitude = -37.8u"°",
    longitude = 145.0u"°",
)
```

## The model

A [`SolarProblem`](@ref) holds the atmosphere and the choice of models. The defaults are a moist atmosphere
with 1 cm of precipitable water, the standard aerosol profile and the [`DaveFurukawaScattering`](@ref) diffuse model. 
The latter model gives diffuse radiation for the ultraviolet only, which is a very small part of the skylight 
(see [Diffuse models](manual/diffuse_models.md)),
so here we choose [`ChandrasekharScattering`](@ref), which includes all wavelengths:

```@example get_started
model = SolarProblem(; diffuse_model = ChandrasekharScattering())
model.precipitable_water
```

## Calculating radiation

[`solar_radiation`](@ref) takes the model and the site, the days of the year and the hours of solar time. Here is the
middle of January (day 15) and July (day 196), every half hour:

```@example get_started
days = [15, 196]
hours = 0:0.5:23.5
out = solar_radiation(model; solar_terrain = terrain, days, hours)
keys(out)
```

The results are vectors with one value per time step, ordered by day and then hour. The irradiance
at 1 pm solar time in January is:

```@example get_started
i = findfirst(==(13.0), out.hour) # the first day is January
out.global_horizontal[i]
```

split into direct and diffuse components:

```@example get_started
(direct = out.direct_horizontal[i], diffuse = out.diffuse_horizontal[i])
```

The spectra are matrices of time steps by wavelengths, in W/m²/nm:

```@example get_started
size(out.direct_spectra)
```

## Daily course

```@example get_started
using CairoMakie
using Main.FigureHelpers

fig = Figure(size = (700, 450))
for (col, (title, day)) in enumerate(("15 January", "15 July") .=> 1:2)
    steps = (day - 1) * length(hours) + 1:day * length(hours)
    ax = Axis(fig[1, col]; title, xlabel = "Solar time (h)", ylabel = col == 1 ? "Irradiance (W m⁻²)" : "",
        limits = (nothing, (0, 1400)))
    lines!(ax, out.hour[steps], ustrip.(out.global_horizontal[steps]); linewidth = 2, label = "global")
    lines!(ax, out.hour[steps], ustrip.(out.direct_horizontal[steps]); linewidth = 2, label = "direct")
    lines!(ax, out.hour[steps], ustrip.(out.diffuse_horizontal[steps]); linewidth = 2, label = "diffuse")
    col == 1 && axislegend(ax; position = :lt)
end
fig
```

## Next steps

- The [Introduction](manual/introduction.md) describes the model and how the code is organised.
- [Diffuse models](manual/diffuse_models.md) shows how to choose the scattering model.
- [Poles, tropics and equator](manual/latitude_examples.md) applies the model at contrasting latitudes.
- [Performance](manual/performance.md) shows how to run many calculations quickly.
