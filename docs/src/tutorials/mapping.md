# Mapping with GADS

The radiation of a place depends on its latitude, its elevation and its atmosphere. This tutorial makes global maps of
clear-sky radiation by calculating it for every cell of a grid, with the aerosols of the
[Global Aerosol Data Set](../manual/aerosols.md#The-Global-Aerosol-Data-Set) (GADS) at each cell. The data are read
with [RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl) and handled as rasters with
[Rasters.jl](https://github.com/rafaqz/Rasters.jl).

## Setup

Install the packages used in this tutorial:

```julia
using Pkg
Pkg.add(["Unitful", "Rasters", "RasterDataSources", "NCDatasets", "CairoMakie"])
Pkg.add(url = "https://github.com/BiophysicalEcology/SolarRadiation.jl")
Pkg.add(url = "https://github.com/BiophysicalEcology/FluidProperties.jl")
```

To download data you will need to specify a folder to put it in. You can do this by assigning the environment variable RASTERDATASOURCES_PATH:

```julia
ENV["RASTERDATASOURCES_PATH"] = joinpath(homedir(), "RasterDataSources") # or "/your/path/here"
```

## Acquiring the data

GADS is available through [RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl) as 
[`GADS`](https://github.com/EcoJulia/RasterDataSources.jl/blob/master/src/gads/gads.jl). It is a NetCDF file 
with the aerosol optical depth `OPTDEPTH` on a 5° grid of longitude and latitude, for 8 relative humidities, 
2 seasons and 25 wavelengths.

```@setup mapping
using Main.FigureHelpers
```

```@example mapping
using Rasters, RasterDataSources, NCDatasets
using SolarRadiation, FluidProperties, Unitful
using CairoMakie

gads = read(Raster(getraster(GADS); name = :OPTDEPTH))
```

The dimensions are `X` and `Y` for longitude and latitude, and the relative humidity `relhum` in %, the `season` (0 is July and
1 is January) and the `wavelength` in nm. A slice of the raster is chosen by naming the dimensions. The example below plots the
optical depth at 550 nm, at 50 % relative humidity, in July and in January. Note that the outlines of the countries on the maps 
are drawn with `countries!`, a helper of these docs that draws the [Natural Earth](https://www.naturalearthdata.com) outlines 
on an axis. Note also the definition of the `depth550` function at the beginning which is used further below. This way of 
specifying elements of a raster comes from the [DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) package on which [Rasters.jl](https://github.com/rafaqz/Rasters.jl) is built.

```@example mapping
depth550(season) = gads[relhum = Near(50.0), season = At(season), wavelength = Near(550.0)]

f = Figure(size = (700, 750))
a1 = Axis(f[1, 1]; title = "Aerosol optical depth at 550 nm, July")
p1 = heatmap!(a1, depth550(0.0); colorrange = (0, 0.6))
countries!(a1, color=(:white, 0.6))
Colorbar(f[1, 2], p1)
a2 = Axis(f[2, 1]; title = "Aerosol optical depth at 550 nm, January")
p2 = heatmap!(a2, depth550(1.0); colorrange = (0, 0.6))
countries!(a2, color=(:white, 0.6))
Colorbar(f[2, 2], p2)
f
```

The elevation of each cell sets its pressure and the elevation correction of the atmosphere. We use the 10 minute
CRU CL 2.0 dataset (New et al. 2002), which is available as `CRUCL2` in 
[RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl/blob/master/src/crucl2/crucl2.jl) and has land only. 
It is read at the centre of each GADS cell, and the sea is at sea level.

```@example mapping
elv = read(RasterStack(getraster(CRUCL2); lazy = true).elv)

template = depth550(0.0)
elevation = rebuild(template, [max(coalesce(elv[X(Near(lon)), Y(Near(lat))], 0.0f0), 0.0f0)
    for lon in lookup(template, X), lat in lookup(template, Y)])

f = Figure(size = (700, 380))
a1 = Axis(f[1, 1]; title = "Elevation of the cells (m)")
p1 = heatmap!(a1, elevation)
countries!(a1, color=(:white, 0.6))
Colorbar(f[1, 2], p1)
f
```

## Aerosols for the model

The model needs the optical depth at each of its 111 wavelengths, so the profile of each cell is
interpolated linearly from the 25 wavelengths of GADS:

```@example mapping
model = SolarProblem()
model_wavelengths = model.wavelengths
gads_wavelengths = collect(lookup(gads, :wavelength)) .* u"nm"

function interpolate_linear(x, y, xnew)
    map(xnew) do q
        q <= first(x) && return first(y)
        q >= last(x) && return last(y)
        k = searchsortedlast(x, q)
        t = (q - x[k]) / (x[k+1] - x[k])
        (1 - t) * y[k] + t * y[k+1]
    end
end

profile(i, j, season) = max.(interpolate_linear(gads_wavelengths,
    Float64.(collect(gads[X(i), Y(j), relhum = Near(50.0), season = At(season)])), model_wavelengths), 0.0)
length(profile(1, 1, 0.0))
```

## Radiation for every cell

For each cell of the grid we make the terrain, with the latitude of the cell and the elevation of the raster, and
a model with the aerosol profile of the cell, and calculate the clear-sky radiation at solar noon with
[`solar_radiation!`](@ref), reusing the arrays for every cell (see [Performance](../manual/performance.md)). The June solstice uses the aerosols of
July, and the December solstice those of January. The diffuse radiation is calculated with
[`ChandrasekharScattering`](@ref), which is needed to include all wavelengths (see [Diffuse models](../manual/diffuse_models.md)) but is slow, so
a single time is used. The pressure of each cell is calculated from its elevation with `atmospheric_pressure` of [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl). The results are the global irradiance and the diffuse fraction at noon.

```@example mapping
hours = [12.0]
model = SolarProblem(; diffuse_model = ChandrasekharScattering())
nmax = model.wavelength_count
out = allocate_output_arrays(length(hours), 1, nmax)
buffers = allocate_buffers(nmax, model.diffuse_model)

function solar_maps(day, season; default_aerosols = false)
    noon_global = zeros(size(template))
    noon_diffuse = zeros(size(template))
    for (j, lat) in enumerate(lookup(template, Y)), (i, lon) in enumerate(lookup(template, X))
        problem = default_aerosols ? model :
            SolarProblem(; diffuse_model = model.diffuse_model, aerosol_optical_depth = profile(i, j, season))
        height = Float64(elevation[i, j]) * u"m"
        terrain = SolarTerrain(;
            elevation = height, albedo = 0.2, latitude = lat * u"°", longitude = lon * u"°",
            horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°",
            atmospheric_pressure = atmospheric_pressure(height),
        )
        solar_radiation!(out, buffers, problem; solar_terrain = terrain, days = [float(day)], hours)
        noon_global[i, j] = ustrip(u"W/m^2", out.global_horizontal[1])
        noon_diffuse[i, j] = ustrip(u"W/m^2", out.diffuse_horizontal[1])
    end
    (noon_global = rebuild(template, noon_global),
        diffuse_fraction = rebuild(template, ifelse.(noon_global .> 1e-3, noon_diffuse ./ noon_global, NaN)))
end

june_solstice, december_solstice = 172, 355
solar_maps(june_solstice, 0.0) # compile
seconds = @elapsed june = solar_maps(june_solstice, 0.0)
december = solar_maps(december_solstice, 1.0)
nothing # hide
```

The global irradiance at solar noon on the solstices:

```@example mapping
f = Figure(size = (700, 750))
a1 = Axis(f[1, 1]; title = "Global irradiance at noon, June solstice (W m⁻²)")
p1 = heatmap!(a1, june.noon_global; colorrange = (0, 1200))
countries!(a1)
Colorbar(f[1, 2], p1)
a2 = Axis(f[2, 1]; title = "Global irradiance at noon, December solstice (W m⁻²)")
p2 = heatmap!(a2, december.noon_global; colorrange = (0, 1200))
countries!(a2)
Colorbar(f[2, 2], p2)
f
```

The sun is below the horizon at noon in the polar night, where the radiation is zero.

## Effect of the aerosols

The difference between the radiation with the aerosols of GADS and with the default profile, which is the same
everywhere (see [Aerosols](../manual/aerosols.md#The-default-profile)):

```@example mapping
june_default = solar_maps(june_solstice, 0.0; default_aerosols = true)
december_default = solar_maps(december_solstice, 1.0; default_aerosols = true)

f = Figure(size = (700, 750))
a1 = Axis(f[1, 1]; title = "GADS minus default aerosols, June solstice (W m⁻²)")
p1 = heatmap!(a1, june.noon_global .- june_default.noon_global; colorrange = (-150, 150), colormap = :balance)
countries!(a1)
Colorbar(f[1, 2], p1)
a2 = Axis(f[2, 1]; title = "GADS minus default aerosols, December solstice (W m⁻²)")
p2 = heatmap!(a2, december.noon_global .- december_default.noon_global; colorrange = (-150, 150), colormap = :balance)
countries!(a2)
Colorbar(f[2, 2], p2)
f
```

The diffuse fraction of the radiation at solar noon on the June solstice is:

```@example mapping
f = Figure(size = (700, 380))
a1 = Axis(f[1, 1]; title = "Diffuse fraction of the global irradiance at noon, June solstice")
p1 = heatmap!(a1, june.diffuse_fraction; colorrange = (0, 0.5))
countries!(a1)
Colorbar(f[1, 2], p1)
f
```

## Timing

Each map is a calculation of one time, at 111 wavelengths, with the [`ChandrasekharScattering`](@ref) diffuse model, for each cell of the grid.
The time taken gives an idea of the speed:

```@example mapping
ncells = length(template)
markdown_table(["Cells", "Seconds", "Milliseconds per cell"], [(ncells, round(seconds; sigdigits = 2), round(1000 * seconds / ncells; sigdigits = 2))])
```

The time is proportional to the number of cells and the number of times. With the default [`DaveFurukawaScattering`](@ref) it would
be more than a thousand times shorter.

## References

New M, Lister D, Hulme M, Makin I. 2002. A high-resolution data set of surface climate over global land areas. Climate Research 21: 1-25.

Koepke P, Hess M, Schult I, Shettle EP. 1997. Global Aerosol Data Set. Max-Planck-Institut für Meteorologie, Report No. 243, Hamburg.
