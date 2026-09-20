# Elevation: the Himalaya

Elevation changes the radiation of a place through the atmosphere. There is less air above a high place, so the pressure and
the Rayleigh scattering are less, and the aerosols, which are mostly in the lower atmosphere, are fewer. This is included in
[`solar_radiation`](@ref) through the `elevation` and the `atmospheric_pressure` of the [`SolarTerrain`](@ref), see the
[Atmosphere](../manual/atmosphere.md#Elevation). This tutorial finds the effect in a region of great relief, the Himalaya of Nepal, with the elevation from
WorldClim (Fick and Hijmans 2017) and the aerosols of the [Global Aerosol Data Set](../manual/aerosols.md#The-Global-Aerosol-Data-Set) (GADS),
which are interpolated between its 5° cells.

## Setup

```julia
using Pkg
Pkg.add(["Unitful", "Rasters", "RasterDataSources", "ArchGDAL", "NCDatasets", "Extents", "CairoMakie"])
Pkg.add(url = "https://github.com/BiophysicalEcology/SolarRadiation.jl")
Pkg.add(url = "https://github.com/BiophysicalEcology/FluidProperties.jl")
```

Data are downloaded to the folder in the environment variable RASTERDATASOURCES_PATH:

```julia
ENV["RASTERDATASOURCES_PATH"] = joinpath(homedir(), "RasterDataSources") # or "/your/path/here"
```

## The elevation data

WorldClim has elevation for the globe at four resolutions: 10, 5 and 2.5 arc minutes, and 30 arc seconds (about 1 km). We
use 2.5 arc minutes (about 4.6 km), a file of 18 MB, because the 30 second file is 340 MB. The calculation below is the same for any resolution.
The area is 3° of longitude by 2° of latitude around Mount Everest (8849 m), from the Gangetic plain to the Tibetan plateau. The cells are averages of
the terrain of their area, so the highest is lower than the highest peak.

```@example elevation
using Rasters, RasterDataSources, ArchGDAL, NCDatasets, Extents
using SolarRadiation, FluidProperties, Unitful
using CairoMakie

area = Extent(X = (85.5, 88.5), Y = (26.5, 28.5))
elevation = read(crop(Raster(getraster(WorldClim{Elevation}, :elev; res = "2.5m")); to = area))
(cells = size(elevation), range_m = extrema(skipmissing(elevation)))
```

```@example elevation
f = Figure(size = (700, 450))
ax = Axis(f[1, 1]; title = "Elevation (m)", xlabel = "Longitude (°)", ylabel = "Latitude (°)")
hm = heatmap!(ax, elevation)
Colorbar(f[1, 2], hm)
f
```

## The effect of elevation alone

At one place the effect of elevation is found by changing the elevation and the pressure, which is given by `atmospheric_pressure` of
[FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl). The code below computes the radiation at noon on the December solstice at 28°N, with the default aerosols, and with [`ChandrasekharScattering`](@ref) for the diffuse radiation:

```@example elevation
model = SolarProblem(; diffuse_model = ChandrasekharScattering())

function terrain(latitude, longitude, height)
    SolarTerrain(;
        elevation = height, albedo = 0.2, latitude = latitude * u"°", longitude = longitude * u"°",
        horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°",
        atmospheric_pressure = atmospheric_pressure(height),
    )
end

heights = 0.0:250.0:8000.0
noon = map(heights) do h
    solar_radiation(model; solar_terrain = terrain(28.0, 87.0, h * u"m"), days = [355], hours = [12.0])
end

fig = Figure(size = (700, 450))
ax = Axis(fig[1, 1]; xlabel = "Elevation (m)", ylabel = "Irradiance at noon (W m⁻²)")
lines!(ax, heights, [ustrip(u"W/m^2", r.global_horizontal[1]) for r in noon]; linewidth = 2, label = "global")
lines!(ax, heights, [ustrip(u"W/m^2", r.direct_horizontal[1]) for r in noon]; linewidth = 2, label = "direct")
lines!(ax, heights, [ustrip(u"W/m^2", r.diffuse_horizontal[1]) for r in noon]; linewidth = 2, label = "diffuse")
axislegend(ax; position = :rc)
fig
```

The direct beam increases with elevation, and the diffuse radiation decreases, because there is less air to scatter the light.

## Aerosols between the GADS cells

The aerosols of GADS are for cells of 5°, so the optical depth at each cell of the elevation data is found by
interpolating it linearly between the four surrounding GADS points, for the wavelengths of the model, in the winter season of GADS
(January) at 50 % relative humidity:

```@example elevation
gads = read(Raster(getraster(GADS); name = :OPTDEPTH))
gads_lons, gads_lats = collect(lookup(gads, X)), collect(lookup(gads, Y))
gads_wavelengths = collect(lookup(gads, :wavelength)) .* u"nm"
gads_winter = Float64.(collect(gads[relhum = Near(50.0), season = At(1.0)])) # lon, lat, wavelength
model_wavelengths = model.wavelengths

function interpolate_linear(x, y, xnew)
    map(xnew) do q
        q <= first(x) && return first(y)
        q >= last(x) && return last(y)
        k = searchsortedlast(x, q)
        t = (q - x[k]) / (x[k+1] - x[k])
        (1 - t) * y[k] + t * y[k+1]
    end
end

function gads_profile(lon, lat)
    i = clamp(searchsortedlast(gads_lons, lon), 1, length(gads_lons) - 1)
    j = clamp(searchsortedlast(gads_lats, lat), 1, length(gads_lats) - 1)
    tx = (lon - gads_lons[i]) / (gads_lons[i+1] - gads_lons[i])
    ty = (lat - gads_lats[j]) / (gads_lats[j+1] - gads_lats[j])
    depth = @views (1 - tx) * (1 - ty) * gads_winter[i, j, :] + tx * (1 - ty) * gads_winter[i+1, j, :] +
        (1 - tx) * ty * gads_winter[i, j+1, :] + tx * ty * gads_winter[i+1, j+1, :]
    max.(interpolate_linear(gads_wavelengths, depth, model_wavelengths), 0.0)
end

aerosol_550 = [gads_profile(lon, lat)[findfirst(>=(550u"nm"), model_wavelengths)] for lon in lookup(elevation, X), lat in lookup(elevation, Y)]

f = Figure(size = (700, 450))
ax = Axis(f[1, 1]; title = "Aerosol optical depth at about 550 nm, January", xlabel = "Longitude (°)", ylabel = "Latitude (°)")
hm = heatmap!(ax, rebuild(elevation, aerosol_550))
Colorbar(f[1, 2], hm)
f
```

## Radiation for every cell

For each cell the radiation at solar noon on the December solstice is calculated three times, with [`solar_radiation!`](@ref) and reused arrays:
at sea level with the default aerosols, at the elevation of the cell with the default aerosols, and at the elevation of the cell with the aerosols of GADS.

```@example elevation
nmax = model.wavelength_count
out = allocate_output_arrays(1, 1, nmax)
buffers = allocate_buffers(nmax, model.diffuse_model)

function noon_radiation(; use_elevation, use_gads)
    result = zeros(size(elevation))
    for (j, lat) in enumerate(lookup(elevation, Y)), (i, lon) in enumerate(lookup(elevation, X))
        height = (use_elevation ? Float64(max(coalesce(elevation[i, j], 0.0f0), 0.0f0)) : 0.0) * u"m"
        problem = use_gads ? SolarProblem(; diffuse_model = model.diffuse_model, aerosol_optical_depth = gads_profile(lon, lat)) : model
        solar_radiation!(out, buffers, problem; solar_terrain = terrain(lat, lon, height), days = [355.0], hours = [12.0])
        result[i, j] = ustrip(u"W/m^2", out.global_horizontal[1])
    end
    rebuild(elevation, result)
end

noon_radiation(use_elevation = true, use_gads = true) # compile
seconds = @elapsed with_both = noon_radiation(use_elevation = true, use_gads = true)
sea_level = noon_radiation(use_elevation = false, use_gads = false)
with_elevation = noon_radiation(use_elevation = true, use_gads = false)
nothing # hide
```

```@example elevation
f = Figure(size = (700, 1000))
a1 = Axis(f[1, 1]; title = "Global irradiance at noon, elevation and GADS aerosols (W m⁻²)")
p1 = heatmap!(a1, with_both)
Colorbar(f[1, 2], p1)
a2 = Axis(f[2, 1]; title = "Effect of elevation: elevation minus sea level (W m⁻²)")
p2 = heatmap!(a2, with_elevation .- sea_level; colormap = :viridis)
Colorbar(f[2, 2], p2)
a3 = Axis(f[3, 1]; title = "Effect of the GADS aerosols: GADS minus default (W m⁻²)")
p3 = heatmap!(a3, with_both .- with_elevation; colormap = :balance, colorrange = (-30, 30))
Colorbar(f[3, 2], p3)
f
```

The elevation increases the radiation, especially on the Tibetan plateau. The GADS aerosols also increase it 
everywhere in the region in January, more in the mountains, because their optical depth is less than that of 
the default profile. The radiation of each cell against its elevation is plotted next:

```@example elevation
fig = Figure(size = (700, 450))
ax = Axis(fig[1, 1]; xlabel = "Elevation (m)", ylabel = "Global irradiance at noon (W m⁻²)")
heights_cells = vec(Float64.(max.(coalesce.(elevation, 0.0f0), 0.0f0)))
scatter!(ax, heights_cells, vec(with_both); markersize = 4, label = "elevation and GADS aerosols")
scatter!(ax, heights_cells, vec(sea_level); markersize = 4, label = "sea level, default aerosols")
axislegend(ax; position = :lt)
fig
```

## Speed

```@example elevation
cells = length(elevation)
(cells = cells, seconds_per_map = round(seconds; sigdigits = 2), milliseconds_per_cell = round(1000 * seconds / cells; sigdigits = 2))
```

The [`ChandrasekharScattering`](@ref) model is used in these examples to get full spectrum diffuse radiation, but this takes most of the time. With the default [`DaveFurukawaScattering`](@ref) it is more than a thousand times shorter, and the terrain of the [Saba tutorial](saba.md) shows a day of radiation for every cell of a much finer grid with the latter, faster scattering algorithm.

## References

Fick SE, Hijmans RJ (2017) WorldClim 2: new 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37: 4302-4315.
