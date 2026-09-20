# Ultraviolet radiation and sunburn

Sunburn (erythema) is caused by ultraviolet radiation, mostly in the UV-B, from 290 to 320 nm. Cripps and Ramsay (1970) irradiated the skin of the back of healthy Caucasian adults with narrow wavebands, and found the
minimal erythema dose (MED), the smallest energy that produced a just perceptible reddening of the skin at 24 hours, at each wavelength. The dose falls steeply with shorter wavelength,
from 1160 mJ cm⁻² at 320 nm to 6.2 mJ cm⁻² at 290 nm, so that short UV-B wavelengths are the most damaging, and the radiation there is strongly reduced by ozone. This tutorial uses these values to find the time to sunburn
from the spectrum of SolarRadiation.jl.

```@setup sunburn
using Main.FigureHelpers
```

## Time to sunburn

The rate at which the minimal dose of energy is accumulated is the spectral irradiance ``E_\lambda`` divided by the minimal energy ``ER_\lambda``,
integrated over the wavelengths, so the time to the minimal erythema is

```math
t = \left[\int_{290\,\mathrm{nm}}^{320\,\mathrm{nm}} \frac{E_\lambda}{ER_\lambda}\, d\lambda\right]^{-1}
```

The first seven wavelengths of the model, 5 nm apart from 290 to 320 nm, are those of the action spectrum. The energies follow the mean 24 hour MED of Table II of Cripps and Ramsay (1970) at nearby wavelengths, in mJ cm⁻². They are
converted to J m⁻².

```@example sunburn
using SolarRadiation, FluidProperties, Unitful
using CairoMakie

uv_wavelengths = SolarProblem().wavelengths[1:7]
minimal_erythemal_energy = uconvert.(u"J/m^2", [6.19, 6.86, 11.6, 25.1, 224.0, 560.0, 1160.0] .* u"mJ/cm^2") # Cripps and Ramsay (1970)

trapezoid(x, y) = sum((x[k+1] - x[k]) * (y[k+1] + y[k]) / 2 for k in 1:length(x)-1)

# `global_spectrum` is the global spectral irradiance at the wavelengths of the model, with units; the result is a time
function minutes_to_erythema(global_spectrum)
    rate = trapezoid(uv_wavelengths, global_spectrum[1:7] ./ minimal_erythemal_energy)
    rate > 0u"1/s" ? uconvert(u"minute", 1 / rate) : Inf * u"minute"
end

uv_wavelengths
```

The [`DaveFurukawaScattering`](@ref) diffuse model is used, which is the default. It was made for the ultraviolet, where it includes the absorption of the
diffuse radiation by ozone, which the model of Chandrasekhar does not, and it is fast. The pressure of a site is calculated from its elevation with `atmospheric_pressure` of [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl). The minutes to sunburn by time of day, at a site in Melbourne, Australia, in summer and winter:

```@example sunburn
model = SolarProblem()
terrain(latitude; elevation = 0.0u"m") = SolarTerrain(;
    elevation, albedo = 0.2, latitude = latitude * u"°", longitude = 0.0u"°",
    horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°", atmospheric_pressure = atmospheric_pressure(elevation),
)

hours = 6:0.25:18
out = solar_radiation(model; solar_terrain = terrain(-37.8), days = [15, 196], hours)
n = length(hours)

fig, ax = figure_axis("Solar time (h)", "Minutes to sunburn"; size = (700, 450), yscale = log10)
for (day, label) in ((1, "15 January"), (2, "15 July"))
    steps = (day - 1) * n + 1:day * n
    minutes = [minutes_to_erythema(out.global_spectra[s, :]) for s in steps]
    lines!(ax, hours[isfinite.(minutes)], ustrip.(u"minute", minutes[isfinite.(minutes)]); linewidth = 2, label)
end
axislegend(ax; position = :ct)
fig
```

The time to sunburn is shortest at solar noon, when the sun is highest. In the middle of summer it is
```@example sunburn
noon = findfirst(==(12.0), hours)
(january = minutes_to_erythema(out.global_spectra[noon, :]), july = minutes_to_erythema(out.global_spectra[n + noon, :]))
```
minutes at noon in January and July, with the clear sky and no shade, for the skin of the healthy Caucasian adults who were measured. The times
are longer for skin that is protected by pigment or by clothing, and depend on the skin type. The calculation also ignores the orientation of the skin
to the sun (it's horizontal surface irradiance) and it ignores reflection from the ground.

## Latitude and season

The ultraviolet radiation is greatest where the sun is high, and where the ozone column, which follows the
[table](../manual/atmosphere.md#Ozone) of the model, is small. The minutes to sunburn at solar noon, by latitude, on the two solstices and the equinox, are calculated next:

```@example sunburn
latitudes = -80:2:80
function noon_minutes(latitude, day; problem = model, elevation = 0.0u"m")
    result = solar_radiation(problem; solar_terrain = terrain(float(latitude); elevation), days = [day], hours = [12.0])
    minutes_to_erythema(result.global_spectra[1, :])
end

fig, ax = figure_axis("Latitude (°)", "Minutes to sunburn at noon"; size = (700, 450), yscale = log10)
for (day, label) in ((172, "June solstice"), (80, "March equinox"), (355, "December solstice"))
    minutes = [noon_minutes(latitude, day) for latitude in latitudes]
    lines!(ax, latitudes, ustrip.(u"minute", minutes); linewidth = 2, label)
end
ylims!(ax, 3, 300)
axislegend(ax; position = :lt)
fig
```

The time is shortest in the tropics and at the summer solstice in the hemisphere of the sun. It is longest at the poles in winter, when the sun is low or below the horizon.
The southern summer is slightly worse than the northern, because the earth is closest to the sun in January.

## Elevation and ozone

The ultraviolet radiation of the direct beam increases with elevation, because there is less air, and less ozone and aerosol, above a high site. The diffuse
ultraviolet radiation of the [`DaveFurukawaScattering`](@ref) model does not change with elevation, because its tables are for sea level, so the change in the
time to sunburn is driven by the direct beam only. The minutes to sunburn at noon in January at 37.8°S are calculated next:

```@example sunburn
using Markdown
heights = (0.0, 1000.0, 2000.0, 3000.0, 4500.0) .* u"m"

Markdown.parse(join(["| Elevation (m) | Minutes to sunburn |"; "| ---: | ---: |";
    ["| $(round(Int, height / u"m")) | $(round(noon_minutes(-37.8, 15; elevation = height) / u"minute"; digits = 1)) |" for height in heights]], Char(10)))
```

The ozone column of the model is a table of the ozone in cm by latitude and month, and it can be replaced. The effect of the depletion of
the ozone by a third, as in the ozone hole, in the same place and time. As with elevation only the direct beam responds, because the tables
of the diffuse ultraviolet radiation are for an ozone column of 0.34 cm, so the effect is understated:

```@example sunburn
depleted = SolarProblem(; ozone_column = SolarProblem().ozone_column .* (2 / 3))
(normal_ozone = noon_minutes(-37.8, 15), ozone_reduced_by_a_third = noon_minutes(-37.8, 15; problem = depleted))
```

## Sunburn over the globe

The radiation depends on the elevation and the aerosols of a place, as well as on its latitude and the season. The elevation of the land is read from the
CRU CL 2.0 dataset (New et al. 2002), and the aerosols from the [Global Aerosol Data Set](../manual/aerosols.md#The-Global-Aerosol-Data-Set) (GADS) at 50 % relative humidity,
with [RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl), as in the [Mapping tutorial](mapping.md), which has more on the data. The calculation is made for every cell of a 5° grid, at solar noon on the
two solstices, with the sea at sea level. It uses the aerosols of July for the June solstice, and of January for the December solstice.

```@example sunburn
using Rasters, RasterDataSources, NCDatasets

elv = read(RasterStack(getraster(CRUCL2); lazy = true).elv)
grid_lons = -177.5:5:177.5
grid_lats = -87.5:5:87.5
cell_heights = [max(coalesce(elv[X(Near(lon)), Y(Near(lat))], 0.0f0), 0.0f0) for lon in grid_lons, lat in grid_lats]

gads = read(Raster(getraster(GADS); name = :OPTDEPTH))
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

aerosols(lon, lat, season) = max.(interpolate_linear(gads_wavelengths,
    Float64.(collect(gads[X(Near(lon)), Y(Near(lat)), relhum = Near(50.0), season = At(season)])), model_wavelengths), 0.0)

nmax = model.wavelength_count
out_cell = allocate_output_arrays(1, 1, nmax)
buffers = allocate_buffers(nmax, model.diffuse_model)

function sunburn_map(day, season)
    result = fill(Inf * u"minute", length(grid_lons), length(grid_lats))
    for (j, lat) in enumerate(grid_lats), (i, lon) in enumerate(grid_lons)
        problem = SolarProblem(; aerosol_optical_depth = aerosols(lon, lat, season))
        height = Float64(cell_heights[i, j]) * u"m"
        solar_radiation!(out_cell, buffers, problem; solar_terrain = terrain(lat; elevation = height), days = [float(day)], hours = [12.0])
        result[i, j] = minutes_to_erythema(out_cell.global_spectra[1, :])
    end
    Raster(replace(ustrip.(u"minute", result), Inf => NaN), (X(grid_lons), Y(grid_lats)))
end

sunburn_map(172, 0.0) # compile
seconds = @elapsed june = sunburn_map(172, 0.0)
december = sunburn_map(355, 1.0)

limits = [5, 10, 20, 30, 60] # minutes, the upper limits of the classes
labels = ["<5", "5-10", "10-20", "20-30", "30-60"]
classes(minutes) = Raster(map(x -> (isnan(x) || x > last(limits)) ? NaN : Float64(searchsortedfirst(limits, x)), parent(minutes)), dims(minutes))
colours = cgrad(:viridis, length(limits), categorical = true)

f = Figure(size = (700, 750))
a1 = Axis(f[1, 1]; title = "Minutes to sunburn at noon, June solstice")
p1 = heatmap!(a1, classes(june); colormap = colours, colorrange = (0.5, length(limits) + 0.5), nan_color = :gray90)
countries!(a1)
Colorbar(f[1, 2], p1; ticks = (1:length(limits), labels))
a2 = Axis(f[2, 1]; title = "Minutes to sunburn at noon, December solstice")
p2 = heatmap!(a2, classes(december); colormap = colours, colorrange = (0.5, length(limits) + 0.5), nan_color = :gray90)
countries!(a2)
Colorbar(f[2, 2], p2; ticks = (1:length(limits), labels))
f
```

The colours are classes of the minutes to sunburn, up to 60 minutes, the time that is of interest for sunburn. The grey cells are those where the time is longer than an hour, which includes where the sun is below the horizon at noon.

The times are shortest in the tropics and subtropics in the summer of each hemisphere, and over high ground, the Andes in December and the Tibetan plateau in June,
where the elevation adds to the effect of the high sun. They become longer toward the winter pole. The differences with longitude are from the elevation and the aerosols of GADS.

## Speed

```@example sunburn
cells = length(cell_heights)
(cells = cells, seconds_per_map = round(seconds; sigdigits = 2), microseconds_per_cell = round(1e6 * seconds / cells; sigdigits = 2))
```

The UV is calculated with the [`DaveFurukawaScattering`](@ref) model at a single time, so the map takes a fraction of a second.

## References

Cripps DJ, Ramsay CA (1970) Ultraviolet action spectrum with a prism-grating monochromator. British Journal of Dermatology 82: 584-592.

New M, Lister D, Hulme M, Makin I (2002) A high-resolution data set of surface climate over global land areas. Climate Research 21: 1-25.
