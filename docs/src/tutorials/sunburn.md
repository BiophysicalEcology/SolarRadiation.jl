# Ultraviolet radiation and sunburn

Sunburn (erythema) is caused by ultraviolet radiation, mostly in the UV-B, from 290 to 320 nm. The erythema reference action spectrum, the relative effectiveness of each wavelength in
causing erythema, and the standard erythema dose (SED) are defined by the International Commission on Illumination and by ISO (CIE S 007/E-1998, ISO 17166:1999). This tutorial uses them to find
the time to sunburn from the spectrum of SolarRadiation.jl. The effectiveness falls by a factor of about 100 from 298 nm to 320 nm, so short UV-B wavelengths are the most damaging,
and the radiation there is strongly reduced by ozone.

```@setup sunburn
using Main.FigureHelpers
```

## Time to sunburn

The erythemal effective irradiance is the spectral irradiance ``E_\lambda`` weighted by the action spectrum ``s_{er}(\lambda)`` and integrated over the wavelengths,

```math
E_{er} = \int_{290\,\mathrm{nm}}^{400\,\mathrm{nm}} E_\lambda\, s_{er}(\lambda)\, d\lambda, \qquad
s_{er}(\lambda) = \begin{cases}
1 & \lambda \le 298\,\mathrm{nm} \\
10^{0.094\,(298 - \lambda)} & 298 < \lambda \le 328\,\mathrm{nm} \\
10^{0.015\,(140 - \lambda)} & 328 < \lambda \le 400\,\mathrm{nm}
\end{cases}
```

with ``\lambda`` in nm. One SED is an effective exposure of 100 J m⁻², so the time to ``n`` SED is ``100\,n / E_{er}`` in seconds. The minimal erythema dose (MED) depends on the person, and the standard
reserves the term for observations of people. It expects the MED of skin types I to IV to be from 1.5 to 6 SED. The time to sunburn here is the time to 1.5 SED, for the most sensitive of these skins,
or the time to 6 SED, for the least sensitive.

The 15 wavelengths of the model from 290 to 400 nm are used. The action spectrum is 1 at the first two, and falls by three orders of magnitude to 400 nm:

```@example sunburn
using SolarRadiation, FluidProperties, Unitful
using CairoMakie

# the erythema reference action spectrum of CIE S 007/E-1998, for a wavelength with units
function erythema_action_spectrum(wavelength)
    λ = wavelength / u"nm"
    λ <= 298 ? 1.0 : λ <= 328 ? 10^(0.094 * (298 - λ)) : λ <= 400 ? 10^(0.015 * (140 - λ)) : 0.0
end

standard_erythema_dose = 100.0u"J/m^2" # 1 SED

trapezoid(x, y) = sum((x[k+1] - x[k]) * (y[k+1] + y[k]) / 2 for k in 1:length(x)-1)

in_uv = findall(<=(400u"nm"), SolarProblem().wavelengths)
uv_wavelengths = SolarProblem().wavelengths[in_uv]
action_spectrum = erythema_action_spectrum.(uv_wavelengths)

# `global_spectrum` is the global spectral irradiance at the wavelengths of the model, with units
effective_irradiance(global_spectrum) = uconvert(u"W/m^2", trapezoid(uv_wavelengths, global_spectrum[in_uv] .* action_spectrum))

# the time to `sed` standard erythema doses
function minutes_to_erythema(global_spectrum, sed = 1.5)
    irradiance = effective_irradiance(global_spectrum)
    irradiance > 0u"W/m^2" ? uconvert(u"minute", sed * standard_erythema_dose / irradiance) : Inf * u"minute"
end

markdown_table(["Wavelength (nm)", "Erythema action spectrum"],
    [(w / u"nm", round(s; sigdigits = 3)) for (w, s) in zip(uv_wavelengths, action_spectrum)])
```

The [`DaveFurukawaScattering`](@ref) diffuse model is used, which is the default. It was made for the ultraviolet, where it includes the absorption of the
diffuse radiation by ozone, which the model of Chandrasekhar does not, and it is fast. It has no diffuse radiation beyond 360 nm, where the action spectrum is less than 0.001 of its maximum. The pressure of a site is calculated from its elevation with `atmospheric_pressure` of [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl). The minutes to sunburn by time of day, at a site in Melbourne, Australia, in summer and winter:

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
colours_of_days = Makie.wong_colors()
for (day, month) in ((1, "15 January"), (2, "15 July")), (sed, linestyle) in ((1.5, :solid), (6.0, :dash))
    steps = (day - 1) * n + 1:day * n
    minutes = [minutes_to_erythema(out.global_spectra[s, :], sed) for s in steps]
    lines!(ax, hours[isfinite.(minutes)], ustrip.(u"minute", minutes[isfinite.(minutes)]);
        linewidth = 2, linestyle, color = colours_of_days[day], label = "$month, $sed SED")
end
axislegend(ax; position = :ct)
fig
```

The time to sunburn is shortest at solar noon, when the sun is highest. The times to 1.5 SED and to 6 SED, the range for skin types I to IV, in the middle of summer and of winter, are:

```@example sunburn
noon = findfirst(==(12.0), hours)
markdown_table(["Day", "Minutes to 1.5 SED", "Minutes to 6 SED"],
    [("15 January", round(minutes_to_erythema(out.global_spectra[noon, :], 1.5) / u"minute"; digits = 1), round(minutes_to_erythema(out.global_spectra[noon, :], 6.0) / u"minute"; digits = 1)),
     ("15 July", round(minutes_to_erythema(out.global_spectra[n + noon, :], 1.5) / u"minute"; digits = 1), round(minutes_to_erythema(out.global_spectra[n + noon, :], 6.0) / u"minute"; digits = 1))])
```
at noon in January and July, with the clear sky and no shade. The times
are longer for skin that is protected by pigment or by clothing. The calculation also ignores the orientation of the skin
to the sun (it's horizontal surface irradiance) and it ignores reflection from the ground.

## Latitude and season

The ultraviolet radiation is greatest where the sun is high, and where the ozone column, which follows the
[table](../manual/atmosphere.md#Ozone) of the model, is small. The minutes to 1.5 SED at solar noon, by latitude, on the two solstices and the equinox, are calculated next:

```@example sunburn
latitudes = -80:2:80
function noon_minutes(latitude, day; problem = model, elevation = 0.0u"m")
    result = solar_radiation(problem; solar_terrain = terrain(float(latitude); elevation), days = [day], hours = [12.0])
    minutes_to_erythema(result.global_spectra[1, :])
end

fig, ax = figure_axis("Latitude (°)", "Minutes to 1.5 SED at noon"; size = (700, 450), yscale = log10)
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
time to sunburn is driven by the direct beam only. The minutes to 1.5 SED at noon in January at 37.8°S are calculated next:

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
markdown_table(["Ozone", "Minutes to 1.5 SED at noon"],
    [("normal", round(noon_minutes(-37.8, 15) / u"minute"; digits = 1)),
     ("reduced by a third", round(noon_minutes(-37.8, 15; problem = depleted) / u"minute"; digits = 1))])
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

limits = [6, 8, 10, 15, 20, 30, 60, 120, 600] # minutes, the upper limits of the classes
labels = ["<6", "6-8", "8-10", "10-15", "15-20", "20-30", "30-60", "60-120", "120-600", ">600"]
nclasses = length(limits) + 1 # the last class is longer than the last limit
classes(minutes) = Raster(map(x -> isnan(x) ? NaN : Float64(searchsortedfirst(limits, x)), parent(minutes)), dims(minutes))
colours = cgrad(:viridis, nclasses, categorical = true)

f = Figure(size = (700, 750))
a1 = Axis(f[1, 1]; title = "Minutes to 1.5 SED at noon, June solstice")
p1 = heatmap!(a1, classes(june); colormap = colours, colorrange = (0.5, nclasses + 0.5), nan_color = :gray90)
countries!(a1, color=(:white, 0.6))
Colorbar(f[1, 2], p1; ticks = (1:nclasses, labels))
a2 = Axis(f[2, 1]; title = "Minutes to 1.5 SED at noon, December solstice")
p2 = heatmap!(a2, classes(december); colormap = colours, colorrange = (0.5, nclasses + 0.5), nan_color = :gray90)
countries!(a2, color=(:white, 0.6))
Colorbar(f[2, 2], p2; ticks = (1:nclasses, labels))
f
```

The colours are classes of the minutes to 1.5 SED, with limits that give more detail where most cells are, from 6 to 30 minutes, and reach the long times of the winter hemisphere. The grey cells are those where the dose is not reached at noon, in the polar night.

The times are shortest in the tropics and subtropics in the summer of each hemisphere, and over high ground, the Andes in December and the Tibetan plateau in June,
where the elevation adds to the effect of the high sun. They become longer toward the winter pole. The differences with longitude are from the elevation and the aerosols of GADS.

## The highest risk over the year

The greatest risk of sunburn at a place is the shortest time to 1.5 SED over the year. The sun is highest at the summer solstice in the tropics of each hemisphere, but
at the equinoxes near the equator, so the shortest time of the four dates of the solstices and the equinoxes is found. The aerosols of July are used for the June
solstice and the September equinox, and those of January for the December solstice and the March equinox:

```@example sunburn
dates = ((80, 1.0), (172, 0.0), (266, 0.0), (355, 1.0)) # the day of the year and the GADS season
year_maps = [sunburn_map(day, season) for (day, season) in dates]
shortest = Raster(map((x...) -> minimum(isnan(v) ? Inf : v for v in x) |> (m -> isinf(m) ? NaN : m), (parent(m) for m in year_maps)...), dims(first(year_maps)))

f = Figure(size = (700, 400))
ax = Axis(f[1, 1]; title = "Shortest time to 1.5 SED over the year, at solar noon (minutes)")
hm = heatmap!(ax, classes(shortest); colormap = colours, colorrange = (0.5, nclasses + 0.5), nan_color = :gray90)
countries!(ax, color=(:white, 0.6))
Colorbar(f[1, 2], hm; ticks = (1:nclasses, labels))
f
```

The shortest time is from 6 to 8 minutes from about 35°S to 25°N, 8 to 10 minutes to about 45°S and 40°N, and it lengthens to 15 to 20 minutes at 65° of latitude and to 30 to 60 minutes
near the poles. The band of the shortest times reaches further from the equator in the southern hemisphere, because the earth is closest to the sun in January. Only the highest ground has less than
6 minutes, the Andes, the plateaus of eastern and southern Africa, the plateau of Mexico and the Tibetan plateau, where there is less air above the ground. This pattern can be compared with the geographic distribution of human skin pigmentation, which is darker where
the ultraviolet radiation is greater (Jablonski and Chaplin 2000). The pigmentation also depends on the history of the populations, and the map is for clear skies at noon, so the comparison is of the radiation
that a skin could receive.

## Comparison with Cripps and Ramsay

Cripps and Ramsay (1970) irradiated the skin of the back of healthy Caucasian adults with narrow wavebands and found the MED at each wavelength, from 1160 mJ cm⁻² at 320 nm to 6.2 mJ cm⁻² at 290 nm
(Table II, the mean at 24 hours). The time to the MED of the spectrum is that at which the sum over the wavelengths of the dose relative to the MED reaches 1,

```math
t = \left[\int_{290\,\mathrm{nm}}^{320\,\mathrm{nm}} \frac{E_\lambda}{\mathrm{MED}_\lambda}\, d\lambda\right]^{-1}
```

using the first seven wavelengths of the model, 5 nm apart:

```@example sunburn
cripps_wavelengths = uv_wavelengths[1:7]
cripps_med = uconvert.(u"J/m^2", [6.19, 6.86, 11.6, 25.1, 224.0, 560.0, 1160.0] .* u"mJ/cm^2")

function minutes_to_cripps_and_ramsay(global_spectrum)
    rate = trapezoid(cripps_wavelengths, global_spectrum[1:7] ./ cripps_med)
    rate > 0u"1/s" ? uconvert(u"minute", 1 / rate) : Inf * u"minute"
end

markdown_table(["Wavelength (nm)", "Cripps and Ramsay, relative to 290 nm", "CIE action spectrum"],
    [(w / u"nm", round(first(cripps_med) / m; sigdigits = 3), round(erythema_action_spectrum(w); sigdigits = 3)) for (w, m) in zip(cripps_wavelengths, cripps_med)])
```

The times to the MED of Cripps and Ramsay and to 1 SED and 1.5 SED at noon, and the effective dose of the CIE spectrum that the Cripps and Ramsay time gives:

```@example sunburn
comparison_cases = (("Melbourne, 15 January", -37.8, 15), ("Melbourne, 15 July", -37.8, 196), ("Equator, 21 March", 0.0, 80))
comparison_rows = map(comparison_cases) do (name, latitude, day)
    result = solar_radiation(model; solar_terrain = terrain(latitude), days = [day], hours = [12.0])
    spectrum = result.global_spectra[1, :]
    minutes = minutes_to_cripps_and_ramsay(spectrum)
    (name, round(minutes / u"minute"; digits = 1), round(minutes_to_erythema(spectrum, 1.0) / u"minute"; digits = 1),
        round(minutes_to_erythema(spectrum, 1.5) / u"minute"; digits = 1),
        round(uconvert(u"J/m^2", effective_irradiance(spectrum) * minutes) / standard_erythema_dose; digits = 2))
end
markdown_table(["Place and day", "Cripps and Ramsay (min)", "1 SED (min)", "1.5 SED (min)", "SED at the Cripps and Ramsay time"], comparison_rows)
```

## Speed

```@example sunburn
cells = length(cell_heights)
markdown_table(["Cells", "Seconds per map", "Microseconds per cell"], [(cells, round(seconds; sigdigits = 2), round(1e6 * seconds / cells; sigdigits = 2))])
```

The UV is calculated with the [`DaveFurukawaScattering`](@ref) model at a single time, so the map takes a fraction of a second.

## References

CIE (1998) Erythema reference action spectrum and standard erythema dose. CIE S 007/E-1998. Also ISO 17166:1999.

Cripps DJ, Ramsay CA (1970) Ultraviolet action spectrum with a prism-grating monochromator. British Journal of Dermatology 82: 584-592.

Jablonski NG, Chaplin G (2000) The evolution of human skin coloration. Journal of Human Evolution 39: 57-106.

New M, Lister D, Hulme M, Makin I (2002) A high-resolution data set of surface climate over global land areas. Climate Research 21: 1-25.
