# Poles, tropics and equator

The functions of this package cover all latitudes, from the polar regions where the sun does not rise for months, to the equator where it
is high all year. This page compares seven sites for a clear-sky year, every hour of solar time on every fifth day and on the
solstices, with the [`ChandrasekharScattering`](@ref) model so that the diffuse radiation is complete:

| Site | Latitude | Elevation | Albedo |
| :--- | -------: | --------: | -----: |
| North Pole | 90°N | sea level | 0.8, sea ice |
| Arctic | 70°N | sea level | 0.2 |
| Mid-latitude | 45°N | sea level | 0.2 |
| Tropic of Cancer | 23.44°N | sea level | 0.2 |
| Equator | 0° | sea level | 0.2 |
| Tropic of Capricorn | 23.44°S | sea level | 0.2 |
| South Pole | 90°S | 2835 m | 0.8, snow |

The pressure of each site is calculated from its elevation with `atmospheric_pressure` of [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl).

```@setup latitudes
using Main.FigureHelpers
using CairoMakie, SolarRadiation, Unitful, Statistics
using FluidProperties: atmospheric_pressure
```

```@example latitudes
sites = [
    (name = "North Pole", latitude = 90.0, elevation = 0.0, albedo = 0.8),
    (name = "Arctic", latitude = 70.0, elevation = 0.0, albedo = 0.2),
    (name = "Mid-latitude", latitude = 45.0, elevation = 0.0, albedo = 0.2),
    (name = "Tropic of Cancer", latitude = 23.44, elevation = 0.0, albedo = 0.2),
    (name = "Equator", latitude = 0.0, elevation = 0.0, albedo = 0.2),
    (name = "Tropic of Capricorn", latitude = -23.44, elevation = 0.0, albedo = 0.2),
    (name = "South Pole", latitude = -90.0, elevation = 2835.0, albedo = 0.8),
]

function terrain((; latitude, elevation, albedo))
    SolarTerrain(;
        elevation = elevation * u"m", albedo, latitude = latitude * u"°", longitude = 0.0u"°",
        horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°",
        atmospheric_pressure = atmospheric_pressure(elevation * u"m"),
    )
end

days = sort(union(1:5:365, [172, 355])) # every fifth day, and the solstices
hours = 0:1:23
model = SolarProblem(; diffuse_model = ChandrasekharScattering())
years = [solar_radiation(model; solar_terrain = terrain(site), days, hours) for site in sites]
nothing # hide
```

The output has one value per time step, in order of day and then hour, so it can be reshaped into hours by days:

```@example latitudes
by_day(x) = reshape(x, length(hours), length(days))
hourly = [by_day(year.global_horizontal) for year in years]
size(hourly[1])
```

## Daily radiation through the year

The radiation on a horizontal surface, summed over each day:

```@example latitudes
daily = [uconvert.(u"MJ/m^2", vec(sum(h; dims = 1)) .* 1u"hr") for h in hourly] # per day

fig, ax = figure_axis("Day of the year", "Global irradiance (MJ m⁻² day⁻¹)"; size = (700, 450))
for (site, total) in zip(sites, daily)
    lines!(ax, days, ustrip.(total); linewidth = 2, label = site.name)
end
Legend(fig[1, 2], ax)
fig
```

- At the **poles** the radiation is zero through the polar night. It starts when the sun rises, at the equinox, and
  is greatest at the solstice, when the sun is up for 24 hours, and falls to zero again in the following autumn. The South
  Pole is higher than the North Pole because the site is at 2835 m, with less air above it, and the earth is closest to the sun in January.
- At the **Arctic** site the polar night is shorter, and midnight sun and long days are in summer.
- The **tropics** have a single peak, in their summer, and the lowest radiation in their winter, when the sun is
  further from the zenith.
- The **equator** has two maxima, around the equinoxes, when the sun is overhead at noon, and lower values
  at the solstices.

The mean, maximum and minimum of the daily radiation of the days calculated, in MJ m⁻² day⁻¹:

```@example latitudes
using Statistics, Markdown
Markdown.parse(join(["| Site | Mean | Maximum | Minimum |"; "| :--- | ---: | ---: | ---: |";
    ["| $(site.name) | $(round(mean(total) / u"MJ/m^2"; digits = 1)) | $(round(maximum(total) / u"MJ/m^2"; digits = 1)) | $(round(minimum(total) / u"MJ/m^2"; digits = 1)) |"
        for (site, total) in zip(sites, daily)]], Char(10)))
```

## Day length

The length of the day is `2 * hour_angle_sunrise`, twice the time between sunrise and solar noon. See
[Sunrise, sunset and day length](solar_geometry.md#Sunrise,-sunset-and-day-length).

```@example latitudes
fig, ax = figure_axis("Day of the year", "Day length (h)"; size = (700, 450), yticks = 0:4:24)
for (site, year) in zip(sites, years)
    lines!(ax, days, 2 .* year.hour_angle_sunrise; linewidth = 2, label = site.name)
end
Legend(fig[1, 2], ax)
fig
```

Poleward of the polar circles at 66.5° the day length is 0 or 24 hours around the solstices. The change in day length
in the tropics is small, and at the equator it is 12 hours all year.

## The course of the day and the year

The irradiance by time of day and day of the year is a picture of the path of the sun:

```@example latitudes
panels = ("North Pole", "Tropic of Cancer", "Equator")
fig = Figure(size = (700, 900))
heatmaps = map(enumerate(panels)) do (row, name)
    i = findfirst(site -> site.name == name, sites)
    ax = Axis(fig[row, 1]; title = name, ylabel = "Solar time (h)", xlabel = row == 3 ? "Day of the year" : "",
        yticks = 0:6:24)
    heatmap!(ax, days, collect(hours), ustrip.(hourly[i])'; colorrange = (0, 1100))
end
Colorbar(fig[1:3, 2], first(heatmaps); label = "Global irradiance (W m⁻²)")
fig
```

At the North Pole the sun is up all day for half of the year, so the irradiance does not vary with the time of day and the pole
has no night in summer (the midnight sun). The image is a band that widens and narrows with the seasons. At the Tropic of Cancer the sun
is overhead at noon on the June solstice, and the equator has two seasons of the highest sun.

## Spectra

Here are the direct and diffuse spectra at solar noon on the December solstice, when the sun is high in the tropics and low at the South Pole, which is at 2835 m and has a bright snow surface:

```@example latitudes
noon = (findfirst(==(355), days) - 1) * length(hours) + findfirst(==(12.0), hours)
selected = ("Equator", "Tropic of Capricorn", "South Pole")
wavelength = ustrip.(u"nm", model.wavelengths)
spectrum(year, x) = ustrip.(u"W/m^2/nm", getproperty(year, x)[noon, :])

fig, ax = figure_axis("Wavelength (nm)", "Irradiance (W m⁻² nm⁻¹)"; size = (700, 450))
for (color, name) in zip((:royalblue, :crimson, :seagreen), selected)
    year = years[findfirst(site -> site.name == name, sites)]
    lines!(ax, wavelength, spectrum(year, :direct_spectra); color, linewidth = 2, label = "$name, direct")
    lines!(ax, wavelength, spectrum(year, :diffuse_spectra); color, linewidth = 2, linestyle = :dash, label = "$name, diffuse")
end
xlims!(ax, 290, 2500)
axislegend(ax; position = :rt)
fig
```

## Direct and diffuse radiation

Here is the proportion of the radiation that is diffuse at solar noon, by latitude, on the two solstices and the equinox:

```@example latitudes
latitudes = -90:2:90
function diffuse_fraction(latitude, day)
    site = (; latitude = float(latitude), elevation = 0.0, albedo = 0.2)
    out = solar_radiation(model; solar_terrain = terrain(site), days = [day], hours = [12.0])
    global_horizontal = out.global_horizontal[1]
    global_horizontal > 1e-3u"W/m^2" ? out.diffuse_horizontal[1] / global_horizontal : NaN
end

fig, ax = figure_axis("Latitude (°)", "Diffuse fraction of global irradiance"; size = (700, 450))
for (day, label) in ((172, "June solstice"), (80, "March equinox"), (355, "December solstice"))
    lines!(ax, latitudes, [diffuse_fraction(latitude, day) for latitude in latitudes]; linewidth = 2, label)
end
axislegend(ax; position = :ct)
fig
```

The proportion is smallest where the sun is high, because the path of the light through the air is short, and it
increases toward the horizon. The lines end in the regions of the polar night.

The radiation at the horizon with the sun below it is skylight (twilight), and is included for zenith angles up to 107°.
