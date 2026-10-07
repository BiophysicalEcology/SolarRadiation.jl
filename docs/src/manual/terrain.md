# Terrain

The terrain of a site is described by a [`SolarTerrain`](@ref). Elevation and pressure change the [atmosphere](atmosphere.md)
above the site and its albedo affects the [diffuse radiation](diffuse_models.md). This page is about the parts of the terrain
that change the geometry: the slope and the aspect of the surface, and the hills around it.

```@setup terrain
using Main.FigureHelpers
using CairoMakie, SolarRadiation, Unitful

function site(; slope = 0.0u"°", aspect = 0.0u"°", horizon_angles = fill(0.0u"°", 24), latitude = -37.8u"°")
    SolarTerrain(;
        elevation = 0.0u"m", albedo = 0.2, latitude, longitude = 0.0u"°",
        horizon_angles, slope, aspect, atmospheric_pressure = 101325.0u"Pa",
    )
end
model = SolarProblem(; diffuse_model = ChandrasekharScattering())
hours = 0:23
simulate(terrain, days) = solar_radiation(model; solar_terrain = terrain, days, hours)
```

## Slope and aspect

The `slope` of the surface is its angle from the horizontal, and its `aspect` is the azimuth it faces, clockwise from north, so
that north is 0°, east is 90°, south is 180° and west is 270°. A slope of 0° is flat, and aspect is ignored.

The angle between the sun and the perpendicular of a sloping surface is (Sellers 1965, eq. 3.15)

```math
Z_{SL} = \arccos\left[\cos Z \cos SL + \sin Z \sin SL \cos(Az_{sun} - Az_{SL})\right]
```

where ``SL`` is the slope, ``Az_{SL}`` is its aspect and ``Az_{sun}`` is the [azimuth of the sun](solar_geometry.md#Azimuth). It is
limited to 90°, when the sun is behind the slope. [`SolarRadiation.slope_zenith_angle`](@ref) calculates it, and the output of
[`solar_radiation`](@ref) has it as `zenith_slope_angle`.

The global irradiance on the slope is the horizontal one, scaled by the ratio of the cosines of the angles (see
[`SolarRadiation.terrain_irradiance`](@ref)) while the sun is above the horizon:

```math
G_{SL} = G\,\frac{\cos Z_{SL}}{\cos Z}
```

and is `global_terrain` in the output. This treats all the radiation as if it came from the sun, which
overestimates the radiation on slopes that face away from it, and underestimates the diffuse radiation that a slope facing
the sun receives from a large part of the sky.

The daily course of the irradiance at 37.8°S on a winter day, when the sun is in the north, for slopes of 30° facing the
four compass directions, is calculated like this:

```@example terrain
winter = 196
fig, ax = figure_axis("Solar time (h)", "Global irradiance on the slope (W m⁻²)"; size = (700, 450))
flat = simulate(site(), [winter])
lines!(ax, flat.hour, ustrip.(flat.global_terrain); color = :black, linewidth = 2, label = "flat")
for (aspect, label) in ((0, "north"), (90, "east"), (180, "south"), (270, "west"))
    out = simulate(site(; slope = 30.0u"°", aspect = aspect * 1.0u"°"), [winter])
    lines!(ax, out.hour, ustrip.(out.global_terrain); linewidth = 2, label)
end
axislegend(ax; position = :lt)
fig
```

The slope facing the sun gets more than the flat surface and the one facing away gets little, and east and west slopes
are lit in the morning and in the afternoon. Over the day, the effect of the slope and the aspect is:

```@example terrain
slopes = 0:10:60
aspects = 0:30:330
daily(terrain, day) = uconvert(u"MJ/m^2", sum(simulate(terrain, [day]).global_terrain) * 1u"hr") # per day
totals = [daily(site(; slope = slope * 1.0u"°", aspect = aspect * 1.0u"°"), winter) for aspect in aspects, slope in slopes]

fig = Figure(size = (700, 450))
ax = Axis(fig[1, 1]; xlabel = "Aspect (°)", ylabel = "Slope (°)", xticks = 0:90:360, title = "Daily global irradiance, 15 July (MJ m⁻² day⁻¹)")
hm = heatmap!(ax, aspects, slopes, ustrip.(totals))
Colorbar(fig[1, 2], hm)
fig
```

## Hills

The horizon angles of the site describe the hills around it. They are the elevation angles of the horizon
in a number of directions at equal steps of azimuth, with the first towards north, so that 24 angles are 15° apart. When
the altitude of the sun is below the horizon angle in its direction, the direct radiation is zero, which
is the shade of the hills (Dozier et al. 1981). The diffuse radiation is not reduced. The angles can be calculated from a digital
elevation model with a GIS. [`SolarRadiation.horizon_angle_at_azimuth`](@ref) finds the angle nearest to the azimuth of the sun.

A ridge of 25° to the east, in the directions from north-east to south-east, delays the sun in the morning:

```@example terrain
horizon_angles = fill(0.0u"°", 24)
horizon_angles[4:8] .= 25.0u"°" # azimuths 45° to 105°
summer = 15

fig, ax = figure_axis("Solar time (h)", "Direct irradiance (W m⁻²)"; size = (700, 400))
for (terrain, label) in ((site(), "open"), (site(; horizon_angles), "ridge to the east"))
    out = simulate(terrain, [summer])
    lines!(ax, out.hour, ustrip.(out.direct_horizontal); linewidth = 2, label)
end
xlims!(ax, 4, 14)
axislegend(ax; position = :lt)
fig
```

The direct radiation is at its minimum value of 10⁻²⁴ W m⁻² nm⁻¹ per wavelength while the sun is hidden by the ridge.

## Sites at other latitudes

The effect of a slope depends on the latitude. The daily irradiance of 30° slopes facing north and south on the December
solstice, relative to that of a flat surface, at five latitudes:

```@example terrain
gain = map((-60.0, -37.8, 0.0, 37.8, 60.0)) do latitude
    tilted(aspect) = daily(site(; latitude = latitude * u"°", slope = 30.0u"°", aspect = aspect * 1.0u"°"), 355)
    level = daily(site(; latitude = latitude * u"°"), 355)
    (latitude = latitude, facing_north = round(tilted(0) / level; digits = 2), facing_south = round(tilted(180) / level; digits = 2))
end

using Markdown
Markdown.parse(join(["| Latitude (°) | Facing north | Facing south |"; "| ---: | ---: | ---: |";
    ["| $(g.latitude) | $(g.facing_north) | $(g.facing_south) |" for g in gain]], Char(10)))
```
