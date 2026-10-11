# Aerosols

Aerosols, the dust, smoke, salt and droplets suspended in the air, scatter and absorb sunlight. Their optical depth at each wavelength,
`aerosol_optical_depth` of the [`SolarProblem`](@ref), attenuates the direct beam and, with the visibility, is the most variable
part of the [atmosphere](atmosphere.md) of the model. It is scaled by the sea-level visibility `mixing_ratio_height` (25 km by default) and adjusted for
[elevation](atmosphere.md#Elevation).

```@setup aerosols
using Main.FigureHelpers
using CairoMakie, SolarRadiation, Unitful
```

## The default profile

The default profile is that of Elterman (1968, 1970). It is based on observations in North America and can differ strongly from
regions elsewhere, such as Australia. Its optical depth falls with wavelength, from 0.27 in the ultraviolet to 0.01 in the infrared:

```@example aerosols
model = SolarProblem(; diffuse_model = ChandrasekharScattering())
wavelength = ustrip.(u"nm", model.wavelengths)
markdown_table(["Wavelength (nm)", "Aerosol optical depth"],
    [(290, round(model.aerosol_optical_depth[1]; digits = 3)),
     (500, round(model.aerosol_optical_depth[findfirst(==(500.0), wavelength)]; digits = 3)),
     (4000, round(model.aerosol_optical_depth[end]; digits = 3))])
```

Visibility changes the amount of aerosol, and the radiation at noon on a summer day at 37.8°S:

```@example aerosols
terrain(latitude) = SolarTerrain(;
    elevation = 0.0u"m", albedo = 0.2, latitude = latitude * u"°", longitude = 0.0u"°",
    horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°", atmospheric_pressure = 101325.0u"Pa",
)

function noon(problem, latitude, day)
    out = solar_radiation(problem; solar_terrain = terrain(latitude), days = [day], hours = [12.0])
    (global_horizontal = out.global_horizontal[1], direct = out.direct_horizontal[1], diffuse = out.diffuse_horizontal[1])
end

watts(x) = round(x / u"W/m^2"; digits = 1)
hazy = [visibility => noon(SolarProblem(; diffuse_model = ChandrasekharScattering(), mixing_ratio_height = visibility), -37.8, 15) for visibility in (5.0, 10.0, 25.0, 50.0, 100.0) .* u"km"]

markdown_table(["Visibility (km)", "Global (W m⁻²)", "Direct (W m⁻²)", "Diffuse (W m⁻²)"],
    [(visibility / u"km", watts(r.global_horizontal), watts(r.direct), watts(r.diffuse)) for (visibility, r) in hazy])
```

A hazy atmosphere reduces the direct beam and so the global irradiance. The diffuse radiation of the model does not change, because it is calculated for a Rayleigh atmosphere (see [Diffuse models](diffuse_models.md)), so the scattering by aerosols does not add to it.

## The Global Aerosol Data Set

The [Global Aerosol Data Set](http://www.lrz.de/~uh234an/www/radaer/gads.html) (GADS, Koepke et al. 1997) gives the
optical properties of aerosols on a grid of 5° by 5° for the globe, for the summer and winter, for 8 values of the relative humidity.
The aerosols are mixtures of components such as water soluble aerosol, soot, mineral dust and sea salt, in proportions that depend on the region. The vertical optical depth at 25 wavelengths from 250 to
4000 nm is available through [RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl), which downloads a
NetCDF file of 4 MB the first time:

```@example aerosols
using RasterDataSources, NCDatasets

gads_path = getraster(GADS)
gads = NCDataset(gads_path)
gads_lon, gads_lat = gads["lon"][:], gads["lat"][:]
gads_relhum, gads_wavelength = gads["relhum"][:], gads["wavelength"][:]
gads_depth = Float64.(coalesce.(gads["OPTDEPTH"][:, :, :, :, :], NaN)) # lon, lat, relative humidity, season, wavelength
close(gads)
size(gads_depth)
```

The file needs a folder for the download, which is set in the environment variable `RASTERDATASOURCES_PATH`.
Its dimensions are longitude, latitude, relative humidity (0, 50, 70, 80, 90, 95, 98 and 99 %), season (Northern Hemisphere
summer in July and winter in January) and the 25 wavelengths.

To use it in the model, the profile for a place is interpolated to the 111 wavelengths of the model. Here the nearest grid cell and relative humidity are
used, the summer and winter values are blended by month, and wavelengths are interpolated linearly:

```@example aerosols
nearest(values, x) = argmin(abs.(values .- x))

function interpolate_linear(x, y, xnew)
    map(xnew) do q
        q <= first(x) && return first(y)
        q >= last(x) && return last(y)
        k = searchsortedlast(x, q)
        t = (q - x[k]) / (x[k+1] - x[k])
        (1 - t) * y[k] + t * y[k+1]
    end
end

function gads_aerosol_optical_depth(lon, lat, relative_humidity, month; wavelengths = wavelength)
    i, j, k = nearest(gads_lon, lon), nearest(gads_lat, lat), nearest(gads_relhum, relative_humidity)
    summer, winter = gads_depth[i, j, k, 1, :], gads_depth[i, j, k, 2, :]
    weight = 0.5 * (1 - cos(2π * (month - 1) / 12)) # 1 in July and 0 in January
    depth = weight .* summer .+ (1 - weight) .* winter
    return max.(interpolate_linear(gads_wavelength, depth, wavelengths), 0.0)
end
nothing # hide
```

This is one of several ways of doing it: the microclimate model of NicheMapR uses the nearest cell, 0 % relative humidity, the mean of the
two seasons and a polynomial fit over wavelength, and MicroclimateMapper.jl has a version that
interpolates between cells and relative humidities (`get_aerosol_optical_depth`).

## Five locations

The optical depth at five locations with different aerosols, in the local summer (January in the south and July in the north,
with a day in the middle of the month) and at a relative humidity of 50 %, compared to the default profile:

```@example aerosols
places = [
    (name = "South Pole", lon = 0.0, lat = -90.0, month = 1, day = 15),
    (name = "Amazon (Manaus)", lon = -60.0, lat = -3.1, month = 1, day = 15),
    (name = "Sahara (Tamanrasset)", lon = 5.5, lat = 22.8, month = 7, day = 196),
    (name = "Beijing", lon = 116.4, lat = 39.9, month = 7, day = 196),
    (name = "Melbourne", lon = 145.0, lat = -37.8, month = 1, day = 15),
]
depths = [gads_aerosol_optical_depth(place.lon, place.lat, 50.0, place.month) for place in places]

fig, ax = figure_axis("Wavelength (nm)", "Aerosol optical depth"; size = (700, 450), xscale = log10, yscale = log10)
lines!(ax, wavelength, model.aerosol_optical_depth; color = :black, linewidth = 3, linestyle = :dash, label = "default")
for (place, depth) in zip(places, depths)
    lines!(ax, wavelength, depth; linewidth = 2, label = place.name)
end
axislegend(ax; position = :lb)
fig
```

The effect on the clear-sky radiation at solar noon at these latitudes, with the aerosols of GADS in place of the default, and
0.2 albedo at sea level, in W m⁻²:

```@example aerosols
rows = map(zip(places, depths)) do (place, depth)
    default = noon(model, place.lat, place.day)
    gads = noon(SolarProblem(; diffuse_model = ChandrasekharScattering(), aerosol_optical_depth = depth), place.lat, place.day)
    (place = place.name, default_global = default.global_horizontal, gads_global = gads.global_horizontal,
        gads_direct = gads.direct, gads_diffuse = gads.diffuse)
end

markdown_table(["Place", "Global, default aerosols", "Global, GADS", "Direct, GADS", "Diffuse, GADS"],
    [(r.place, watts(r.default_global), watts(r.gads_global), watts(r.gads_direct), watts(r.gads_diffuse)) for r in rows])
```

The profile has to be for the same wavelengths as the model. The [Mapping tutorial](../tutorials/mapping.md) does this
for all the cells of the globe.

## Relative humidity

Aerosols take up water and swell as the air becomes humid, which increases their optical depth. For Beijing:

```@example aerosols
beijing = places[4]
fig, ax = figure_axis("Wavelength (nm)", "Aerosol optical depth"; size = (700, 450), xscale = log10, yscale = log10)
for relhum in (0.0, 50.0, 80.0, 95.0, 99.0)
    lines!(ax, wavelength, gads_aerosol_optical_depth(beijing.lon, beijing.lat, relhum, beijing.month); linewidth = 2, label = "$(Int(relhum)) %")
end
axislegend(ax; position = :lb)
fig
```
