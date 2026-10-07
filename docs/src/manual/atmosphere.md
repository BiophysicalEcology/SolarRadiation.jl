# Atmosphere

On its way to the ground the direct beam is attenuated by Rayleigh scattering by air molecules, scattering by aerosols
(dust, smoke, droplets), and absorption by ozone and water vapour. The amount of each is described by a vertical optical
depth at each wavelength, ``{}_\lambda\tau``, and is multiplied by the relative air mass ``m(Z_a)`` along the slanted path of the sun.
The direct irradiance on a horizontal surface is

```math
I_\lambda = S_\lambda \left(\frac{a}{r}\right)^2 \cos Z_a \exp\left[-{}_\lambda\tau(Z_a)\right]
```

where ``Z_a`` is the apparent zenith angle after [refraction](solar_geometry.md#Refraction-and-air-mass).

```@setup atmosphere
using Main.FigureHelpers
using CairoMakie, SolarRadiation, Unitful
using FluidProperties: atmospheric_pressure

function site(latitude; elevation = 0.0u"m", albedo = 0.2)
    SolarTerrain(;
        elevation, albedo, latitude, longitude = 0.0u"°",
        horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°",
        atmospheric_pressure = atmospheric_pressure(elevation),
    )
end
```

## Optical depths

The optical depth at wavelength ``\lambda`` and elevation ``h`` is the sum of the contributions of the constituents, for the
air mass of the sun (eq. 14 of McCullough and Porter 1971):

```math
{}_\lambda\tau(Z_a) = m(Z_a)\left[\frac{P}{P_0}\,{}_\lambda\tau_R A_R(h)
    + \frac{25}{MR_0}\,{}_\lambda\tau_A A_A(h)
    + \frac{X}{0.34}\,{}_\lambda\tau_O A_O(h)\right]
    + {}_\lambda\tau_W \sqrt{m(Z_a)\, w\, A_W(h)}
```

The sea-level optical depths ``{}_\lambda\tau`` are the tables of the [`SolarProblem`](@ref). They are for a standard atmosphere with
a pressure ``P_0`` of 101.3 kPa, a sea-level visibility ``MR_0`` (the meteorological range at 0.55 μm) of 25 km, a total
ozone column ``X`` of 0.34 cm and 1 mm of precipitable water ``w``, and are scaled to a specific location as follows:

- the Rayleigh (molecular) depth ``\tau_R`` by the ratio of the pressure ``P`` of the site to ``P_0``;
- the aerosol depth ``\tau_A`` by the ratio of 25 km to the visibility `mixing_ratio_height` of the model;
- the ozone depth ``\tau_O`` by the ratio of the ozone column of the latitude and month to 0.34 cm;
- the water vapour depth ``\tau_W`` by the square root of the product of the air mass and the precipitable water
  `precipitable_water` of the model, in mm. Water vapour is absorbing, so it is not proportional to the amount.
  Gates and Harrop (1963) define the coefficients for ``w`` in mm; McCullough and Porter (1971) took them as for 1 cm;
- the absorption by the uniformly mixed gases, O₂ near 1.27 μm and CO₂ near 2.0 μm, with the coefficients ``a_M`` of Bird
  and Riordan (1986), as ``\tau_M = 1.41\,a_M M / (1 + 118.3\,a_M M)^{0.45}``, with ``M = m(Z_a)\,P/P_0``. It depends on the
  air mass and pressure, not on the water vapour.

The factors ``A_R``, ``A_A`` and ``A_O`` adjust each depth to the elevation of the site, and ``A_W`` is 1.
The total optical depth is limited to 80. [`SolarRadiation.spectral_optical_depth`](@ref) does the calculation
at one wavelength.

The default tables of ``\tau_R``, ``\tau_A``, ``\tau_O`` and ``\tau_W`` are from Elterman (1968, 1970) for Rayleigh, ozone and
aerosol scattering, and from Table II of Gates and Harrop (1963) for water vapour:

```@example atmosphere
model = SolarProblem(; diffuse_model = ChandrasekharScattering())
wavelength = model.wavelengths

fig = Figure(size = (700, 600))
panels = (
    ("Rayleigh (molecular)", model.rayleigh_optical_depth),
    ("Ozone", model.ozone_optical_depth),
    ("Water vapour", model.water_optical_depth),
    ("Aerosol", model.aerosol_optical_depth),
)
for (i, (title, depth)) in enumerate(panels)
    ax = Axis(fig[div(i - 1, 2) + 1, mod(i - 1, 2) + 1]; title, xlabel = "Wavelength (nm)", ylabel = "Optical depth")
    lines!(ax, ustrip.(wavelength), depth; linewidth = 2)
    title == "Water vapour" && ylims!(ax, 0, 1)
end
fig
```

Ozone absorbs ultraviolet light strongly, below about 330 nm, and Rayleigh scattering is strongest at short wavelengths, while
the water vapour bands are in the infrared. The water vapour samples in the strong bands near 1.4, 1.9 and 2.7 μm are
opaque (80, the limit), and are cut off in the figure.

## Elevation

The vertical distribution of each constituent differs, and the factors ``A_R``, ``A_A`` and ``A_O`` are polynomials in
the elevation fitted to the profiles of Elterman (1968), see [`elevation_correction`](@ref).
Water vapour has no standard profile, so its factor is 1. The pressure of the site is set separately in the terrain, and can be calculated
from the elevation with `atmospheric_pressure` of [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl).

```@example atmosphere
elevations = range(-400, 6000, length = 200) .* u"m"
factors = elevation_correction.(elevations)

fig, ax = figure_axis("Elevation (m)", "Correction factor")
for name in (:molecular, :aerosol, :ozone)
    lines!(ax, ustrip.(elevations), getproperty.(factors, name); linewidth = 2, label = string(name))
end
axislegend(ax; position = :rt)
fig
```

## Ozone

The ozone column varies with latitude and season, and follows the table of Robinson (1966) in [`SolarProblem`](@ref) with
values every 10° of latitude for each month. [`SolarRadiation.ozone_depth_lookup`](@ref) takes the value for the nearest
band and the month of the day of the year.

```@example atmosphere
latitudes = -90:10:90
ozone = model.ozone_column # latitude bands by month

fig, ax = figure_axis("Latitude (°)", "Ozone column (cm)")
months = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
for (i, month) in enumerate(months)
    lines!(ax, latitudes, ozone[:, i]; linewidth = 2, label = month, color = i, colorrange = (1, 12), colormap = :twilight)
end
Legend(fig[1, 2], ax; nbanks = 1, framevisible = false)
fig
```

## From the top of the atmosphere to the ground

The spectrum at the top of the atmosphere is shaped by each constituent of the atmosphere on the way to the ground. To see the effect of each one
we calculate the direct beam with only some of them, by setting the optical depth tables of the others to a negligible value: first with only
Rayleigh scattering, then with ozone added, then aerosols, and finally water vapour, which is the complete direct beam. The
diffuse radiation, from [`ChandrasekharScattering`](@ref), is added to give the global radiation. The spectra are for solar noon on a summer day, at a sea-level
site at 37.8°S and at 5000 m, on a horizontal surface, in W m⁻² nm⁻¹, with a logarithmic wavelength axis.

```@example atmosphere
none = fill(1e-9, length(wavelength)) # not zero: a total optical depth of exactly zero gives no direct beam
direct_only(; kw...) = SolarProblem(; diffuse_model = NoScattering(), kw...)
steps = (
    "Rayleigh scattering" => direct_only(; ozone_optical_depth = none, aerosol_optical_depth = none, water_optical_depth = none),
    "+ ozone" => direct_only(; aerosol_optical_depth = none, water_optical_depth = none),
    "+ aerosols" => direct_only(; water_optical_depth = none),
    "+ water vapour (the direct beam)" => direct_only(),
)

geometry = solar_geometry(model.solar_geometry_model, -37.8u"°"; day_of_year = 15, hour_angle = 0.0u"rad")
top_of_atmosphere = model.solar_spectral_irradiance ./ 1000 .* geometry.sun_distance_factor .* cos(geometry.zenith_angle)

noon_run(problem, terrain) = solar_radiation(problem; solar_terrain = terrain, days = [15], hours = [12.0])
spectrum(x) = x[1, :]

function attenuation_figure(height)
    terrain = site(-37.8u"°"; elevation = height)
    plain(v) = ustrip.(u"W/m^2/nm", v) # numbers for the figure
    nm = ustrip.(u"nm", wavelength)
    peak = maximum(plain(top_of_atmosphere))
    fig = Figure(size = (900, 500), fontsize = 16)
    ax = Axis(fig[1, 1]; xlabel = "Wavelength (nm)", ylabel = "Irradiance (W m⁻² nm⁻¹)", xscale = log10,
        title = "Elevation $(height)", limits = ((290, 4000), (0, 1.5 * peak)),
        xticks = ([300, 400, 500, 700, 1000, 2000, 4000], ["300", "400", "500", "700", "1000", "2000", "4000"]))
    for (label, (from, to), x) in (("UV", (290, 400), 300), ("visible", (400, 700), 430), ("near infrared", (700, 4000), 760))
        vspan!(ax, from, to; color = (:gray, label == "visible" ? 0.04 : 0.12))
        text!(ax, x, 1.46 * peak; text = label, align = (:left, :top), fontsize = 20)
    end
    lines!(ax, nm, plain(top_of_atmosphere); color = :black, linewidth = 3, label = "top of the atmosphere")
    for (label, problem) in steps
        lines!(ax, nm, plain(spectrum(noon_run(problem, terrain).direct_spectra)); linewidth = 2, label)
    end
    lines!(ax, nm, plain(spectrum(noon_run(model, terrain).global_spectra)); linewidth = 2, linestyle = :dash, color = :violet, label = "global (direct + diffuse)")
    lines!(ax, nm, plain(spectrum(noon_run(model, terrain).diffuse_spectra)); linewidth = 2, color = :violet, label = "diffuse")
    axislegend(ax; position = :rt, framevisible = false, labelsize = 16)
    fig
end
attenuation_figure(0.0u"m")
```

```@example atmosphere
attenuation_figure(5000.0u"m")
```

Rayleigh scattering removes most of the radiation at the shortest wavelengths, because its optical depth falls rapidly with wavelength, roughly as
the inverse fourth power. Ozone removes the ultraviolet below about 330 nm, and the aerosols remove a little at all wavelengths, more in the
blue. The water vapour makes the gaps in the infrared. Where there is less air above the site, at 5000 m, more of each part of the spectrum reaches the ground. The
diffuse radiation is greatest in the blue and the ultraviolet, where Rayleigh scattering is strongest.

The values at each wavelength are integrated with the trapezoidal rule, over the wavelengths, to give the irradiance on a horizontal surface, in W m⁻². Each constituent takes some away
from the top of the atmosphere:

```@example atmosphere
trapezoid(y) = sum((wavelength[k+1] - wavelength[k]) * (y[k+1] + y[k]) / 2 for k in 1:length(wavelength)-1)
terrain = site(-37.8u"°")
full = noon_run(model, terrain)
watts(x) = round(uconvert(u"W/m^2", x) / u"W/m^2"; digits = 1)
rows = [("top of the atmosphere" => trapezoid(top_of_atmosphere)); [label => trapezoid(spectrum(noon_run(problem, terrain).direct_spectra)) for (label, problem) in steps];
    "diffuse" => full.diffuse_horizontal[1]; "global" => full.global_horizontal[1]]

using Markdown
Markdown.parse(join(["| Case | Irradiance (W m⁻²) |"; "| :--- | ---: |"; ["| $label | $(watts(value)) |" for (label, value) in rows]], Char(10)))
```

The ultraviolet-B, from 290 to 315 nm, is the part of the spectrum where the atmosphere has the largest effect. It is integrated in the same way in the [sunburn tutorial](../tutorials/sunburn.md).

## Sensitivity

The precipitable water and the elevation change the irradiance at noon, in W m⁻²:

```@example atmosphere
noon(model, terrain) = solar_radiation(model; solar_terrain = terrain, days = [15], hours = [12.0]).global_horizontal[1]

water = [w => noon(SolarProblem(; diffuse_model = ChandrasekharScattering(), precipitable_water = w), site(-37.8u"°")) for w in (0.1u"cm", 1.0u"cm", 2.0u"cm", 4.0u"cm")]
height = [h => noon(model, site(-37.8u"°"; elevation = h)) for h in (0.0u"m", 1000.0u"m", 2000.0u"m", 4000.0u"m")]
Markdown.parse(join(["| Varied | Value | Global irradiance at noon (W m⁻²) |"; "| :--- | ---: | ---: |";
    ["| precipitable water | $(w / u"cm") cm | $(watts(g)) |" for (w, g) in water];
    ["| elevation | $(h / u"m") m | $(watts(g)) |" for (h, g) in height]], Char(10)))
```

Higher sites receive more radiation because there is less air above them, less Rayleigh scattering and fewer aerosols.
