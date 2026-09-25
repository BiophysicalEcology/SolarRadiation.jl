# Photosynthetically active radiation

Photosynthetically active radiation (PAR) is the part of the solar spectrum from 400 to 700 nm that plants use for photosynthesis. It is
measured as the photon flux density (PPFD), in micromoles of photons per square metre per second, because photosynthesis is driven by the number of photons and not
by their energy. Because the model calculates the spectrum, both are found by integrating it over the wavelengths from 400 to 700 nm.

```@setup par
using Main.FigureHelpers
```

## Calculating PAR

The energy flux of PAR, in W m⁻², is the integral of the spectral irradiance ``E_\lambda`` over 400 to 700 nm. The photon flux is the integral of the
number of photons of each wavelength, which is ``E_\lambda\,\lambda/(hc)`` for a photon of energy ``hc/\lambda``, with the Planck constant ``h``
and the speed of light ``c``, divided by the Avogadro constant ``N_A`` to give moles:

```math
\mathrm{PPFD} = \frac{1}{N_A}\int_{400\,\mathrm{nm}}^{700\,\mathrm{nm}} \frac{E_\lambda\,\lambda}{h c}\, d\lambda
```

```@example par
using SolarRadiation, FluidProperties, Unitful
using CairoMakie

trapezoid(x, y) = sum((x[k+1] - x[k]) * (y[k+1] + y[k]) / 2 for k in 1:length(x)-1)

# `spectrum` is the spectral irradiance at each of the `wavelengths`, with units. The constants are those of Unitful:
# the Planck constant `h`, the speed of light `c0` and the Avogadro constant `Na`
function par(spectrum, wavelengths)
    band = findall(w -> 400u"nm" <= w <= 700u"nm", wavelengths)
    λ, E = wavelengths[band], spectrum[band]
    energy = uconvert(u"W/m^2", trapezoid(λ, E))
    photon_flux = uconvert(u"μmol/m^2/s", trapezoid(λ, E .* λ ./ (Unitful.h * Unitful.c0)) / Unitful.Na)
    (; energy, photon_flux)
end

wavelengths = SolarProblem().wavelengths
band = wavelengths[findall(w -> 400u"nm" <= w <= 700u"nm", wavelengths)]
markdown_table(["Wavelengths in the band", "First (nm)", "Last (nm)"], [(length(band), first(band) / u"nm", last(band) / u"nm")])
```

## PAR through the day

Here we use the [`ChandrasekharScattering`](@ref) model, which includes the diffuse radiation of the visible wavelengths, at a site in
Melbourne, Australia, in summer and winter:

```@example par
model = SolarProblem(; diffuse_model = ChandrasekharScattering())
terrain = SolarTerrain(;
    elevation = 0.0u"m", albedo = 0.2, latitude = -37.8u"°", longitude = 0.0u"°",
    horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°", atmospheric_pressure = 101325.0u"Pa",
)
hours = 0:0.5:23.5
out = solar_radiation(model; solar_terrain = terrain, days = [15, 196], hours)

spectrum(x, step) = x[step, :]
total(step) = par(spectrum(out.global_spectra, step), wavelengths)
direct(step) = par(spectrum(out.direct_spectra, step), wavelengths)
diffuse(step) = par(spectrum(out.diffuse_spectra, step), wavelengths)

ppfd(f, steps) = ustrip.(u"μmol/m^2/s", [f(s).photon_flux for s in steps]) # numbers for the figure

fig = Figure(size = (700, 450))
axes = Axis[]
for (col, (title, day)) in enumerate(("15 January" => 1, "15 July" => 2))
    steps = (day - 1) * length(hours) + 1:day * length(hours)
    ax = Axis(fig[1, col]; title, xlabel = "Solar time (h)", ylabel = col == 1 ? "PPFD (µmol m⁻² s⁻¹)" : "")
    push!(axes, ax)
    lines!(ax, out.hour[steps], ppfd(total, steps); linewidth = 2, label = "total")
    lines!(ax, out.hour[steps], ppfd(direct, steps); linewidth = 2, label = "direct")
    lines!(ax, out.hour[steps], ppfd(diffuse, steps); linewidth = 2, label = "diffuse")
    col == 2 && axislegend(ax; position = :rt)
end
linkyaxes!(axes...)
ylims!(axes[1], 0, nothing)
fig
```

At solar noon the PAR is:

```@example par
noon = findfirst(==(12.0), hours)
january = (energy = total(noon).energy, photon_flux = total(noon).photon_flux,
    fraction_of_global = total(noon).energy / out.global_horizontal[noon],
    diffuse_fraction = diffuse(noon).energy / total(noon).energy)

markdown_table(["Quantity", "Value"],
    [("PAR (W m⁻²)", round(january.energy / u"W/m^2"; digits = 1)),
     ("PPFD (µmol m⁻² s⁻¹)", round(january.photon_flux / u"μmol/m^2/s"; digits = 0)),
     ("PAR as a fraction of the global radiation", round(january.fraction_of_global; digits = 3)),
     ("Diffuse fraction of the PAR", round(january.diffuse_fraction; digits = 3)),
     ("Photons per joule of PAR (µmol J⁻¹)", round(uconvert(u"μmol/J", january.photon_flux / january.energy) / u"μmol/J"; digits = 2))])
```

## Compared with a fixed fraction of the global radiation

Where the spectrum is not available, PAR is usually estimated from the global solar radiation, which is measured at many weather stations, by assuming
that a fixed fraction of it, commonly about 45 %, is PAR, and that each joule of PAR is 4.57 micromoles of photons (Thimijan and Heins 1983).
Because the spectrum of the model depends on the sun and the atmosphere, we can see how good this assumption is:

```@example par
assumed_fraction = 0.45 # of the global radiation
assumed_conversion = 4.57u"μmol/J" # photons per joule of PAR
assumed_ppfd(global_irradiance) = uconvert(u"μmol/m^2/s", assumed_fraction * assumed_conversion * global_irradiance)

global_irradiance = out.global_horizontal
lit = findall(>(20u"W/m^2"), global_irradiance) # steps with a high enough sun
model_ppfd = [total(s).photon_flux for s in eachindex(global_irradiance)]

fig = Figure(size = (700, 650))
top, bottom = Axis[], Axis[]
for (col, (title, day)) in enumerate(("15 January" => 1, "15 July" => 2))
    steps = filter(s -> s in lit, (day - 1) * length(hours) + 1:day * length(hours))
    ax1 = Axis(fig[1, col]; title, ylabel = col == 1 ? "PPFD (µmol m⁻² s⁻¹)" : "")
    push!(top, ax1)
    lines!(ax1, out.hour[steps], ustrip.(model_ppfd[steps]); linewidth = 2, label = "from the spectrum")
    lines!(ax1, out.hour[steps], ustrip.(assumed_ppfd.(global_irradiance[steps])); linewidth = 2, linestyle = :dash, label = "45 % of global × 4.57")
    col == 2 && axislegend(ax1; position = :rt)
    ax2 = Axis(fig[2, col]; xlabel = "Solar time (h)", ylabel = col == 1 ? "Assumed / spectrum" : "")
    push!(bottom, ax2)
    lines!(ax2, out.hour[steps], assumed_ppfd.(global_irradiance[steps]) ./ model_ppfd[steps]; linewidth = 2)
    hlines!(ax2, [1.0]; color = :gray, linestyle = :dash)
end
linkyaxes!(top...)
linkyaxes!(bottom...)
fig
```

The fraction of the global radiation that is PAR, and the photons for each watt, are not constant. Through these two days they are:

```@example par
fraction = [total(s).energy / global_irradiance[s] for s in lit]
photons_per_watt = [uconvert(u"μmol/J", total(s).photon_flux / total(s).energy) for s in lit]
(fraction_low, fraction_high), (photons_low, photons_high) = extrema(fraction), extrema(photons_per_watt)

using Markdown
Markdown.parse(join(["| Quantity | Minimum | Maximum |"; "| :--- | ---: | ---: |";
    "| PAR as a fraction of the global radiation | $(round(fraction_low; digits = 3)) | $(round(fraction_high; digits = 3)) |";
    "| Photons per joule of PAR (µmol J⁻¹) | $(round(photons_low / u"μmol/J"; digits = 2)) | $(round(photons_high / u"μmol/J"; digits = 2)) |"], Char(10)))
```

The 45 % is a rule of thumb from measurements, and it lies within the range of the fraction of the model, so on these days the fixed fraction is
within −7 to +6 % of the PAR of the model. At the top of the atmosphere, in the solar spectrum of the model, the fraction is lower:

```@example par
in_band = findall(w -> 400u"nm" <= w <= 700u"nm", wavelengths)
solar = SolarProblem().solar_spectral_irradiance
round(trapezoid(wavelengths[in_band], solar[in_band]) / trapezoid(wavelengths, solar); digits = 3)
```

The atmosphere raises the fraction because water vapour removes infrared light, so the fraction, and the error of a fixed one, depend on the
precipitable water (see the table below). The comparison is between two approximations, and the size of the difference depends on the fraction
that is assumed.

The spectrum also separates the direct and the diffuse PAR, which a fixed fraction of the global radiation cannot do. They are needed for the
photosynthesis of the sunlit and shaded leaves of a canopy.

## Daily PAR through the year

The photosynthesis of a plant depends on the daily total of the photons, in mol m⁻² day⁻¹. This is calculated for every tenth day at five latitudes with
hourly steps, and integrating the photon flux over the hours:

```@example par
latitudes = (60.0, 37.8, 0.0, -37.8, -60.0)
days = collect(1:10:361)
hours_hourly = 0:1:23

function daily_par(latitude)
    site = SolarTerrain(;
        elevation = 0.0u"m", albedo = 0.2, latitude = latitude * u"°", longitude = 0.0u"°",
        horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°", atmospheric_pressure = 101325.0u"Pa",
    )
    result = solar_radiation(model; solar_terrain = site, days, hours = hours_hourly)
    photons = [par(spectrum(result.global_spectra, s), wavelengths).photon_flux for s in 1:length(result.hour)]
    daily(x) = uconvert.(u"mol/m^2", vec(sum(reshape(x, length(hours_hourly), length(days)); dims = 1)) .* 1.0u"hr") # per day
    (spectrum = daily(photons), assumed = daily(assumed_ppfd.(result.global_horizontal)))
end
seconds = @elapsed yearly = [daily_par(latitude) for latitude in latitudes]

fig, ax = figure_axis("Day of the year", "Daily PAR (mol m⁻² day⁻¹)"; size = (700, 450))
for (latitude, daily) in zip(latitudes, yearly)
    lines!(ax, days, ustrip.(daily.spectrum); linewidth = 2, label = "$(latitude)°")
end
axislegend(ax; position = :cb, nbanks = 2)
fig
```

The daily totals given by the fixed fraction of the global radiation differ from those of the spectrum, as a percentage of the latter, by:

```@example par
fig, ax = figure_axis("Day of the year", "Assumed minus spectrum (%)"; size = (700, 450))
for (latitude, daily) in zip(latitudes, yearly)
    difference = ifelse.(daily.spectrum .> 0.05u"mol/m^2", 100 .* (daily.assumed .- daily.spectrum) ./ daily.spectrum, NaN)
    lines!(ax, days, difference; linewidth = 2, label = "$(latitude)°")
end
hlines!(ax, [0.0]; color = :gray, linestyle = :dash)
axislegend(ax; position = :cb, nbanks = 2)
fig
```

The summer maximum, in the high latitudes, is as large as in the tropics, because the days are long, and the winter minimum is close to zero
at 60 degrees. In the tropics the daily total varies little through the year.

## Effect of the atmosphere and elevation

PAR is changed by the aerosols, by the water vapour, which absorbs at the red end of PAR, and by elevation, where the pressure is calculated with `atmospheric_pressure` of [FluidProperties.jl](https://github.com/BiophysicalEcology/FluidProperties.jl). Below we compute the PAR at noon in January at Melbourne
for other atmospheres, in µmol m⁻² s⁻¹, from the spectrum and from the fixed fraction of the global radiation with the error of the latter:

```@example par
function noon_par(problem, terrain)
    result = solar_radiation(problem; solar_terrain = terrain, days = [15], hours = [12.0])
    model_ppfd = par(spectrum(result.global_spectra, 1), wavelengths).photon_flux
    assumed = assumed_ppfd(result.global_horizontal[1])
    (spectrum = round(u"μmol/m^2/s", model_ppfd; digits = 0), assumed = round(u"μmol/m^2/s", assumed; digits = 0),
        assumed_error_percent = round(100 * (assumed - model_ppfd) / model_ppfd; digits = 1))
end

with(; kw...) = SolarProblem(; diffuse_model = ChandrasekharScattering(), kw...)
high(height) = SolarTerrain(;
    elevation = height, albedo = 0.2, latitude = -37.8u"°", longitude = 0.0u"°",
    horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°", atmospheric_pressure = atmospheric_pressure(height),
)
cases = [
    "default" => noon_par(model, terrain),
    "hazy (5 km)" => noon_par(with(; mixing_ratio_height = 5.0u"km"), terrain),
    "very clear (100 km)" => noon_par(with(; mixing_ratio_height = 100.0u"km"), terrain),
    "dry (0.1 cm of water)" => noon_par(with(; precipitable_water = 0.1u"cm"), terrain),
    "humid (4 cm of water)" => noon_par(with(; precipitable_water = 4.0u"cm"), terrain),
    "elevation 2000 m" => noon_par(model, high(2000.0u"m")),
]

Markdown.parse(join(["| Case | From the spectrum | Fixed fraction | Error of the fixed fraction (%) |"; "| :--- | ---: | ---: | ---: |";
    ["| $name | $(round(Int, r.spectrum / u"μmol/m^2/s")) | $(round(Int, r.assumed / u"μmol/m^2/s")) | $(r.assumed_error_percent) |" for (name, r) in cases]], Char(10)))
```

## Speed

```@example par
steps = length(days) * length(hours_hourly) * length(latitudes)
markdown_table(["Steps", "Seconds", "Milliseconds per step"], [(steps, round(seconds; sigdigits = 2), round(1000 * seconds / steps; sigdigits = 2))])
```

The time is that of the [`ChandrasekharScattering`](@ref) diffuse model, which calculates the diffuse radiation of the visible wavelengths, so that PAR is complete.

## References

Thimijan RW, Heins RD (1983) Photometric, radiometric, and quantum light units of measure: a review of procedures for interconversion. HortScience 18: 818-822.
