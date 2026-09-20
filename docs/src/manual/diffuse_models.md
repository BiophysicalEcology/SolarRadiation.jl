# Diffuse models

Some of the radiation removed from the direct beam by scattering reaches the ground as diffuse skylight, ``D_\lambda``, so that
the global radiation is ``G_\lambda = I_\lambda + D_\lambda``. The diffuse radiation depends on the wavelength, the solar zenith
angle, the optical properties of the atmosphere and the reflectance of the ground below it, which scatters some of the
radiation again. 

Diffuse irradiance comes from all directions and thus is absorbed proportional to the total surface area of an
object facing the sky whereas direct irradiance is absorbed proportional to the "sillhouette area" (these areas are calculated for
different shaped objects by the [BiophysicalGeometry.jl](https://github.com/BiophysicalEcology/BiophysicalGeometry.jl) package).

There are three ways to calculate diffuse irradiance in SolarRadiation.jl, chosen with the `diffuse_model` of the 
[`SolarProblem`](@ref):

| Model | Wavelengths | Basis | Cost |
| :---- | :---------- | :---- | :--- |
| [`NoScattering`](@ref) | none | | none |
| [`DaveFurukawaScattering`](@ref) | ultraviolet, 290 to 360 nm | tables of Dave and Furukawa (1966) | low |
| [`ChandrasekharScattering`](@ref) | all | X and Y functions of Chandrasekhar (1960) | high |

They are subtypes of [`AbstractDiffuseModel`](@ref). The default is [`DaveFurukawaScattering`](@ref), which is fast but
gives diffuse radiation only in the ultraviolet, a small part of the skylight, so that with the default the global irradiance is
mostly the direct beam. [`ChandrasekharScattering`](@ref) includes all wavelengths, and is used where the diffuse
radiation matters, at more than a thousand times the cost per time step.

```@setup diffuse
using Main.FigureHelpers
using CairoMakie, SolarRadiation, Unitful

site = SolarTerrain(;
    elevation = 0.0u"m", albedo = 0.2, latitude = -37.8u"°", longitude = 0.0u"°",
    horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°",
    atmospheric_pressure = 101325.0u"Pa",
)
```

## Dave and Furukawa

For the ultraviolet region, where ozone absorbs as well as scatters and the Rayleigh calculation below is not adequate,
McCullough and Porter (1971) use the tables that Dave and Furukawa (1966) computed for a Rayleigh atmosphere with ozone. The
diffuse radiation is

```math
D_\lambda = \frac{S_\lambda}{\pi}\left(\frac{a}{r}\right)^2
    \left[F_d(Z, \lambda) + \frac{F_d'}{Q}(Z, \lambda)\, Q(\lambda)\right],
    \qquad Q(\lambda) = \frac{A(\lambda)}{1 - A(\lambda)\,\bar s(\lambda)}
```

where ``F_d`` is the downward flux of scattered radiation, of all orders of scattering, at the ground, due to the illumination of the atmosphere by the
direct solar beam, and ``F_d'/Q`` is the quantity that gives the contribution to it from the illumination of the atmosphere from below by the radiation
reflected by the ground (Dave and Furukawa 1966, eqs. 38 and 44). The fluxes are on a horizontal surface for an incident beam of ``\pi`` units of flux normal
to the beam, which is why they are multiplied by ``S_\lambda/\pi``. ``Q`` is the factor for the repeated reflections between the ground, with albedo ``A``, and the
atmosphere (their eq. 31). In it ``\bar s`` is their ``S^b``, from Table B: the fraction of the radiation reflected by the ground that the atmosphere scatters back to it.

Dave and Furukawa tabulate these at 16 wavelengths and five solar zenith angles, 0°, 30°, 60°, 75° and 85°, so the tables of the model, for 11 wavelengths and zenith
angles in steps of 5°, are interpolated from them, and the nearest zenith angle is used. There are no tables for longer wavelengths, so the diffuse radiation is zero there.

The tables are for the sea-level surface, at 1000 mb, and a total ozone column of 341 m atm-cm (0.34 cm). McCullough and Porter (1971) note that for another elevation the
fluxes should be found from the optical depth ``\tau_R + \tau_O`` at that elevation. This is not done (yet), so the diffuse ultraviolet radiation does not change with the elevation
or the ozone column of the model, although the direct beam does.

### More on Dave and Furukawa's tables

The light reaching the ground has two parts: the direct beam, and skylight, which is sunlight that the air has scattered on the way. In the ultraviolet the
air also absorbs, through ozone, so the skylight is taken from tables computed for a realistic atmosphere. The three tables are:

- **`sky_irradiance` (``F_d``): skylight from the sun's own light.** How much of the sun's beam ends up on the ground as skylight, after being scattered any number
  of times and partly absorbed by ozone. It is largest when the sun is high, falls towards zero as the sun nears the horizon, and rises steeply
  with wavelength, from minimal at 290 nm to maximal by 330 nm, because ozone blocks the shortest ultraviolet.
- **`ground_reflected` (``F_d'/Q``): skylight from light bounced off the ground.** Some sunlight reflects off the ground, and the air above scatters part of it back down
  as extra skylight. This table says how much extra skylight there is for a given amount of reflection, and ``Q`` scales it by how reflective the ground is, so bright
  ground such as snow adds much more diffuse ultraviolet than dark ground.
- **`spherical_albedo` (``\bar s``): the bounce-back fraction.** Of the light reflected up by the ground, this is the 
  fraction that the air scatters back down. Light can bounce between the ground and the sky several times, and the bounce-back 
  fraction sets how much the repeated bouncing adds to the ground-reflection term.

Thus, the diffuse ultraviolet radiation is the strength of the sun times the skylight from the direct sunlight, plus the skylight from ground reflection, which
depends on how reflective the ground and the augmentation caused by the bouncing between ground and sky.

The tables are fields of the model, so different ones can be used:

```@example diffuse
DaveFurukawaScattering()
```

The three fields are `sky_irradiance` (``F_d``), `ground_reflected` (``F_d'/Q``) and `spherical_albedo` (``\bar s = S^b``, the spherical albedo of the atmosphere).
For example, without the repeated reflection between the ground and the atmosphere, which is the factor ``1/(1 - A\bar s)``:

```@example diffuse
no_reflection = DaveFurukawaScattering(; spherical_albedo = zeros(11))
nothing # hide
```

## Chandrasekhar

For the other wavelengths, McCullough and Porter (1971) use the solution of Chandrasekhar (1960) for the transfer of radiation in a
plane-parallel Rayleigh atmosphere. Writing ``I_{0,\lambda}`` for the extraterrestrial irradiance on the horizontal (eq. 15):

```math
D_\lambda = I_{0,\lambda}\left[\frac{\gamma_l + \gamma_r}{2\left(1 - A(\lambda)\,\bar s\right)} - e^{-{}_\lambda\tau_R\, m(Z_a)}\right]
```

where ``\gamma_l`` and ``\gamma_r`` are functions of the zenith angle and the Rayleigh optical depth ``{}_\lambda\tau_R`` for the
two directions of polarization of the scattered light, and ``\bar s`` depends on ``{}_\lambda\tau_R`` only. They are obtained
from Chandrasekhar's ``X`` and ``Y`` functions, which are found by iteration from the fourth approximation of Chandrasekhar (1960), Sec. 59, using the
method in Sec. 60, and stop when successive ``Y`` differ by less than 2 parts in 10 000 or after 15 iterations.
This is the most costly part of the calculation, and it is done for every wavelength and time.

[`scattered_radiation`](@ref) returns the functions for an optical thickness (up to 2) at 101 values of the cosine of the
zenith angle:

```@example diffuse
scattering = scattered_radiation(0.5)
keys(scattering)
```

```@example diffuse
μ = range(0, 1, length = 101)

fig, ax = figure_axis("Cosine of the zenith angle", "γ")
lines!(ax, μ, scattering.polarization_l; linewidth = 2, label = "γₗ")
lines!(ax, μ, scattering.polarization_r; linewidth = 2, label = "γᵣ")
axislegend(ax; position = :lt)
fig
```

The functions are used where the Rayleigh optical depth of the wavelength is at least 0.03 and are zero elsewhere. With the default
tables this is

```@example diffuse
last_wavelength = SolarProblem().wavelengths[findlast(>=(0.03), SolarProblem().rayleigh_optical_depth)]
```

## Comparison

The diffuse spectrum at noon in summer at the site, for the three models, and for all wavelengths and the ultraviolet:

```@example diffuse
models = (NoScattering(), DaveFurukawaScattering(), ChandrasekharScattering())
results = map(models) do diffuse_model
    solar_radiation(SolarProblem(; diffuse_model); solar_terrain = site, days = [15], hours = [12.0])
end
wavelength = ustrip.(u"nm", SolarProblem().wavelengths)

fig = Figure(size = (700, 380))
for (col, xlimit) in enumerate(((290, 1000), (290, 400)))
    ax = Axis(fig[1, col]; xlabel = "Wavelength (nm)", ylabel = col == 1 ? "Diffuse irradiance (W m⁻² nm⁻¹)" : "",
        limits = (xlimit, nothing))
    for (model, result) in zip(models, results)
        lines!(ax, wavelength, ustrip.(u"W/m^2/nm", result.diffuse_spectra[1, :]);
            linewidth = 2, label = string(nameof(typeof(model))))
    end
    col == 1 && axislegend(ax; position = :rt)
end
fig
```

From 330 to 360 nm the two calculations agree closely. At shorter wavelengths the values of Dave and Furukawa are lower, because their
tables include the absorption of the diffuse radiation by ozone, which the Rayleigh atmosphere of Chandrasekhar does not, and
only [`ChandrasekharScattering`](@ref) gives diffuse radiation beyond 360 nm, which is nearly all of it. The diffuse and
global irradiance integrated over all wavelengths, in W m⁻², are:

```@example diffuse
markdown_table(["Model", "Diffuse (W m⁻²)", "Global (W m⁻²)"],
    [(nameof(typeof(model)), round(result.diffuse_horizontal[1] / u"W/m^2"; digits = 2), round(result.global_horizontal[1] / u"W/m^2"; digits = 1))
        for (model, result) in zip(models, results)])
```

The models differ greatly in cost. The time to calculate four days of hourly radiation, after compilation, is:

```@example diffuse
days = [15, 105, 196, 288]
hours = 0:23
timings = map(models) do diffuse_model
    model = SolarProblem(; diffuse_model)
    solar_radiation(model; solar_terrain = site, days, hours) # compile
    @elapsed solar_radiation(model; solar_terrain = site, days, hours)
end
markdown_table(["Model", "Time (ms)"], [(nameof(typeof(model)), round(1000 * seconds; sigdigits = 2)) for (model, seconds) in zip(models, timings)])
```

## Adding a diffuse model

The diffuse model is chosen by dispatch, so another model is a new subtype of [`AbstractDiffuseModel`](@ref) with a method of
`SolarRadiation.diffuse_irradiance` for it. The method is called for each wavelength with the wavelength index, the Rayleigh
optical depth (dimensionless), and the parameters of the time step, and returns the irradiance in W m⁻² nm⁻¹. Its parameters can be
fields of the model. The example new method below takes the diffuse irradiance to be a fraction of the extraterrestrial radiation 
on the horizontal:

```@example diffuse
struct FractionOfExtraterrestrial{F} <: AbstractDiffuseModel
    fraction::F
end

function SolarRadiation.diffuse_irradiance(model::FractionOfExtraterrestrial, wavelength_index, rayleigh_optical_depth, params, buffers)
    (; solar_spectral_irradiance, sun_distance_factor, cosine_zenith) = params
    return model.fraction * solar_spectral_irradiance[wavelength_index] * sun_distance_factor * cosine_zenith / 1000
end

model = SolarProblem(; diffuse_model = FractionOfExtraterrestrial(0.01))
only(solar_radiation(model; solar_terrain = site, days = [15], hours = [12.0]).diffuse_horizontal)
```

A model that needs working arrays can extend `allocate_buffers`, as [`ChandrasekharScattering`](@ref) does.
