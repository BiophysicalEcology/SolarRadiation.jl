# Data tables

The atmosphere and the sun are described by tables at 111 wavelengths. They are the defaults of the fields of a
[`SolarProblem`](@ref) and are stored in the package as constants, which are not exported and have to be qualified with the
name of the module, as in `SolarRadiation.DEFAULT_WAVELENGTHS`. Each can be replaced by giving the keyword argument.

```@setup tables
using Main.FigureHelpers
using CairoMakie, SolarRadiation, Unitful
```

| Constant | Field | Contents | Source |
| :------- | :---- | :------- | :----- |
| `DEFAULT_WAVELENGTHS` | `wavelengths` | 111 wavelengths from 290 to 4000 nm | McCullough and Porter (1971) |
| `DEFAULT_SOLAR_SPECTRAL_IRRADIANCE` | `solar_spectral_irradiance` | extraterrestrial solar spectrum | McCullough and Porter (1971) |
| `DEFAULT_RAYLEIGH_OPTICAL_DEPTH` | `rayleigh_optical_depth` | molecular scattering, at 101.3 kPa | Elterman (1968, 1970) |
| `DEFAULT_OZONE_OPTICAL_DEPTH` | `ozone_optical_depth` | ozone absorption, for a column of 0.34 cm | Elterman (1968, 1970) |
| `DEFAULT_AEROSOL_OPTICAL_DEPTH` | `aerosol_optical_depth` | aerosols, for a visibility of 25 km | Elterman (1968, 1970) |
| `DEFAULT_WATER_OPTICAL_DEPTH` | `water_optical_depth` | water vapour, for 1 cm of precipitable water | Gates and Harrop (1963) |
| `DEFAULT_OZONE_COLUMN` | `ozone_column` | ozone column in cm by 10° latitude band and month | Robinson (1966) |
| `DEFAULT_DIFFUSE_SKY_IRRADIANCE`, `DEFAULT_DIFFUSE_GROUND_REFLECTED`, `DEFAULT_SPHERICAL_ALBEDO` | fields of [`DaveFurukawaScattering`](@ref) | ultraviolet diffuse radiation | Dave and Furukawa (1966) |

## Wavelengths and the solar spectrum

The wavelengths are closer together where the spectrum changes quickly, in the ultraviolet and visible, and further apart in
the infrared. The spectrum is the irradiance of a plane perpendicular to the rays at one astronomical unit from the sun. It is
stored as ten times the tabulated values, with a nominal unit of W m⁻² nm⁻¹, and the calculation divides by 1000, so that the irradiances
it gives are in W m⁻² nm⁻¹.

```@example tables
wavelength = ustrip.(u"nm", SolarRadiation.DEFAULT_WAVELENGTHS)
spectrum = ustrip.(u"W/m^2/nm", SolarRadiation.DEFAULT_SOLAR_SPECTRAL_IRRADIANCE) ./ 1000

fig = Figure(size = (700, 380))
ax1 = Axis(fig[1, 1]; xlabel = "Wavelength (nm)", ylabel = "Solar spectral irradiance (W m⁻² nm⁻¹)")
lines!(ax1, wavelength, spectrum; linewidth = 2)
ax2 = Axis(fig[1, 2]; xlabel = "Wavelength (nm)", ylabel = "Step to the next wavelength (nm)")
scatter!(ax2, wavelength[1:end-1], diff(wavelength); markersize = 6)
fig
```

The integral of the spectrum over the 111 wavelengths, with the trapezoidal rule, is the extraterrestrial irradiance
of the plane in the range 290 to 4000 nm:

```@example tables
sum(diff(wavelength) .* (spectrum[1:end-1] .+ spectrum[2:end]) ./ 2) * u"W/m^2"
```

## Diffuse radiation in the ultraviolet

The tables of Dave and Furukawa (1966) are for the 11 ultraviolet wavelengths, up to 360 nm, and for zenith angles from 0 to 90° in steps of
5°:

```@example tables
(sky_irradiance = size(SolarRadiation.DEFAULT_DIFFUSE_SKY_IRRADIANCE),
    ground_reflected = size(SolarRadiation.DEFAULT_DIFFUSE_GROUND_REFLECTED),
    spherical_albedo = size(SolarRadiation.DEFAULT_SPHERICAL_ALBEDO))
```

```@example tables
fig = Figure(size = (700, 380))
for (col, (title, table)) in enumerate(("Scattered from the direct beam" => SolarRadiation.DEFAULT_DIFFUSE_SKY_IRRADIANCE,
        "Reflected from the ground" => SolarRadiation.DEFAULT_DIFFUSE_GROUND_REFLECTED))
    ax = Axis(fig[1, col]; title, xlabel = "Wavelength (nm)", ylabel = col == 1 ? "Zenith angle (°)" : "", yreversed = true)
    heatmap!(ax, wavelength[1:11], 0:5:90, log10.(table); colormap = :viridis)
end
fig
```

The tables are shown on a logarithmic scale.

## Providing your own data

Any of the tables can be replaced when the model is made, as long as it has a value at each wavelength, or a different set of
`wavelengths` and `wavelength_count` is given with tables of the same length. For example the aerosol optical depth
can be the 0.01 of a very clear atmosphere at all wavelengths:

```@example tables
model = SolarProblem(; aerosol_optical_depth = fill(0.01, 111))
model.aerosol_optical_depth[1:3]
```

The [Aerosols](aerosols.md) page shows how to make profiles from the Global Aerosol Data Set. The ozone column table
is described in [Atmosphere](atmosphere.md#Ozone).
