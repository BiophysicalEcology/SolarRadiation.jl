# Solar geometry

The position of the sun, and the irradiance at the top of the atmosphere, follow from the orbit of the earth, the day of
the year, the time of day and the latitude. The equations are those of McCullough and Porter (1971).

```@setup geometry
using Main.FigureHelpers
using CairoMakie, SolarRadiation, Unitful
```

## Extraterrestrial irradiance

The spectral irradiance on a plane at the top of the atmosphere, at an angle ``Z`` of the sun from directly overhead (the zenith
angle), is

```math
I_\lambda = S_\lambda \left(\frac{a}{r}\right)^2 \cos Z
```

where ``S_\lambda`` is the solar spectral irradiance at one astronomical unit from the sun, ``a`` is the semi-major axis
of the orbit of the earth, and ``r`` is the distance from the earth to the sun. No direct radiation is received when
``Z > 90^\circ``. The spectrum is described in [Data tables](data_tables.md).

## Position of the sun

The distance factor is approximated by (eq. 2)

```math
\left(\frac{a}{r}\right)^2 \approx 1 + 2\epsilon \cos(\omega d)
```

with ``\omega = 2\pi/365`` (see [`orbital_angular_frequency`](@ref)), the eccentricity of the orbit ``\epsilon`` (default = 0.0167238
but this changes over millenia) and the day of the year ``d``. The cosine of the zenith angle is (eq. 3)

```math
\cos Z = \cos\phi \cos\delta \cos h + \sin\phi \sin\delta
```

where ``\phi`` is the latitude, ``\delta`` the solar declination and ``h`` the solar hour angle. The declination, the
the angle between Earth's equatorial plane and a line drawn from the center of the Earth to the center of the sun, 
is (eqs. 4 and 5)

```math
\delta = \arcsin(0.39784993 \sin\zeta), \qquad
\zeta = \omega(d - 80) + 2\epsilon\left[\sin(\omega d) - \sin(80\,\omega)\right]
```

where ``\zeta`` is the ecliptic longitude of the earth in its orbit, and day 80 is the March equinox.
[`solar_geometry`](@ref) calculates these for a given latitude for the [`McCulloughPorterSolarGeometry`](@ref SolarRadiation.McCulloughPorterSolarGeometry) model. For example, for a southern hemisphere latitude (that of Melbourne, Australia) 
at noon in winter:

```@example geometry
model = SolarProblem().solar_geometry_model
geometry = solar_geometry(model, -37.8u"°"; day_of_year = 172, hour_angle = 0.0u"rad")
```

The declination and zenith angles outputs are in radians but the input angles can be entered as degrees or radians. 
At midday on the equator, the seasonal course of the declination and the distance factor are:

```@example geometry
days = 1:365
noon = [solar_geometry(model, 0.0u"°"; day_of_year = d, hour_angle = 0.0u"rad") for d in days]

fig = Figure(size = (700, 320))
ax1 = Axis(fig[1, 1]; xlabel = "Day of the year", ylabel = "Declination (°)")
lines!(ax1, days, rad2deg.(getproperty.(noon, :solar_declination)); linewidth = 2)
ax2 = Axis(fig[1, 2]; xlabel = "Day of the year", ylabel = "Sun distance factor")
lines!(ax2, days, getproperty.(noon, :sun_distance_factor); linewidth = 2)
fig
```

## Solar time and the hour angle

The hour angle is the angle of the sun from solar noon, 15° per hour (eq. 6):

```math
h = \frac{\pi}{12}(t - t_{sn})
```

where ``t`` is the solar time in hours, and the time of solar noon ``t_{sn}`` is 12 plus a longitude correction, in hours, that is
4 minutes for each degree of longitude from the meridian of the time zone. [`hour_angle`](@ref) returns
the angle and the time of solar noon:

```@example geometry
hour_angle(14.0, 0.5)
```

The daily course of the zenith angle at a mid-latitude southern hemisphere site on the solstices and equinox is:

```@example geometry
hours = 0:0.1:24
fig, ax = figure_axis("Solar time (h)", "Zenith angle (°)"; yreversed = true)
for (day, label) in ((355, "December solstice"), (80, "March equinox"), (172, "June solstice"))
    zenith = [solar_geometry(model, -37.8u"°"; day_of_year = day, hour_angle = hour_angle(t)[1]).zenith_angle for t in hours]
    lines!(ax, hours, rad2deg.(zenith); linewidth = 2, label)
end
hlines!(ax, [90]; color = :gray, linestyle = :dash)
axislegend(ax; position = :cb)
fig
```

The dashed line is the horizon: the sun is up when the zenith angle is below 90°.

## Sunrise, sunset and day length

The sun rises and sets at the hour angles ``\pm H`` where the zenith angle is 90° (eq. 7):

```math
H = \arccos(-\tan\delta\tan\phi)
```

Where ``|\tan\delta\tan\phi| \ge 1`` there is no sunrise or sunset. Either the sun stays above the horizon all day
(polar day, ``H = \pi``) or below it (polar night, ``H = 0``), at latitudes poleward of about 66.5° near the solstices.
[`SolarRadiation.sunrise_hour_angle`](@ref) returns ``H``, and the time of sunrise before solar noon in hours, ``H_-``, which is
also in the output of [`solar_radiation`](@ref) as `hour_angle_sunrise`. The day is ``2H_-`` hours long.

```@example geometry
lats = -90:0.5:90
fig, ax = figure_axis("Latitude (°)", "Day length (h)"; yticks = 0:4:24)
for (day, label) in ((172, "June solstice"), (80, "March equinox"), (355, "December solstice"))
    δ = solar_geometry(model, 0.0u"°"; day_of_year = day, hour_angle = 0.0u"rad").solar_declination
    daylength = [2 * SolarRadiation.sunrise_hour_angle(δ, lat * u"°").H₋ for lat in lats]
    lines!(ax, lats, daylength; linewidth = 2, label)
end
axislegend(ax; position = :cb)
fig
```

At the equinox the day is 12 hours long everywhere. The flat tops are polar day and night. The
[polar examples](latitude_examples.md) show what this means for the radiation.

## Azimuth

The azimuth of the sun, measured clockwise from north, is obtained from

```math
\tan Az = \frac{\sin h}{\cos\phi\tan\delta - \sin\phi\cos h}
```

with the quadrant determined by the sign of the hour angle: [`SolarRadiation.solar_azimuth_angle`](@ref) returns
0 to 360°.

## Refraction and air mass

The atmosphere bends the light of a sun near the horizon, so that its apparent zenith angle ``Z_a`` is smaller than the
true one. The correction is applied for ``Z \ge 88^\circ`` in [`SolarRadiation.refraction_correction`](@ref). The refraction is
``R = 16' + 7.5'(Z - 88^\circ)`` where McCullough and Porter (1971) give ``10'`` per degree, so the code corrects slightly
less at low sun angles than the paper. The relative optical air mass is (Rozenberg 1966)

```math
m(Z_a) = \left[\cos Z_a + 0.025 \exp(-11 \cos Z_a)\right]^{-1}
```

which is close to ``\sec Z_a`` for a high sun and finite at the horizon, where it is 40.

```@example geometry
zenith = range(0, 90, length = 361)
airmass = [SolarRadiation.optical_air_mass(deg2rad(z)) for z in zenith]

fig, ax = figure_axis("Zenith angle (°)", "Relative air mass"; yscale = log10)
lines!(ax, zenith, airmass; linewidth = 2, label = "Rozenberg")
lines!(ax, zenith[1:end-20], 1 ./ cosd.(zenith[1:end-20]); linewidth = 2, linestyle = :dash, label = "sec Z")
axislegend(ax; position = :lt)
fig
```

## Twilight

For zenith angles from 88° to 107°, radiation is scattered to the ground although the sun is below the horizon. This
skylight is the only radiation received, and is calculated from a regression of the illuminance of the twilight sky
(Rozenberg 1966), converted to irradiance with 0.0146 W m⁻² per lux (Diem 1966):

```math
G = D = 0.0146 \times 10^{41.34615384 - 0.423076923 Z}
```

with ``Z`` in degrees. It is returned by [`SolarRadiation.twilight_irradiance`](@ref).

```@example geometry
zenith = 85:0.1:110
twilight = [something(SolarRadiation.twilight_irradiance(z * u"°"), 0.0u"W/m^2") for z in zenith]

fig, ax = figure_axis("Zenith angle (°)", "Twilight irradiance (W m⁻²)"; yscale = log10)
lines!(ax, zenith[twilight .> 0u"W/m^2"], ustrip.(twilight[twilight .> 0u"W/m^2"]); linewidth = 2)
fig
```

The skylight is about 27 W m⁻² when the sun is at the horizon and falls by a factor of ten every 2.4°, so it is
less than 0.01 W m⁻² when the sun is more than 7° below it.
