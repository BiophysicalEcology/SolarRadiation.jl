# Dates and calendars

The times of the calculation can be given as numbers, the days of the year and the hours of solar time, or as dates and
hours from the `Dates` standard library.

```@setup dates
using Main.FigureHelpers
using CairoMakie, SolarRadiation, Unitful, Dates

terrain = SolarTerrain(;
    elevation = 0.0u"m", albedo = 0.2, latitude = -37.8u"°", longitude = 144.96u"°",
    horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°", atmospheric_pressure = 101325.0u"Pa",
)
model = SolarProblem()
```

## Numbers

With the numeric interface the `days` are days of the year and the `hours` are hours of solar time. The `year` is used only to find
the month for the [ozone column](atmosphere.md#Ozone), and `days_in_year` is 365 unless it is changed:

```@example dates
out = solar_radiation(model; solar_terrain = terrain, days = [15, 196], hours = 6:2:18, year = 2024, days_in_year = 366)
markdown_table(["Day of the year", "Hour", "Global irradiance (W m⁻²)"],
    [(out.day_of_year[i], out.hour[i], round(out.global_horizontal[i] / u"W/m^2"; digits = 1)) for i in eachindex(out.hour)])
```

The defaults are the middle days of the twelve months, and every hour of the day.

## Dates

With `dates` the day of the year, the year and the length of the year come from the date, and the `hours` are `Period`s:

```@example dates
dates = [Date(2024, 1, 15), Date(2024, 7, 15)]
out = solar_radiation(model; solar_terrain = terrain, dates, hours = Hour(6):Hour(2):Hour(18))
markdown_table(["Day of the year", "Hour", "Global irradiance (W m⁻²)"],
    [(out.day_of_year[i], out.hour[i], round(out.global_horizontal[i] / u"W/m^2"; digits = 1)) for i in eachindex(out.hour)])
```

`Date` and `DateTime` are both accepted. The hours must be periods, such as `Hour(12)` or `Minute(30)`. Bare numbers are read as
milliseconds, so that with `hours = 0:23` all the times are within the first few milliseconds of midnight:

```@example dates
out = solar_radiation(model; solar_terrain = terrain, dates = [Date(2024, 1, 15)], hours = 0:3)
markdown_table(["Step", "Hour"], [(i, out.hour[i]) for i in eachindex(out.hour)])
```

### Leap years

The year length of the first date is used, so the days of the year and the position of the earth in its orbit follow the
calendar. The 21st of June is day 173 of a leap year, and 172 in other years:

```@example dates
markdown_table(["Year", "Day of the year of 21 June"], [(year, Dates.dayofyear(Date(year, 6, 21))) for year in (2023, 2024)])
```

Calendars with other year lengths, such as the 360-day and no-leap calendars of climate models, can be used with the
date types of the CFTime.jl package, for example `DateTimeNoLeap`, and the length of the year comes from the date in the same way.

## Solar time and time zones

The hours are solar time unless a longitude correction is given. Solar noon is then at 12:00, and the sun is highest at that time.
`longitude_correction` (or `timezone_offset` with the `dates`) is the number of hours added to 12:00 to give the time of solar noon
on the clock. It is 4 minutes for each degree of longitude between the site and the meridian of the time zone. It is not the offset from UTC.

Melbourne is at 144.96°E in a zone of UTC+10, whose standard meridian is 150°E. Solar noon on the clock is then about 12:20:

```@example dates
correction = (150.0 - 144.96) / 15 # hours
clock = solar_radiation(model; solar_terrain = terrain, days = [15], hours = 0:0.25:23.75, longitude_correction = correction)
peak = clock.hour[argmax(clock.global_horizontal)]

markdown_table(["Quantity", "Hours"],
    [("Longitude correction", round(correction; digits = 2)),
     ("Solar noon on the clock", round(clock.hour_solar_noon[1]; digits = 2)),
     ("Time of the greatest irradiance", round(peak; digits = 2))])
```

The clock time of the greatest irradiance is that of solar noon, within the 15 minute time step. Daylight saving time,
if it applies, is added to the clock hours by the user.
