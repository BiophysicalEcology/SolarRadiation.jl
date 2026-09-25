# Performance

The calculation for a step has to be made at each of 111 wavelengths, and with [`ChandrasekharScattering`](@ref) the diffuse model
dominates the time. For a few places the calculation is fast enough as it is. For many places or times, such as the cells of
a map, it is worth reusing memory.

```@setup performance
using Main.FigureHelpers
using CairoMakie, SolarRadiation, Unitful, Test

function site(latitude)
    SolarTerrain(;
        elevation = 0.0u"m", albedo = 0.2, latitude = latitude * u"°", longitude = 0.0u"°",
        horizon_angles = fill(0.0u"°", 24), slope = 0.0u"°", aspect = 0.0u"°", atmospheric_pressure = 101325.0u"Pa",
    )
end
```

## Reusing memory

[`solar_radiation`](@ref) allocates the arrays of the results and the working arrays for every call. [`solar_radiation!`](@ref) writes
into arrays that were made with [`allocate_output_arrays`](@ref) and [`allocate_buffers`](@ref), so they can be used for any
number of calls with the same number of steps. The buffers depend on the diffuse model, because
[`ChandrasekharScattering`](@ref) needs arrays for its iterations.

```@example performance
model = SolarProblem()
days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349] .* 1.0 # a vector of numbers
hours = collect(0.0:1.0:23.0)
nsteps = length(days) * length(hours)
nmax = model.wavelength_count

out = allocate_output_arrays(nsteps, length(days), nmax)
buffers = allocate_buffers(nmax, model.diffuse_model)

solar_radiation!(out, buffers, model; solar_terrain = site(-37.8), days, hours)
out.global_horizontal[1:3]
```

`solar_radiation!` returns `out`, which the next call overwrites, so results that are needed later have to be copied.
The buffers and the arrays of one call must not be used by another at the same time.

## Speed and allocations

The type of the results is inferred by the compiler, and the calculation allocates little. The memory allocated by a call,
in bytes, does not depend on the number of steps:

```@example performance
terrain = site(-37.8)
allocated(out, buffers, model, terrain, days, hours) = @allocated solar_radiation!(out, buffers, model; solar_terrain = terrain, days, hours)
allocated(out, buffers, model, terrain, days, hours) # compile

long_days = collect(1.0:365.0)
long_out = allocate_output_arrays(length(long_days) * length(hours), length(long_days), nmax)
markdown_table(["Steps", "Allocated (bytes)"],
    [(nsteps, allocated(out, buffers, model, terrain, days, hours)),
     (length(long_days) * length(hours), allocated(long_out, buffers, model, terrain, long_days, hours))])
```

and the compiler infers the result, which can be tested with `@inferred`:

```@example performance
@inferred(solar_radiation!(out, buffers, model; solar_terrain = terrain, days, hours)) isa NamedTuple
```

The time per step, in microseconds, of the three diffuse models:

```@example performance
function microseconds_per_step(diffuse_model)
    model = SolarProblem(; diffuse_model)
    buffers = allocate_buffers(nmax, diffuse_model)
    solar_radiation!(out, buffers, model; solar_terrain = terrain, days, hours) # compile
    1e6 * @elapsed(solar_radiation!(out, buffers, model; solar_terrain = terrain, days, hours)) / nsteps
end
markdown_table(["Diffuse model", "Microseconds per step"],
    [(nameof(typeof(diffuse_model)), round(microseconds_per_step(diffuse_model); sigdigits = 2))
        for diffuse_model in (NoScattering(), DaveFurukawaScattering(), ChandrasekharScattering())])
```

Steps when the sun is below the horizon are cheap, so the times depend on the hours and the latitude.

## Many sites

The work of a map is to do this for each cell. The same arrays can be used for each, one after the other:

```@example performance
latitudes = -80.0:10.0:80.0
function annual_radiation(latitudes)
    map(latitudes) do latitude
        solar_radiation!(long_out, buffers, model; solar_terrain = site(latitude), days = long_days, hours)
        uconvert(u"MJ/m^2", sum(long_out.global_horizontal) * 1u"hr") # per year
    end
end
annual_radiation(latitudes[1:1]) # compile
seconds = @elapsed totals = annual_radiation(latitudes)

fig, ax = figure_axis("Latitude (°)", "Annual global irradiance (MJ m⁻²)"; size = (700, 400))
lines!(ax, latitudes, ustrip.(totals); linewidth = 2)
fig
```

```@example performance
markdown_table(["Sites", "Seconds", "Seconds per site"], [(length(latitudes), round(seconds; sigdigits = 2), round(seconds / length(latitudes); sigdigits = 2))])
```

The buffers are modified by the calculation, so a threaded calculation needs a set of arrays and buffers for each thread. For
example with `Threads.@spawn` for chunks of the sites:

```julia
chunks = Iterators.partition(eachindex(sites), cld(length(sites), Threads.nthreads()))
tasks = map(chunks) do chunk
    Threads.@spawn begin
        out = allocate_output_arrays(nsteps, ndays, nmax)
        buffers = allocate_buffers(nmax, model.diffuse_model)
        map(chunk) do i
            solar_radiation!(out, buffers, model; solar_terrain = terrains[i], days, hours)
            sum(out.global_horizontal) # reduce in the task, as `out` is reused
        end
    end
end
results = reduce(vcat, fetch.(tasks))
```
