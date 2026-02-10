module SolarRadiation

using Unitful
using ModelParameters
using SpecialFunctions, StaticArrays, Dates

export SolarProblem, SolarTerrain

export scattered_radiation
export elevation_correction
export solar_geometry, hour_angle
export solar_radiation, solar_radiation!
export allocate_output_arrays, allocate_buffers

include("constants.jl")
include("elevation_correction.jl")
include("landscape.jl")
include("scattered_uv.jl")
include("solar_geometry.jl")
include("solar_radiation.jl")

end
