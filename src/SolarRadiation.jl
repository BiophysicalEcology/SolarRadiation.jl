module SolarRadiation

using Unitful
using ModelParameters
using SpecialFunctions, StaticArrays, Dates

export SolarProblem, SolarTerrain, SpectralParams
export AbstractDiffuseModel, NoScattering, DaveFurukawaScattering, ChandrasekharScattering

export scattered_radiation
export elevation_correction
export solar_geometry, hour_angle, orbital_angular_frequency
export solar_radiation, solar_radiation!
export allocate_output_arrays, allocate_buffers

include("constants.jl")
include("elevation_correction.jl")
include("diffuse/abstract.jl")
include("diffuse/no_scattering.jl")
include("diffuse/dave_furukawa.jl")
include("diffuse/chandrasekhar.jl")
include("landscape.jl")
include("solar_geometry.jl")
include("solar_radiation.jl")

end
