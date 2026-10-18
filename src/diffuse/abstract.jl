"""
    AbstractDiffuseModel

Supertype of diffuse (scattered) irradiance models: [`NoScattering`](@ref),
[`DaveFurukawaScattering`](@ref), [`ChandrasekharScattering`](@ref).
"""
abstract type AbstractDiffuseModel end

"""
    diffuse_irradiance(model::AbstractDiffuseModel, wavelength_index, rayleigh_optical_depth, params::SpectralParams, buffers)

Diffuse spectral irradiance at `wavelength_index` given the Rayleigh optical depth at that wavelength.
"""
function diffuse_irradiance end

"""
    allocate_buffers(nmax, diffuse_model)

Working arrays for `solar_radiation!`, reusable across calls. Models needing solver
state add it to the returned NamedTuple.
"""
allocate_buffers(nmax, ::AbstractDiffuseModel) = allocate_spectral_buffers(nmax)

function allocate_spectral_buffers(nmax)
    global_integral = fill(0.0u"W/m^2", nmax)
    rayleigh_integral = fill(0.0u"W/m^2", nmax)
    direct_integral = fill(0.0u"W/m^2", nmax)
    diffuse_integral = fill(0.0u"W/m^2", nmax)
    global_spectrum = global_integral * u"1/nm"
    rayleigh_spectrum = global_integral * u"1/nm"
    direct_spectrum = global_integral * u"1/nm"
    diffuse_spectrum = global_integral * u"1/nm"
    return (; global_integral, rayleigh_integral, direct_integral, diffuse_integral,
             global_spectrum, rayleigh_spectrum, direct_spectrum, diffuse_spectrum)
end
