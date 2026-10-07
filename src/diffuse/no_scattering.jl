"""
    NoScattering()

No diffuse irradiance.
"""
struct NoScattering <: AbstractDiffuseModel end

diffuse_irradiance(::NoScattering, wavelength_index, rayleigh_optical_depth, params, buffers) = 0.0u"W/m^2/nm"
