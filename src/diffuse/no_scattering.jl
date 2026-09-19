"""
    NoScattering()

No diffuse irradiance.
"""
struct NoScattering <: AbstractDiffuseModel end

diffuse_irradiance(::NoScattering, n, λτR, params, buffers) = 0.0u"W/m^2/nm"
