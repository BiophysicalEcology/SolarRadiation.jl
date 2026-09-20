"""
    DaveFurukawaScattering(; sky_irradiance, ground_reflected, spherical_albedo)

Diffuse irradiance from tabulated values, UV wavelengths only (first 11 intervals).
The tables are for a sea-level surface (1000 mb) and a total ozone column of 0.34 cm.

# Keywords

- `sky_irradiance`: ``F_d``, the downward flux of scattered radiation, of all orders of scattering, due to the direct solar beam,
  by wavelength (11) and zenith angle (19, from 0° to 90° in steps of 5°), for an incident flux of ``\\pi``, so in units of ``S_\\lambda/\\pi``.
- `ground_reflected`: ``F_d'/Q``, the quantity that gives the contribution to ``F_d`` from illumination of the atmosphere from
  below by ground-reflected radiation, with ``Q = A/(1 - A\\bar s)``. Same layout as `sky_irradiance`.
- `spherical_albedo`: ``\\bar s``, by wavelength. It is ``S^b`` of Table B of Dave & Furukawa (1966), the fraction of the
  ground-reflected radiation that the atmosphere scatters back to the ground, the spherical albedo of the atmosphere.

# References

Dave, J. V. & Furukawa, P. M. (1966). Scattered radiation in the ozone absorption bands
at selected levels of a terrestrial, Rayleigh atmosphere. Meteorological Monographs 7(29).
"""
@kwdef struct DaveFurukawaScattering{FD,FDQ,SA} <: AbstractDiffuseModel
    sky_irradiance::FD = DEFAULT_DIFFUSE_SKY_IRRADIANCE
    ground_reflected::FDQ = DEFAULT_DIFFUSE_GROUND_REFLECTED
    spherical_albedo::SA = DEFAULT_SPHERICAL_ALBEDO
end

function diffuse_irradiance(model::DaveFurukawaScattering, wavelength_index, rayleigh_optical_depth, params, buffers)
    n = wavelength_index
    n > 11 && return 0.0u"W/m^2/nm"
    ar², A, z, Sλ = params.sun_distance_factor, params.albedo, params.zenith_angle, params.solar_spectral_irradiance
    (; sky_irradiance, ground_reflected, spherical_albedo) = model
    B = ustrip(u"°", z) / 5
    k = trunc(Int, B) + 1 + (B % 1 > 0.5)
    Q = A / (1.0 - A * spherical_albedo[n])  # eq. 31 in Dave & Furukawa 1966
    Dλ = (Sλ[n] / π) * (sky_irradiance[n, k] + ground_reflected[n, k] * Q) / 1000.0
    return Dλ * ar²
end
