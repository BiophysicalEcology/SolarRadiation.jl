"""
    DaveFurukawaScattering(; sky_irradiance, ground_reflected, single_scattering_albedo)

Diffuse irradiance from tabulated values, UV wavelengths only (first 11 intervals).

## Keywords

- `sky_irradiance`: radiation scattered from the direct beam, by wavelength and zenith angle.
- `ground_reflected`: ground-reflected radiation rescattered downward, by wavelength and zenith angle.
- `single_scattering_albedo`: molecular scattering function by wavelength.

## References

Dave, J. V. & Furukawa, P. M. (1966). Scattered radiation in the ozone absorption bands
at selected levels of a terrestrial, Rayleigh atmosphere. Meteorological Monographs 7(29).
"""
@kwdef struct DaveFurukawaScattering{FD,FDQ,SSA} <: AbstractDiffuseModel
    sky_irradiance::FD = DEFAULT_DIFFUSE_SKY_IRRADIANCE
    ground_reflected::FDQ = DEFAULT_DIFFUSE_GROUND_REFLECTED
    single_scattering_albedo::SSA = DEFAULT_SINGLE_SCATTERING_ALBEDO
end

function diffuse_irradiance(model::DaveFurukawaScattering, wavelength_index, rayleigh_optical_depth, params, buffers)
    n = wavelength_index
    n > 11 && return 0.0u"W/m^2/nm"
    ar², A, z, Sλ = params.sun_distance_factor, params.albedo, params.zenith_angle, params.solar_spectral_irradiance
    (; sky_irradiance, ground_reflected, single_scattering_albedo) = model
    B = ustrip(u"°", z) / 5
    k = trunc(Int, B) + 1 + (B % 1 > 0.5)
    Q = A / (1.0 - A * single_scattering_albedo[n])  # eq. 31 in Dave & Furukawa 1966
    Dλ = (Sλ[n] / π) * (sky_irradiance[n, k] + ground_reflected[n, k] * Q) / 1000.0
    return Dλ * ar²
end
