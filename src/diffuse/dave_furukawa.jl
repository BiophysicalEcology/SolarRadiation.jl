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
    wavelength_index > 11 && return 0.0u"W/m^2/nm"
    (; sky_irradiance, ground_reflected, single_scattering_albedo) = model
    (; solar_spectral_irradiance, sun_distance_factor, zenith_angle, albedo) = params

    n = wavelength_index
    Fd = sky_irradiance
    Fd′_Q = ground_reflected
    s̄ = single_scattering_albedo
    Sλ = solar_spectral_irradiance
    ar² = sun_distance_factor
    Z = zenith_angle
    A = albedo

    B = uconvert(NoUnits, Z / 5u"°") # zenith angle in steps of 5°, the rows of the tables
    k = trunc(Int, B) + 1 + (B % 1 > 0.5)
    Q = A / (1 - A * s̄[n]) # eq. 31 in Dave & Furukawa (1966)
    Dλ = (Sλ[n] / π) * ar² * (Fd[n, k] + Fd′_Q[n, k] * Q) # eq. 16 in McCullough & Porter (1971)
    return Dλ
end
