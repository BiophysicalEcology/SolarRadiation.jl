# Validates the water vapour optical depths against RRTMGP band-integrated direct transmission.
# Not part of the test suite. Run with `julia --project=test/RRTMGP test/RRTMGP/water_bands.jl`;
# the RRTMGP lookup tables are downloaded on first use.
#
# RRTMGP `ClearSkyRadiation` (no aerosol) on its `standard_atmosphere` columns, surface albedo 0,
# per-band fluxes. Surface over top-of-atmosphere flux in a band is the direct transmission for
# bands above 625 nm, where Rayleigh diffuse is negligible. SolarRadiation is run with the same
# precipitable water and ozone column, no aerosol, weighted by its own spectrum within each band.
using SolarRadiation, RRTMGP, NCDatasets, Unitful, Printf
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
const SR = SolarRadiation

FT = Float64
g, Md, Mw, NA = 9.81, 0.02897, 0.018015, 6.02214076e23

function band_transmission(prof, μ, lookups)
    nlay, ncol = size(prof.p_lay)
    grid_params = RRTMGP.RRTMGPGridParams(FT; context=ClimaComms.context(), domain_nlay=nlay, ncol)
    vmr_wm = zeros(FT, lookups.ngas_lw)
    for (gas, v) in prof.well_mixed_vmr
        vmr_wm[lookups.idx_gases_lw[gas]] = v
    end
    layerdata = zeros(FT, 4, nlay, ncol)
    layerdata[2, :, :] .= prof.p_lay
    layerdata[3, :, :] .= prof.t_lay
    vmr = RRTMGP.VolumeMixingRatios.VmrGM(copy(prof.vmr_h2o), copy(prof.vmr_o3), vmr_wm)
    as = RRTMGP.AtmosphericStates.AtmosphericState(zeros(FT, ncol), copy(prof.lat), layerdata,
        copy(prof.p_lev), copy(prof.t_lev), copy(prof.t_sfc), vmr, nothing, nothing)
    bcs_lw = RRTMGP.BCs.LwBCs(ones(FT, lookups.nbnd_lw, ncol), nothing)
    albedo = zeros(FT, lookups.nbnd_sw, ncol)
    bcs_sw = RRTMGP.BCs.SwBCs(fill(FT(μ), ncol), fill(FT(1361), ncol), albedo, nothing, copy(albedo))
    s = RRTMGP.RRTMGPSolver(grid_params, RRTMGP.ClearSkyRadiation(false), RRTMGP.default_parameters(FT),
        bcs_lw, bcs_sw, as; lookups, spectral_fluxes=true)
    RRTMGP.update_fluxes!(s)
    F = Array(RRTMGP.spectral_sw_flux_dn(s))
    return F[1, 1, :] ./ F[end, 1, :], Array(RRTMGP.sw_band_bounds(s))
end

# precipitable water (cm) and ozone column (atm-cm) of an RRTMGP profile
function columns(prof)
    Δp = prof.p_lev[1:end-1, 1] .- prof.p_lev[2:end, 1]
    r = prof.vmr_h2o[:, 1] .* (Mw / Md)
    return sum(r ./ (1 .+ r) .* Δp) / g / 10, sum(prof.vmr_o3[:, 1] .* Δp) / (g * Md) * NA / 2.687e20 / 1000
end

λ = ustrip.(u"nm", SR.DEFAULT_WAVELENGTHS)
Sλ = ustrip.(SR.DEFAULT_SOLAR_SPECTRAL_IRRADIANCE)
w = zeros(length(λ)) # trapezoid weights
for i in 2:length(λ)
    w[i-1] += (λ[i] - λ[i-1]) / 2
    w[i] += (λ[i] - λ[i-1]) / 2
end

# SolarRadiation direct transmission at each wavelength, no aerosol
function sr_transmission(τW, μ, pw_cm, o3_cm)
    m = SR.optical_air_mass(acos(μ))
    ef = SR.elevation_correction(0.0u"m")
    map(eachindex(λ)) do n
        τ = SR.spectral_optical_depth(n, 101325.0u"Pa", 25.0u"km", SR.DEFAULT_RAYLEIGH_OPTICAL_DEPTH,
            SR.DEFAULT_OZONE_OPTICAL_DEPTH, zeros(length(λ)), τW, SR.DEFAULT_MIXED_GAS_ABSORPTION, ef, o3_cm, m,
            pw_cm * u"cm")
        exp(-τ.total)
    end
end
weighted(T, idx) = sum(w[idx] .* Sλ[idx] .* T[idx]) / sum(w[idx] .* Sλ[idx])

lookups = RRTMGP.solve(RRTMGP.standard_atmosphere(FT; nlay=60)).solver.lookups
cases = map(Iterators.product((:tropical, :midlatitude_summer, :subarctic_winter), (1.0, 0.7, 0.4, 0.2))) do (kind, μ)
    prof = RRTMGP.standard_atmosphere(FT; kind, nlay=60)
    T, bounds = band_transmission(prof, μ, lookups)
    (; kind, μ, cols=columns(prof), T, bounds)
end |> vec
bounds = first(cases).bounds
nm_lo, nm_hi = 1e7 ./ bounds[2, :], 1e7 ./ bounds[1, :]
bands = [(b, findall(i -> nm_lo[b] <= λ[i] < nm_hi[b], eachindex(λ))) for b in axes(bounds, 2) if 625 <= nm_lo[b] < 3700]
nir = reduce(vcat, last.(bands))

τ_opaque = Float64.(SR.DEFAULT_WATER_OPTICAL_DEPTH)
opaque = findall(==(SR.MAX_OPTICAL_DEPTH), τ_opaque)
τ_transparent = copy(τ_opaque); τ_transparent[opaque] .= 0
τ_half = copy(τ_opaque); τ_half[opaque] .= SR.MAX_OPTICAL_DEPTH / 2

for (label, τW) in ("band centres transparent" => τ_transparent, "band centres opaque (default)" => τ_opaque)
    println("\n$label: band direct transmission, SolarRadiation − RRTMGP")
    for (b, idx) in bands
        @printf("%5.0f–%5.0f nm |", nm_lo[b], nm_hi[b])
        for c in cases
            @printf(" %+.3f", weighted(sr_transmission(τW, c.μ, c.cols...), idx) - c.T[b])
        end
        println()
    end
end

println("\n625–3700 nm direct transmission relative to RRTMGP, and change when opaque is halved")
for c in cases
    Trr = sum(sum(w[idx] .* Sλ[idx]) * c.T[b] for (b, idx) in bands) / sum(w[nir] .* Sλ[nir])
    rel(τW) = 100 * (weighted(sr_transmission(τW, c.μ, c.cols...), nir) / Trr - 1)
    @printf("  %-19s μ=%.1f  transparent %+5.1f%%  opaque %+5.1f%%  halved %+.1e%%\n",
            c.kind, c.μ, rel(τ_transparent), rel(τ_opaque), rel(τ_half) - rel(τ_opaque))
end
