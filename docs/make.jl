using Documenter
using DocumenterVitepress
using SolarRadiation
using CairoMakie
using Unitful

# Don't output huge svgs for Makie plots
CairoMakie.activate!(type = "png")

# RasterDataSources downloads (GADS, CRU CL2) go here unless the path is set
get!(ENV, "RASTERDATASOURCES_PATH", joinpath(@__DIR__, "data"))
mkpath(ENV["RASTERDATASOURCES_PATH"])

# Helpers for the figures, loaded in the examples with `using Main.FigureHelpers`
include("figure_helpers.jl")

makedocs(
    modules = [SolarRadiation],
    sitename = "SolarRadiation.jl",
    authors = "Michael Kearney, Rafael Schouten et al.",
    clean = true,
    doctest = false,
    checkdocs = :exports,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/BiophysicalEcology/SolarRadiation.jl", # this must be the full URL!
        devbranch = "main",
        devurl = "dev";
    ),
    source = "src",
    build = "build",
    warnonly = true,
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/BiophysicalEcology/SolarRadiation.jl",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true,
)
