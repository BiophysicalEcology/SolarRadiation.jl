module FigureHelpers

using CairoMakie
using GeoInterface
using Markdown
using NaturalEarth: naturalearth

export figure_axis, countries!, markdown_table

"""
    figure_axis(xlabel, ylabel; size=(700, 500), kw...)

A `Figure` and `Axis` with minor ticks and grid lines, as used for the figures of the manual.
"""
function figure_axis(xlabel, ylabel; size=(700, 500), kw...)
    fig = Figure(; size)
    ax = Axis(fig[1, 1];
        xlabel, ylabel,
        xminorticksvisible=true, yminorticksvisible=true,
        xminorgridvisible=true, yminorgridvisible=true,
        xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5),
        kw...,
    )
    return fig, ax
end

"""
    markdown_table(header, rows)

A Markdown table with the column names `header` and one row for each element of `rows`, which has one value for each column.
The first column is aligned to the left and the others to the right.
"""
function markdown_table(header, rows)
    line(cells) = "| " * join(string.(cells), " | ") * " |"
    alignment = "| :--- | " * join(fill("---:", length(header) - 1), " | ") * " |"
    return Markdown.parse(join([line(header); alignment; [line(row) for row in rows]], Char(10)))
end

const COUNTRY_OUTLINES = Ref{Vector{Point2f}}()

# Points of every ring or line, with NaN between them so that one `lines!` call draws them all
function _outline_points!(points, geometry)
    trait = GeoInterface.geomtrait(geometry)
    if trait isa Union{GeoInterface.LineStringTrait,GeoInterface.LinearRingTrait}
        for p in GeoInterface.getpoint(geometry)
            GeoInterface.y(p) > -89.9 ? push!(points, Point2f(GeoInterface.x(p), GeoInterface.y(p))) : push!(points, Point2f(NaN, NaN)) # not the edge of Antarctica at the pole
        end
        push!(points, Point2f(NaN, NaN))
    elseif trait isa GeoInterface.PolygonTrait
        foreach(ring -> _outline_points!(points, ring), GeoInterface.getring(geometry))
    else
        foreach(part -> _outline_points!(points, part), GeoInterface.getgeom(geometry))
    end
    return points
end

"""
    countries!(ax; color=(:black, 0.6), linewidth=0.6)

Draw the outlines of the countries, and so the coasts of the continents, from the Natural Earth 110 m
data, on the axis `ax` of a map in degrees of longitude and latitude.
"""
function countries!(ax; color=(:black, 0.6), linewidth=0.6)
    if !isassigned(COUNTRY_OUTLINES)
        points = Point2f[]
        foreach(f -> _outline_points!(points, GeoInterface.geometry(f)), naturalearth("admin_0_countries", 110))
        COUNTRY_OUTLINES[] = points
    end
    return lines!(ax, COUNTRY_OUTLINES[]; color, linewidth)
end

end
