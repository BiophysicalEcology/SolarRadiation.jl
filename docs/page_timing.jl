# Optional timing of the pages of the manual, loaded by make.jl when the environment variable DOCS_PAGE_TIMES is set.
# Documenter runs the examples of one page after another. The time of a page is the time since the previous page ended,
# so the first page also includes the start of the build.
module PageTiming

using Documenter

const TIMES = Pair{String,Float64}[]
const LAST = Ref(time())

function Documenter.collect_named_anchors!(page::Documenter.Page, doc::Documenter.Document)
    now = time()
    push!(TIMES, page.source => now - LAST[])
    LAST[] = now
    return invoke(Documenter.collect_named_anchors!, Tuple{Any,Any}, page, doc)
end

"Print the time of each page, longest first, and write it to `page_times.txt`."
function report(file = joinpath(@__DIR__, "page_times.txt"))
    sorted = sort(TIMES; by = last, rev = true)
    lines = ["$(rpad(first(t), 34)) $(lpad(round(last(t); digits = 1), 8)) s" for t in sorted]
    push!(lines, "$(rpad("total", 34)) $(lpad(round(sum(last, TIMES); digits = 1), 8)) s")
    write(file, join(lines, "\n") * "\n")
    println("\nTime of each page:\n", join(lines, "\n"))
end

end
