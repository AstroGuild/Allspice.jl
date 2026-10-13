# Write the reference pages as plain Markdown with each `@docs` entry replaced by its
# docstring, for a reader review that has no access to the source.
#
#     julia --project=docs docs/reference_text.jl <output directory>

using Allspice, Markdown, REPL

function reference_text(outdir::AbstractString)
    src = joinpath(@__DIR__, "src", "reference")
    mkpath(outdir)
    for page in filter(endswith(".md"), readdir(src))
        io = IOBuffer()
        inblock = false
        for line in eachline(joinpath(src, page))
            if startswith(line, "```@docs")
                inblock = true
            elseif inblock && startswith(line, "```")
                inblock = false
            elseif inblock
                name = strip(line)
                isempty(name) && continue
                binding = Docs.Binding(Allspice, Symbol(replace(name, "Allspice." => "")))
                println(io, "---\n\n### `", name, "`\n\n", Markdown.plain(Docs.doc(binding)), "\n")
            else
                println(io, line)
            end
        end
        write(joinpath(outdir, page), take!(io))
    end
    return outdir
end

reference_text(only(ARGS))
