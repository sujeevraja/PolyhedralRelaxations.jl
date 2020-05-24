using Documenter, PolyhedralRelaxations

makedocs(
    sitename = "PolyhedralRelaxations",

    format = Documenter.HTML(
        # See https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    strict = true,
    authors = "Sujeevraja Sanjeevi and Kaarthik Sundar",
    pages = [
        "Introduction" => "index.md",
        "Quick Start Guide" => "quickstart.md",
    ],
)

# deploydocs(
#     repo   = "github.com/JuliaOpt/JuMP.jl.git",
# )
