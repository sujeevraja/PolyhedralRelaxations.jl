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
        "API" => "api.md",
        "Interfacing with JuMP" => "interface.md",
        "Reference" => "reference.md"
    ],
)

deploydocs(
    repo   = "github.com/sujeevraja/PolyhedralRelaxations.jl.git",
)
