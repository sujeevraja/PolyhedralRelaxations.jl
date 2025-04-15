using Documenter, PolyhedralRelaxations

makedocs(
    sitename = "PolyhedralRelaxations",
    format = Documenter.HTML(
        mathengine = Documenter.MathJax(),
        assets = [
            asset(
                "https://fonts.googleapis.com/css?family=Fira+Sans|Source+Code+Pro&display=swap",
                class = :css,
            ),
        ],
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    authors = "Sujeevraja Sanjeevi and Kaarthik Sundar",
    pages = [
        "Introduction" => "index.md",
        "API" => "api.md",
        "Interfacing with JuMP" => "interface.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(repo = "github.com/sujeevraja/PolyhedralRelaxations.jl.git")
