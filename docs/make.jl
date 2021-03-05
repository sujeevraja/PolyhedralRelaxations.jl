using Documenter, PolyhedralRelaxations
using DocumenterTools: Themes

Themes.compile(joinpath(@__DIR__,"src/assets/pr.scss"), joinpath(@__DIR__,"src/assets/themes/documenter-light.css"))

makedocs(
    sitename = "PolyhedralRelaxations",
    format = Documenter.HTML(
        mathengine = Documenter.MathJax(),
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    strict = true,
    authors = "Sujeevraja Sanjeevi and Kaarthik Sundar",
    pages = [
        "Introduction" => "index.md",
        "API" => "api.md",
        "Interfacing with JuMP" => "interface.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(repo = "github.com/sujeevraja/PolyhedralRelaxations.jl.git")
