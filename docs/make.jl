using Documenter, PolyhedralRelaxations
using DocumenterTools: Themes

Themes.compile(
    joinpath(@__DIR__, "src/assets/pr-light.scss"),
    joinpath(@__DIR__, "src/assets/themes/documenter-light.css"),
)

Themes.compile(
    joinpath(@__DIR__, "src/assets/pr-dark.scss"),
    joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"),
)

makedocs(
    sitename = "PolyhedralRelaxations",
    format = Documenter.HTML(
        mathengine = Documenter.MathJax(),
        assets=[asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css)],
        prettyurls = get(ENV, "CI", nothing) == "true",
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
