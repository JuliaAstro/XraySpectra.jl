using Documenter
using Documenter.Remotes: GitHub
using XraySpectra

makedocs(;
    modules = [XraySpectra],
    checkdocs = :exports,
    sitename = "XraySpectra.jl",
    authors = "James Lee and contributors",
    repo = GitHub("JuliaAstro/XraySpectra.jl"),
    format = Documenter.HTML(;
        canonical = "https://juliaastro.org/XraySpectra/stable/",
        edit_link = "main",
        warn_outdated = true,
    ),
    pages = [
        "Home" => "index.md",
        "Loading OGIP data" => "loading.md",
        "Responses" => "responses.md",
        "Folding" => "folding.md",
        "Channel rebinning" => "rebinning.md",
        "Building the documentation" => "development.md",
        "API reference" => "reference.md",
    ],
)

deploydocs(;
    repo = "github.com/JuliaAstro/XraySpectra.jl.git",
    devbranch = "main",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"],
)
