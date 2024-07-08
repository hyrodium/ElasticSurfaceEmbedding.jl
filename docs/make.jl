using ElasticSurfaceEmbedding
using Documenter
using DemoCards



# Create demo with DemoCards.jl
gallery_demopage, gallery_cb, gallery_assets = makedemos(joinpath("gallery"))

# Standard Documenter.jl process
DocMeta.setdocmeta!(ElasticSurfaceEmbedding, :DocTestSetup, :(using ElasticSurfaceEmbedding); recursive = true)
makedocs(;
    modules = [ElasticSurfaceEmbedding],
    authors = "hyrodium <hyrodium@gmail.com> and contributors",
    repo = "https://github.com/hyrodium/ElasticSurfaceEmbedding.jl/blob/{commit}{path}#{line}",
    sitename = "ElasticSurfaceEmbedding.jl",
    format = Documenter.HTML(;
        prettyurls = true,
        canonical = "https://hyrodium.github.io/ElasticSurfaceEmbedding.jl",
        assets = ["assets/favicon.ico", "assets/custom.css", gallery_assets],
        repolink = "https://github.com/hyrodium/ElasticSurfaceEmbedding.jl"
    ),
    pages = [
        "Home" => "index.md",
        "Craft" => "craft.md",
        "Numerical computation" => "run-julia.md",
        "Symbolic computation" => "run-wolfram.md",
        "Gallery" => gallery_demopage,
        "API" => "api.md",
    ],
)

# Postprocess for demos
gallery_cb()

# Deploy docs
deploydocs(; repo = "github.com/hyrodium/ElasticSurfaceEmbedding.jl")
