using ElasticSurfaceEmbedding
using Documenter

DocMeta.setdocmeta!(ElasticSurfaceEmbedding, :DocTestSetup, :(using ElasticSurfaceEmbedding); recursive=true)

makedocs(;
    modules=[ElasticSurfaceEmbedding],
    authors="hyrodium <hyrodium@gmail.com> and contributors",
    repo="https://github.com/hyrodium/ElasticSurfaceEmbedding.jl/blob/{commit}{path}#{line}",
    sitename="ElasticSurfaceEmbedding.jl",
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://hyrodium.github.io/ElasticSurfaceEmbedding.jl",
        assets=["assets/custom.css"],
    ),
    pages=[
        "Home" => "index.md",
        "Craft" => "craft.md",
        "Numerical computation" => "run-julia.md",
        "Symbolic computation" => "run-wolfram.md",
        "Theoretical Framework" => "theory.md",
        "Gallery" => "gallery.md",
        # "Function Reference" => "functionreference.md",
    ],
)

deploydocs(;
    repo="github.com/hyrodium/ElasticSurfaceEmbedding.jl",
)
