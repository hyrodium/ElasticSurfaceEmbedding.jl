using ElasticSurfaceEmbedding
using Documenter

DocMeta.setdocmeta!(ElasticSurfaceEmbedding, :DocTestSetup, :(using ElasticSurfaceEmbedding); recursive=true)

makedocs(;
    modules=[ElasticSurfaceEmbedding],
    authors="hyrodium <hyrodium@gmail.com> and contributors",
    repo="https://github.com/hyrodium/ElasticSurfaceEmbedding.jl/blob/{commit}{path}#{line}",
    sitename="ElasticSurfaceEmbedding.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://hyrodium.github.io/ElasticSurfaceEmbedding.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Craft" => "craft.md",
        "Computing" => "computing.md",
        "Theory" => "theory.md",
        "Gallery" => "gallery.md",
        # "Function Reference" => "functionreference.md",
    ],
)

deploydocs(;
    repo="github.com/hyrodium/ElasticSurfaceEmbedding.jl",
)
