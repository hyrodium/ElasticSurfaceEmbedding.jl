# Default output directory
OUT_DIR = joinpath(homedir(),"ElasticSurfaceEmbedding-Result")

"""
    config_dir(dir)

Set the output directory.
The default is `~/ElasticSurfaceEmbedding-Result`.
"""
function config_dir(dir)
    _dir = expanduser(dir)
    mkpath(_dir)
    global OUT_DIR = _dir
end

function _check_filename(name::String)
    invalidcharacters = [' ', '&', '\\', '/', '.', '<', '>', '|', ':', ';', '*', '?', '=', '%', '$', '"', '~']
    if length(invalidcharacters ∩ name) ≠ 0
        error("The name[$(name)] must not consists of following chars: ", invalidcharacters)
    end
end

"""
    settings(name::String; canvas=(10,10), mesh=(10, 1), unit=100, maximumstrain=0.0, colorbarsize=0.2)

Initial settings for the given shape, or load results if the computed piecies of surface exists.
"""
function settings(
    name::String;
    canvas::Tuple{<:Real,<:Real}=(10,10),
    mesh::Tuple{Int,Int}=(10, 1),
    unit::Real=100,
    maximumstrain::Real=0.0,
    colorbarsize::Real=0.2,
)
    _check_filename(name)
    global NAME = name
    global DIR = joinpath(OUT_DIR, NAME)
    global YLIMS = (-canvas[2]/2, canvas[2]/2)
    global XLIMS = (-canvas[1]/2, canvas[1]/2)
    global MESH = mesh
    global UNIT = (unit, "pt")
    global MAXIMUMSTRAIN = maximumstrain
    global COLORBARSIZE = colorbarsize
    if isTheShapeComputed()
        dict = LoadResultDict()
        println(_tree_as_string(dict["result"]))
    end
end

function toJSON(k::KnotVector)
    return k.vector
end
function toJSON(P::BSplineSpace)
    return Dict("degree" => degree(P), "knotvector" => toJSON(knotvector(P)))
end
function toJSON(P::Array{P,1} where {P<:BSplineSpace})
    return convert(Array{Any,1}, toJSON.(P))
end
function toJSON(a::Array{Float64})
    return convert(Array{Any,1}, a[:])
end
function toJSON(M::BSplineManifold{2})
    return Dict("bsplinespaces" => toJSON(collect(bsplinespaces(M))), "controlpoints" => toJSON(controlpoints(M)))
end

function JSONtoBSplineSpace(jP::Dict)
    p = jP["degree"]
    return BSplineSpace{p}(KnotVector(Vector{Float64}(jP["knotvector"])))
end
function JSONtoBSplineSpaces(jPs::Array)
    P1 = JSONtoBSplineSpace(jPs[1])
    P2 = JSONtoBSplineSpace(jPs[2])
    return (P1,P2)
end
function JSONtoControlPoints(ja, dims)
    n = prod(dims)
    d = length(ja) ÷ n
    return reshape(convert(Array{Float64}, ja), dims..., d)
end
function JSONtoBSplineManifold(dict::Dict)
    P = JSONtoBSplineSpaces(dict["bsplinespaces"])
    dims = dim.(P)
    𝒂 = JSONtoControlPoints(dict["controlpoints"], dims)
    return BSplineManifold(𝒂, P)
end

function NodeSeries(tree::Dict, node)
    Nodes = [node]
    while Nodes[end] ≠ "0"
        push!(Nodes, tree[Nodes[end]]["parent"])
    end
    return Nodes
end

function _tree_as_string(tree::Dict)
    serieses = Array{Int,1}[]
    for key in keys(tree)
        push!(serieses, (s -> parse(Int, s)).(reverse(NodeSeries(tree, key))))
    end
    sort!(serieses)
    lowstrings = String[]
    n = length(serieses)
    for i in 1:n
        l = length(serieses[i])
        key = string(serieses[i][end])
        comment = tree[key]["comment"]
        if l == 2
            lowstring = key * ": " * comment
            push!(lowstrings, lowstring)
        elseif l ≥ 3
            lowstring = "  "^(l - 3) * "└─" * key * ": " * comment
            push!(lowstrings, lowstring)
            for j in 1:(i-1)
                chars = collect(lowstrings[end-j])
                if chars[2(l-3)+1] == ' '
                    lowstrings[end-j] = join(chars[1:2(l-3)]) * "│" * join(chars[2(l-3)+2:end])
                elseif chars[2(l-3)+1] == '└'
                    lowstrings[end-j] = join(chars[1:2(l-3)]) * "├" * join(chars[2(l-3)+2:end])
                    break
                else
                    break
                end
            end
        end
    end
    outsting = ""
    for s in lowstrings
        outsting = outsting * s * "\n"
    end
    return outsting
end

function isTheShapeComputed()
    return isfile(joinpath(DIR,NAME*".json"))
end

function latest_index()
    dict = LoadResultDict()
    latest_index(dict)
end

function latest_index(dict)
    result_nums = [parse(Int, i) for i in keys(dict["result"])]
    return maximum(result_nums)
end

function _realparent(index::Int)
    if index==0
        return latest_index()
    else
        return index
    end
end

function loadM(; index=0)
    if !isTheShapeComputed()
        error("Result file doesn't exists")
    end
    dict = LoadResultDict()
    index = _realparent(index)
    M = JSONtoBSplineManifold(dict["result"][string(index)]["BSplineManifold"])
    return M
end

function LoadResultDict()::Dict
    json = read(joinpath(DIR, "$(NAME).json"), String)
    dict = JSON.parse(json)
    v = dict["meta"]["version"]
    JSON_VERSION = VersionNumber(v["major"], v["minor"], v["patch"])

    if JSON_VERSION.major == ESE_VERSION.major
        return dict
    else
        error("The version of computed shape is not supported.")
    end
end

function _export(M::BSplineManifold{2}, parent::Int; comment = "", maximumstrain = MAXIMUMSTRAIN)
    if isTheShapeComputed()
        dict = LoadResultDict()
        index = latest_index(dict) + 1
    else
        meta = Dict{String,Any}("version" => ESE_VERSION, "generated by" => "ElasticSurfaceEmbedding.jl")
        dict = Dict{String,Any}("meta" => meta)
        dict["result"] = Dict{String,Any}()
        index = 1
    end

    dict["result"][string(index)] = Dict{String,Any}("parent" => string(parent))
    dict["result"][string(index)]["BSplineManifold"] = toJSON(M)
    dict["result"][string(index)]["comment"] = comment

    # Save as json
    mkpath(DIR)
    write(joinpath(DIR, NAME*".json"), JSON.json(dict, 4))

    # Compute maximum strain for graphics export
    if maximumstrain == 0.0
        MS = ComputeMaximumStrain(index = index)
        MaximumStrain = max(-MS[1], MS[2])
    else
        MaximumStrain = maximumstrain
    end

    # Save graphics
    ExportFiles(M, MaximumStrain, index)

    # Send messages
    path_png_append = joinpath(DIR, "append", "$(NAME)-$(index)_append.png")
    message = _tree_as_string(dict["result"])
    println(message)
    _send_file_to_slack(path_png_append, comment = "```\n" * message * "```")
end

function ExportFiles(
    M::BSplineManifold{2},
    MaximumStrain::Real,
    index;
    Name::String = NAME,
    Dir = DIR,
    Xlims = XLIMS,
    Ylims = YLIMS,
    Mesh = MESH,
    Unit = UNIT,
    Colorbarsize = COLORBARSIZE,
)
    mkpath(joinpath(DIR, "nurbs"))
    mkpath(joinpath(DIR, "strain"))
    mkpath(joinpath(DIR, "colorbar"))
    mkpath(joinpath(DIR, "append"))

    Width = (Xlims[2] - Xlims[1]) * Unit[1]
    Height = (Ylims[2] - Ylims[1]) * Unit[1]

    normalized_strain(u¹, u²) = E⁽⁰⁾₁₁_cont(M, u¹, u²) / MaximumStrain # bounded in -1 to 1

    aa = 5 # magnification parameter for antialias

    path_svg_nurbs = joinpath(Dir, "nurbs", "$(Name)-$(index)_Bspline.svg")
    path_png_nurbs = joinpath(Dir, "nurbs", "$(Name)-$(index)_Bspline.png")
    path_png_strain = joinpath(Dir, "strain", "$(Name)-$(index)_strain.png")
    path_png_colorbar = joinpath(Dir, "colorbar", "$(Name)-$(index)_colorbar.png")
    path_png_append = joinpath(Dir, "append", "$(Name)-$(index)_append.png")

    colorfunc(u¹,u²) = normalized_strain(u¹,u²) * RGB(0.5, -0.5, -0.5) + RGB(0.5, 0.5, 0.5) # red to cyan

    save_svg(path_svg_nurbs, M, xlims=Xlims, ylims=Ylims, mesh=Mesh, unitlength=Int(Unit[1]))
    save_png(path_png_nurbs, M, xlims=Xlims, ylims=Ylims, mesh=Mesh, unitlength=Int(Unit[1]))
    save_png(path_png_strain, M, colorfunc, xlims=Xlims, ylims=Ylims, unitlength=Int(aa*Unit[1]))
    ColorBar(max=MaximumStrain, filename=path_png_colorbar, width=aa*Colorbarsize*Width)

    img_nurbs = load(path_png_nurbs)
    img_strain = load(path_png_strain)
    img_colorbar = load(path_png_colorbar)

    img_nurbs = convert(Array{RGBA{Float64},2}, img_nurbs)
    img_strain = convert(Array{RGBA{Float64},2}, img_strain)
    img_colorbar = convert(Array{RGBA{Float64},2}, img_colorbar)

    size_nurbs = size(img_nurbs)
    size_strain = size(img_strain)
    size_colorbar = size(img_colorbar)

    img_nurbs_white_background = img_nurbs ./ RGB(1, 1, 1)
    img_strain_white_background = img_strain ./ RGB(1, 1, 1)
    Δ = collect(size_strain) - collect(size_colorbar)

    img_offset_colorbar = OffsetArray(img_colorbar, Δ...)
    img_strain_with_colorbar = copy(img_strain_white_background)
    img_strain_with_colorbar[axes(img_offset_colorbar)...] = img_offset_colorbar ./ img_strain_with_colorbar[axes(img_offset_colorbar)...]
    img_strain_with_colorbar = [RGB(mean(img_strain_with_colorbar[5i-4:5i, 5j-4:5j])) for i in 1:size_nurbs[1], j in 1:size_nurbs[2]]
    # img_strain_with_colorbar = imresize(img_strain_with_colorbar, (800,800)) # could be coded like this, but the previous one is better for anti-alias
    img_append = hcat(img_nurbs_white_background, img_strain_with_colorbar)

    save(path_png_append, img_append)
end


"""
    computed_shapes()

Get names of computed shapes in the output directory (default: `~/ElasticSurfaceEmbedding-Result`)
"""
function computed_shapes()
    shapes = String[]
    for name in readdir(OUT_DIR)
        dir = joinpath(OUT_DIR, name)
        !isdir(dir) && continue
        !isfile(joinpath(dir, name*".json")) && continue
        push!(shapes, name)
    end
    return shapes
end


# TODO
# export ReDraw
# function ReDraw()
#     return
# end
