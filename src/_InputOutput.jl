"""
A macro to define parametrized shape of surface
"""
macro parametric_mapping(ex)
    expr = toJSON(ex)
    if startswith(expr, "function 𝒑₍₀₎(u)\n") || startswith(expr, "𝒑₍₀₎(u) =")
        global EXPR = expr
    else
        error("symbol of parametric mapping must be 𝒑₍₀₎(u)")
    end
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
    global UP = canvas[2]/2
    global DOWN = -canvas[2]/2
    global RIGHT = canvas[1]/2
    global LEFT = -canvas[1]/2
    global MESH = mesh
    global UNIT = (unit, "pt")
    global MAXIMUMSTRAIN = maximumstrain
    global COLORBARSIZE = colorbarsize
    if isTheShapeComputed()
        dict = LoadResultDict()
        println(TreeString(dict["result"]))
        global EXPR = dict["expr"]
    elseif !(@isdefined EXPR)
        error("use @parametric_mapping or input name of computed shape")
    end
    eval(:($(Meta.parse(EXPR))))
    return
end

function toJSON(ex::Expr)::String
    return replace(string(ex), r"#=.*=#\n" => s"")
end
function toJSON(k::Knots)
    return k.vector
end
function toJSON(P::FastBSplineSpace)
    return Dict("degree" => degree(P), "knots" => toJSON(knots(P)))
end
function toJSON(P::Array{P,1} where {P<:FastBSplineSpace})
    return convert(Array{Any,1}, toJSON.(P))
end
function toJSON(a::Array{Float64})
    return convert(Array{Any,1}, a[:])
end
function toJSON(M::AbstractBSplineManifold)
    return Dict("bsplinespaces" => toJSON(collect(bsplinespaces(M))), "controlpoints" => toJSON(controlpoints(M)))
end

function JSONtoBSplineSpace(jP::Dict)
    return FastBSplineSpace(jP["degree"], Knots(Vector{Float64}(jP["knots"])))
end
function JSONtoBSplineSpaces(jPs::Array)
    return [JSONtoBSplineSpace(jP) for jP in jPs]
end
function JSONtoControlPoints(ja, dims)
    n = prod(dims)
    d = length(ja) ÷ n
    return reshape(convert(Array{Float64}, ja), dims..., d)
end
function JSONtoBSplineSurface(dict::Dict)
    P = JSONtoBSplineSpaces(dict["bsplinespaces"])
    dims = dim.(P)
    𝒂 = JSONtoControlPoints(dict["controlpoints"], dims)
    return BSplineSurface(P, 𝒂)
end

function NodeSeries(tree::Dict, node)
    Nodes = [node]
    while Nodes[end] ≠ "0"
        push!(Nodes, tree[Nodes[end]]["parent"])
    end
    return Nodes
end

function TreeString(tree::Dict)
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

function latest_index(; dict::Union{Dict,Nothing} = nothing)
    if isnothing(dict)
        dict = LoadResultDict()
    end
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
    if EXPR ≠ dict["expr"]
        println(EXPR)
        println(dict["expr"])
        # error("The definition of 𝒑₍₀₎(u) has been changed")
    end
    index = _realparent(index)
    M = JSONtoBSplineSurface(dict["result"][string(index)]["BSplineManifold"])
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

function _export(M::AbstractBSplineManifold, parent::Int; comment = "", maximumstrain = MAXIMUMSTRAIN)
    if isTheShapeComputed()
        dict = LoadResultDict()
        index = latest_index(dict=dict) + 1
    else
        meta = Dict{String,Any}("version" => ESE_VERSION, "generated by" => "ElasticSurfaceEmbedding.jl")
        dict = Dict{String,Any}("expr" => EXPR, "meta" => meta)
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
    message = TreeString(dict["result"])
    println(message)
    _send_file_to_slack(path_png_append, comment = "```\n" * message * "```")
end

function ExportFiles(
    M::AbstractBSplineManifold,
    MaximumStrain::Real,
    index;
    Name::String = NAME,
    Dir = DIR,
    Up = UP,
    Down = DOWN,
    Right = RIGHT,
    Left = LEFT,
    Mesh = MESH,
    Unit = UNIT,
    Colorbarsize = COLORBARSIZE,
)
    mkpath(joinpath(DIR, "nurbs"))
    mkpath(joinpath(DIR, "strain"))
    mkpath(joinpath(DIR, "colorbar"))
    mkpath(joinpath(DIR, "append"))

    Width = (Right - Left) * Unit[1]
    Height = (Up - Down) * Unit[1]

    normalized_strain(u) = E⁽⁰⁾₁₁_cont(M, u) / MaximumStrain # bounded in -1 to 1

    aa = 5 # magnification parameter for antialias

    path_svg_nurbs = joinpath(Dir, "nurbs", "$(Name)-$(index)_Bspline.svg")
    path_png_nurbs = joinpath(Dir, "nurbs", "$(Name)-$(index)_Bspline.png")
    path_png_strain = joinpath(Dir, "strain", "$(Name)-$(index)_strain.png")
    path_png_colorbar = joinpath(Dir, "colorbar", "$(Name)-$(index)_colorbar.png")
    path_png_append = joinpath(Dir, "append", "$(Name)-$(index)_append.png")

    colorfunc(u) = normalized_strain(u) * RGB(0.5, -0.5, -0.5) + RGB(0.5, 0.5, 0.5) # red to cyan

    save_svg(path_svg_nurbs, M, up = Up, down = Down, right = Right, left = Left, mesh = Mesh, unitlength = Int(Unit[1]))
    save_png(path_png_nurbs, M, up = Up, down = Down, right = Right, left = Left, mesh = Mesh, unitlength = Int(Unit[1]))
    save_png(path_png_strain, M, colorfunc, up = Up, down = Down, right = Right, left = Left, unitlength = Int(aa * Unit[1]))
    ColorBar(max = MaximumStrain, filename = path_png_colorbar, width = aa * Colorbarsize * Width)

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
