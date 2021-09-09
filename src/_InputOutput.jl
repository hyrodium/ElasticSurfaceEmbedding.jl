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
    settings(name::String; up::Real=5, down::Real=-5, right::Real=5, left::Real=-5, mesh::Tuple{Int, Int}=(10,  1), unit::Real=100, slack::Bool=true, maximumstrain::Real=0.0, colorbarsize::Float64=0.2)

Initial settings for the given shape
"""
function settings(
    name::String;
    up::Real = 5,
    down::Real = -5,
    right::Real = 5,
    left::Real = -5,
    mesh::Tuple{Int,Int} = (10, 1),
    unit::Real = 100,
    slack::Bool = true,
    maximumstrain::Real = 0.0,
    colorbarsize::Float64 = 0.2,
)
    _check_filename(name)
    global NAME = name
    global DIR = OUT_DIR * NAME
    global UP = up
    global DOWN = down
    global RIGHT = right
    global LEFT = left
    global MESH = mesh
    global UNIT = (unit, "pt")
    global SLACK = slack
    global MAXIMUMSTRAIN = maximumstrain
    global COLORBARSIZE = colorbarsize
    if isTheShapeComputed()
        dict = LoadResultDict()
        println(TreeString(dict["Result"]))
        global EXPR = dict["Expr"]
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
    return FastBSplineSpace(jP["degree"], Knots(jP["knots"]))
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
    a = JSONtoControlPoints(dict["controlpoints"], dims)
    return BSplineSurface(P, a)
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
    return isfile(DIR * "/" * NAME * ".json")
end

function NewestIndex(; dict::Union{Dict,Nothing} = nothing)
    if isnothing(dict)
        dict = LoadResultDict()
    end
    result_nums = [parse(Int, i) for i in keys(dict["Result"])]
    return maximum(result_nums)
end

function Parent(index::Union{Int,Nothing})
    dict = LoadResultDict()
    if index == 0
        return NewestIndex()
    elseif isnothing(index)
        return NewestIndex()
    else
        return index
    end
end

function loadM(; index = 0, dict::Union{Dict,Nothing} = nothing)
    if !isTheShapeComputed()
        error("Result file doesn't exists")
    end
    dict = LoadResultDict()
    if EXPR ≠ dict["Expr"]
        println(EXPR)
        println(dict["Expr"])
        # error("The definition of 𝒑₍₀₎(u) has been changed")
    end
    index = Parent(index)
    M = JSONtoBSplineSurface(dict["Result"][string(index)]["BSplineManifold"])
    return M
end

function LoadResultDict()::Dict
    f = open(DIR * "/" * NAME * ".json", "r")
    dict = JSON.parse(f)
    close(f)
    v = dict["Version"]
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
        index = NewestIndex(dict = dict) + 1
    else
        dict = Dict{String,Any}("Expr" => EXPR, "Version" => ESE_VERSION)
        dict["Result"] = Dict{String,Any}()
        index = 1
    end

    dict["Result"][string(index)] = Dict{String,Any}("parent" => string(parent))
    dict["Result"][string(index)]["BSplineManifold"] = toJSON(M)
    dict["Result"][string(index)]["comment"] = comment

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
    path_png_append = DIR * "/append/" * NAME * "-" * string(index) * "_append.png"
    message = TreeString(dict["Result"])
    println(message)
    if SLACK
        _send_file_to_slack(path_png_append, comment = "```\n" * message * "```")
    end
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
    mkpath(DIR * "/nurbs")
    mkpath(DIR * "/strain")
    mkpath(DIR * "/colorbar")
    mkpath(DIR * "/append")

    Width = (Right - Left) * Unit[1]
    Height = (Up - Down) * Unit[1]

    normalized_strain(u) = E⁽⁰⁾₁₁_cont(M, u) / MaximumStrain # bounded in -1 to 1

    aa = 5 # magnification parameter for antialias

    path_svg_nurbs = Dir * "/nurbs/" * Name * "-" * string(index) * "_Bspline.svg"
    path_png_nurbs = Dir * "/nurbs/" * Name * "-" * string(index) * "_Bspline.png"
    path_png_strain = Dir * "/strain/" * Name * "-" * string(index) * "_strain.png"
    path_png_colorbar = Dir * "/colorbar/" * Name * "-" * string(index) * "_colorbar.png"
    path_png_append = Dir * "/append/" * Name * "-" * string(index) * "_append.png"

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
Show and return the names of computed shapes in the directory (default: ~/ElasticSurfaceEmbedding-Result)
"""
function computed_shapes()
    shapes = readdir(ElasticSurfaceEmbedding.OUT_DIR)
    println(shapes)
    return shapes
end


# TODO
# export ReDraw
# function ReDraw()
#     return
# end
