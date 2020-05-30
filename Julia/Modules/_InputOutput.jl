using ParametricDraw

export @ParametricMapping
macro ParametricMapping(ex)
    expr=toJSON(ex)
    if startswith(expr,"function 𝒑₍₀₎(u)\n") || startswith(expr,"𝒑₍₀₎(u) =")
        global EXPR=expr
    else
        error("symbol of parametric mapping must be 𝒑₍₀₎(u)")
    end
end

function CheckAsFileName(name::String)
    invalidcharacters=[' ', '&', '\\', '/', '.', '<', '>', '|', ':', ';', '*', '?', '=', '%', '$', '"', '~']
    if length(invalidcharacters ∩ name) ≠ 0
        error("The name[$(name)] must not consists of following chars: ", invalidcharacters)
    end
end

export Settings
function Settings(name::String; up::Real=5, down::Real=-5, right::Real=5, left::Real=-5, mesh::Tuple{Int,Int}=(10,1), unit::Real=100, slack::Bool=true, maximumstrain::Real=0.0, colorbarsize::Float64=0.2)
    CheckAsFileName(name)
    global NAME=name
    global DIR=OUT_DIR*NAME
    global UP=up
    global DOWN=down
    global RIGHT=right
    global LEFT=left
    global MESH=mesh
    global UNIT=(unit,"pt")
    global SLACK=slack
    global MAXIMUMSTRAIN=maximumstrain
    global COLORBARSIZE=colorbarsize
    if isTheShapeComputed()
        dict=LoadResultDict()
        println(TreeString(dict["Result"]))
        global EXPR=dict["Expr"]
    elseif !(@isdefined EXPR)
        error("use @ParametricMapping or input name of computed shape")
    end
    eval(:(@everywhere $(Meta.parse(EXPR))))
    return nothing
end

function toJSON(ex::Expr) ::String
    return replace(string(ex), r"#=.*=#\n" => s"")
end
function toJSON(k::Knots)
    return k.vector
end
function toJSON(P::FastBSplineSpace)
    return Dict("degree"=>degree(P), "knots"=>toJSON(knots(P)))
end
function toJSON(P::Array{P,1} where P <: FastBSplineSpace)
    return convert(Array{Any,1}, toJSON.(P))
end
function toJSON(𝒂::Array{Float64})
    return convert(Array{Any,1}, 𝒂[:])
end
function toJSON(M::FastBSplineManifold)
    return Dict("bsplinespaces"=>toJSON(M.bsplinespaces), "controlpoints"=>toJSON(M.controlpoints))
end

function JSONtoBSplineSpace(jP::Dict)
    return FastBSplineSpace(jP["degree"],Knots(jP["knots"]))
end
function JSONtoBSplineSpaces(jPs::Array)
    return [JSONtoBSplineSpace(jP) for jP ∈ jPs]
end
function JSONtoControlPoints(j𝒂, dims)
    n = prod(dims)
    d = length(j𝒂) ÷ n
    return reshape(convert(Array{Float64},j𝒂),dims...,d)
end
function JSONtoFastBSplineManifold(dict::Dict)
    P = JSONtoBSplineSpaces(dict["bsplinespaces"])
    dims = dim.(P)
    𝒂 = JSONtoControlPoints(dict["controlpoints"], dims)
    return FastBSplineManifold(P,𝒂)
end

function NodeSeries(tree::Dict,node)
    Nodes=[node]
    while Nodes[end] ≠ "0"
        push!(Nodes, tree[Nodes[end]]["parent"])
    end
    return Nodes
end

function TreeString(tree::Dict)
    serieses=Array{Int,1}[]
    for key in keys(tree)
        push!(serieses,(s->parse(Int,s)).(reverse(NodeSeries(tree,key))))
    end
    sort!(serieses)
    lowstrings=String[]
    n = length(serieses)
    for i in 1:n
        l=length(serieses[i])
        key=string(serieses[i][end])
        comment=tree[key]["comment"]
        if l == 2
            lowstring=key*": "*comment
            push!(lowstrings,lowstring)
        elseif l ≥ 3
            lowstring="  "^(l-3)*"└─"*key*": "*comment
            push!(lowstrings,lowstring)
            for j in 1:(i-1)
                chars=collect(lowstrings[end-j])
                if chars[2(l-3)+1]==' '
                    lowstrings[end-j]=join(chars[1:2(l-3)])*"│"*join(chars[2(l-3)+2:end])
                elseif chars[2(l-3)+1]=='└'
                    lowstrings[end-j]=join(chars[1:2(l-3)])*"├"*join(chars[2(l-3)+2:end])
                    break
                else
                    break
                end
            end
        end
    end
    outsting=""
    for s in lowstrings
        outsting=outsting*s*"\n"
    end
    return outsting
end


function isTheShapeComputed()
    return isfile(DIR*"/"*NAME*".json")
end

function NewestIndex(; dict::Union{Dict,Nothing}=nothing)
    if dict isa Nothing
        dict=LoadResultDict()
    end
    result_nums = [parse(Int, i) for i in keys(dict["Result"])]
    return maximum(result_nums)
end

function Parent(index::Union{Int,Nothing})
    dict=LoadResultDict()
    if index==0
        return NewestIndex()
    elseif index==nothing
        return NewestIndex()
    else
        return index
    end
end

function loadM(; index=0, dict::Union{Dict,Nothing}=nothing)
    if !isTheShapeComputed()
        error("Result file doesn't exists")
    end
    dict=LoadResultDict()
    if EXPR ≠ dict["Expr"]
        println(EXPR)
        println(dict["Expr"])
        # error("The definition of 𝒑₍₀₎(u) has been changed")
    end
    index=Parent(index)
    M=JSONtoFastBSplineManifold(dict["Result"][string(index)]["FastBSplineManifold"])
    return M
end

function LoadResultDict() :: Dict
    f = open(DIR*"/"*NAME*".json","r")
    dict = JSON.parse(f)
    close(f)
    v = dict["Version"]
    JSON_VERSION = VersionNumber(v["major"],v["minor"],v["patch"])

    if JSON_VERSION.major == ESE_VERSION.major
        return dict
    else
        error("The version of computed shape is not supported.")
    end
end

function SaveResultDict(dict::Dict)
    mkpath(DIR)
    open(DIR*"/"*NAME*".json","w") do f
        JSON.print(f, dict, 4)
    end
    PrintResultDict(dict)
    return nothing
end

function PrintResultDict(dict::Dict; slack=SLACK)
    treestring = TreeString(dict["Result"])
    println(treestring)
    if slack
        SlackString("```\n"*treestring*"```")
    end
end

function Export(M::FastBSplineManifold, parent::Int; comment="", maximumstrain=MAXIMUMSTRAIN)
    if isTheShapeComputed()
        dict = LoadResultDict()
        index = NewestIndex(dict=dict) +1
    else
        dict=Dict{String,Any}("Expr" => EXPR, "Version" => ESE_VERSION)
        dict["Result"] = Dict{String,Any}()
        index = 1
    end

    dict["Result"][string(index)] = Dict{String,Any}("parent" => string(parent))
    dict["Result"][string(index)]["FastBSplineManifold"] = toJSON(M)
    dict["Result"][string(index)]["comment"] = comment

    SaveResultDict(dict)

    if maximumstrain == 0.0
        MS = ComputeMaximumStrain(index=index)
        MaximumStrain = max(-MS[1],MS[2])
    else
        MaximumStrain = maximumstrain
    end

    if distributed
        # @spawnat 1 ExportFiles(M, MaximumStrain, index)
        ExportFiles(M, MaximumStrain, index)
    else
        ExportFiles(M, MaximumStrain, index)
    end

    return nothing
end

function ExportFiles(M::FastBSplineManifold, MaximumStrain::Real, index; Name::String=NAME, Dir=DIR, Up=UP, Down=DOWN, Right=RIGHT, Left=LEFT, Mesh=MESH, Unit=UNIT, Slack::Bool=SLACK, Colorbarsize=COLORBARSIZE)
    mkpath(DIR*"/nurbs")
    mkpath(DIR*"/strain")
    mkpath(DIR*"/colorbar")
    mkpath(DIR*"/append")

    dict=LoadResultDict()
    𝒂 = M.controlpoints
    P₁,P₂ = P = M.bsplinespaces
    p₁,p₂ = p = degree.(P)
    k₁,k₂ = k = knots.(P)
    D₁,D₂ = D = k₁[1+p₁]..k₁[end-p₁],k₂[1+p₂]..k₂[end-p₂]

    Width = (Right-Left)*Unit[1]
    Height = (Up-Down)*Unit[1]

    rgb(u) = E⁽⁰⁾₁₁(M,u)*[1,-1,-1]/(2*MaximumStrain) .+0.5 # Red to Cyan

    aa = 5 # magnification parameter for antialias

    path_svg_nurbs=Dir*"/nurbs/"*Name*"-"*string(index)*"_Bspline.svg"
    path_png_nurbs=Dir*"/nurbs/"*Name*"-"*string(index)*"_Bspline.png"
    path_png_strain=Dir*"/strain/"*Name*"-"*string(index)*"_strain.png"
    path_png_colorbar=Dir*"/colorbar/"*Name*"-"*string(index)*"_colorbar.png"
    path_png_append=Dir*"/append/"*Name*"-"*string(index)*"_append.png"

    DrawBSpline(M, filename=path_svg_nurbs, up=Up, down=Down, right=Right, left=Left, mesh=Mesh, unitlength=Int(Unit[1]))
    DrawBSpline(M, filename=path_png_nurbs, up=Up, down=Down, right=Right, left=Left, mesh=Mesh, unitlength=Int(Unit[1]))
    ParametricColor(u->𝒑₍ₜ₎(M,u), D, rgb=rgb,  filename=path_png_strain, up=Up, down=Down, right=Right, left=Left, mesh=tuple(10*[Mesh...]...), unit=aa*Unit[1])
    ColorBar(max=MaximumStrain, filename=path_png_colorbar, width=aa*Colorbarsize*Width)

    img_nurbs = load(path_png_nurbs)
    img_strain = load(path_png_strain)
    img_colorbar = load(path_png_colorbar)

    img_nurbs = convert(Array{RGB{Float64},2},img_nurbs)
    img_strain = convert(Array{RGBA{Float64},2},img_strain)
    img_colorbar = convert(Array{RGBA{Float64},2},img_colorbar)

    size_nurbs = size(img_nurbs)
    size_strain = size(img_strain)
    size_colorbar = size(img_colorbar)

    img_strain_white_background = img_strain ./ RGB(1,1,1)
    Δ = collect(size_strain) - collect(size_colorbar)

    img_offset_colorbar = OffsetArray(img_colorbar, Δ...)
    img_strain_with_colorbar = copy(img_strain_white_background)
    img_strain_with_colorbar[axes(img_offset_colorbar)...] = img_offset_colorbar ./ img_strain_with_colorbar[axes(img_offset_colorbar)...]
    img_strain_with_colorbar = [RGB(mean(img_strain_with_colorbar[5i-4:5i,5j-4:5j])) for i in 1:800, j in 1:800]
    # img_strain_with_colorbar = imresize(img_strain_with_colorbar, (800,800)) # could be coded like this, but previous one is better for anti-alias
    img_append = hcat(img_nurbs, img_strain_with_colorbar)

    save(path_png_append, img_append)

    if Slack
        SlackFile(path_png_append, comment=Name*"-"*string(index))
    end
end


export ComputedShapes
"""
    ComputedShapes()

show and return the computed shapes.
"""
function ComputedShapes()
    shapes = Base.Filesystem.readdir(ElasticSurfaceEmbedding.OUT_DIR)
    println(shapes)
    return shapes
end


# TODO
# export ReDraw
# function ReDraw()
#     return nothing
# end
