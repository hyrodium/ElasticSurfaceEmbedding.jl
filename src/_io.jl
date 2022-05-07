# Default output directory
OUT_DIR = joinpath(homedir(),"ElasticSurfaceEmbedding-Result")

mutable struct Step{T<:BSplineManifold{2}}
    manifold::T
    comment::String
    pinned::Bool
    function Step(manifold::BSplineManifold{2},comment,pinned=false)
        new{typeof(manifold)}(manifold,comment,pinned)
    end
end

struct AllSteps
    steps::Vector{Tuple{Step,Int}}
    function AllSteps()
        new(Vector{Tuple{Step,Int}}())
    end
end

function addstep!(allsteps::AllSteps, step::Step, parent::Int)
    push!(allsteps.steps, (step, parent))
    return allsteps
end

function parent_id(allsteps, id)
    allsteps.steps[id][2]
end

function children_ids(allsteps, id)
    findall(step->step[2]==id, allsteps.steps)
end

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
end

function nodeseries(allsteps, i)
    series = [i]
    while i ≠ 0
        i = parent_id(allsteps, i)
        pushfirst!(series, i)
    end
    return series
end

function _tree_as_string(allsteps::AllSteps)
    n = length(allsteps.steps)
    serieses = [nodeseries(allsteps, i) for i in 1:n]
    sort!(serieses)
    lowstrings = String[]
    for i in 1:n
        l = length(serieses[i])
        key = serieses[i][end]
        step = allsteps.steps[key][1]
        pinned = step.pinned
        comment = "📌 "^pinned * step.comment
        if l == 2
            lowstring = "$(key): " * comment
            push!(lowstrings, lowstring)
        elseif l ≥ 3
            lowstring = "  "^(l - 3) * "└─$(key): " * comment
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

function Base.show(io::IO, allsteps::AllSteps)
    print(_tree_as_string(allsteps))
end

function _realparent(allsteps, index::Int)
    if index==0
        return length(allsteps.steps)
    else
        return index
    end
end

function loadM(allsteps; index=0)
    if index == 0
        index = length(allsteps.steps)
    end
    M = allsteps.steps[index][1].manifold
    return M
end

function _export(M::BSplineManifold{2}, index::Int; comment="", maximumstrain=MAXIMUMSTRAIN)
    mkpath(DIR)

    # Compute maximum strain for graphics export
    if iszero(maximumstrain)
        MS = _compute_minmax_strain(M)
        maximumstrain = max(-MS[1], MS[2])
    end

    # Save graphics
    ExportFiles(M, maximumstrain, index)
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

    normalized_strain(u¹, u²) = E⁽⁰⁾₁₁(M, u¹, u²) / MaximumStrain # bounded in -1 to 1

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
    _colorbar(max=MaximumStrain, filename=path_png_colorbar, width=aa*Colorbarsize*Width)

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

# TODO
# export ReDraw
# function ReDraw()
#     return
# end
