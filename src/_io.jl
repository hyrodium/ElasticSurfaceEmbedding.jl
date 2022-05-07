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

function _check_filename(name::String)
    invalidcharacters = [' ', '&', '\\', '/', '.', '<', '>', '|', ':', ';', '*', '?', '=', '%', '$', '"', '~']
    if length(invalidcharacters ∩ name) ≠ 0
        error("The name[$(name)] must not consists of following chars: ", invalidcharacters)
    end
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
    print(io, _tree_as_string(allsteps))
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

function export_all_steps(
        dir,
        allsteps::AllSteps;
        maximumstrain = 0,
        xlims = (-5,5),
        ylims = (-5,5),
        mesh = (10,1),
        unit = (100,"mm"),
        colorbarsize = 0.3,
    )
    for (step, i) in allsteps.steps
        M = step.manifold
        _export_one_step(dir, M, i, maximumstrain=maximumstrain, xlims=xlims, ylims=ylims, mesh=mesh, unit=unit, colorbarsize=colorbarsize)
    end
end

function _export_one_step(
        dir,
        M::BSplineManifold{2},
        index::Integer;
        maximumstrain = 0,
        xlims = (-5,5),
        ylims = (-5,5),
        mesh = (10,1),
        unit = (100,"mm"),
        colorbarsize = 0.3,
    )
    if maximumstrain ≤ 0
        MS = _compute_minmax_strain(M)
        maximumstrain = max(-MS[1], MS[2])
    end

    width = (xlims[2] - xlims[1]) * unit[1]

    normalized_strain(u¹, u²) = E⁽⁰⁾₁₁(M, u¹, u²) / maximumstrain # bounded in -1 to 1

    aa = 5 # magnification parameter for antialias

    path_svg_bspline = joinpath(dir, "bspline-$(index).svg")
    path_png_bspline = joinpath(dir, "bspline-$(index).png")
    path_png_strain = joinpath(dir, "strain-$(index).png")
    path_png_colorbar = joinpath(dir, "colorbar-$(index).png")
    path_png_append = joinpath(dir, "append-$(index).png")

    colorfunc(u¹,u²) = normalized_strain(u¹,u²) * RGB(0.5, -0.5, -0.5) + RGB(0.5, 0.5, 0.5) # red to cyan

    save_svg(path_svg_bspline, M, xlims=xlims, ylims=ylims, mesh=mesh, unitlength=Int(unit[1]))
    save_png(path_png_bspline, M, xlims=xlims, ylims=ylims, mesh=mesh, unitlength=Int(unit[1]))
    save_png(path_png_strain, M, colorfunc, xlims=xlims, ylims=ylims, unitlength=Int(aa*unit[1]))
    _colorbar(max=maximumstrain, filename=path_png_colorbar, width=aa*colorbarsize*width)

    img_bspline = load(path_png_bspline)
    img_strain = load(path_png_strain)
    img_colorbar = load(path_png_colorbar)

    img_bspline = convert(Array{RGBA{Float64},2}, img_bspline)
    img_strain = convert(Array{RGBA{Float64},2}, img_strain)
    img_colorbar = convert(Array{RGBA{Float64},2}, img_colorbar)

    size_bspline = size(img_bspline)
    size_strain = size(img_strain)
    size_colorbar = size(img_colorbar)

    img_bspline_white_background = img_bspline ./ RGB(1, 1, 1)
    img_strain_white_background = img_strain ./ RGB(1, 1, 1)
    Δ = collect(size_strain) - collect(size_colorbar)

    img_offset_colorbar = OffsetArray(img_colorbar, Δ...)
    img_strain_with_colorbar = copy(img_strain_white_background)
    img_strain_with_colorbar[axes(img_offset_colorbar)...] = img_offset_colorbar ./ img_strain_with_colorbar[axes(img_offset_colorbar)...]
    img_strain_with_colorbar = [RGB(mean(img_strain_with_colorbar[5i-4:5i, 5j-4:5j])) for i in 1:size_bspline[1], j in 1:size_bspline[2]]
    # img_strain_with_colorbar = imresize(img_strain_with_colorbar, (800,800)) # could be coded like this, but the previous one is better for anti-alias
    img_append = hcat(img_bspline_white_background, img_strain_with_colorbar)

    save(path_png_append, img_append)
end
