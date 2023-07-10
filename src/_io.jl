mutable struct Step{T<:BSplineManifold{2}}
    manifold::T
    comment::String
    info::Dict
    function Step(manifold::BSplineManifold{2}, comment, info)
        new{typeof(manifold)}(manifold, comment, info)
    end
end

struct StepTree
    steps::Vector{Step}
    parents::Vector{Int}
    pinned::Vector{Bool}
    function StepTree()
        new(Vector{Step}(), Vector{Int}(), Vector{Bool}())
    end
end

function addstep!(steptree::StepTree, step::Step, parent::Int)
    push!(steptree.steps, step)
    push!(steptree.parents, parent)
    push!(steptree.pinned, false)
    return steptree
end

function parent_id(steptree, id)
    steptree.parents[id]
end

function nodeseries(steptree, i)
    series = [i]
    while i ≠ 0
        i = parent_id(steptree, i)
        pushfirst!(series, i)
    end
    return series
end

function _tree_as_string(steptree::StepTree)
    n = length(steptree.steps)
    serieses = [nodeseries(steptree, i) for i in 1:n]
    sort!(serieses)
    lowstrings = String[]
    for i in 1:n
        l = length(serieses[i])
        key = serieses[i][end]
        step = steptree.steps[key]
        pinned = steptree.pinned[key]
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

function Base.show(io::IO, steptree::StepTree)
    print(io, _tree_as_string(steptree))
end

function _validindex(steptree, index::Int)
    if index == 0
        return length(steptree.steps)
    else
        return index
    end
end

function loadM(steptree; index = 0)
    if index == 0
        index = length(steptree.steps)
    end
    M = steptree.steps[index].manifold
    return M
end

function export_all_steps(
    dir,
    steptree::StepTree;
    maximumstrain = 0,
    xlims = (-2, 2),
    ylims = (-2, 2),
    mesh = (10, 1),
    unitlength::Tuple{<:Real,<:AbstractString} = (50, "mm"),
    colorbarsize = 0.3,
)
    mkpath(dir)
    for i in eachindex(steptree.steps)
        M = steptree.steps[i].manifold
        export_one_step(
            dir,
            M,
            i,
            maximumstrain = maximumstrain,
            xlims = xlims,
            ylims = ylims,
            mesh = mesh,
            unitlength = unitlength,
            colorbarsize = colorbarsize,
        )
    end
    export_pinned_steps(dir, steptree, xlims = xlims, ylims = ylims, mesh = mesh, unitlength = unitlength)
    write(joinpath(dir, "log.txt"), _tree_as_string(steptree))
end

function export_one_step(
    dir,
    M::BSplineManifold{2},
    index::Integer;
    maximumstrain = 0,
    xlims = (-2, 2),
    ylims = (-2, 2),
    mesh = (10, 1),
    unitlength::Tuple{<:Real,<:AbstractString} = (50, "mm"),
    colorbarsize = 0.3,
)
    if maximumstrain ≤ 0
        MS = _compute_minmax_strain(M)
        maximumstrain = max(-MS[1], MS[2])
    end
    aa = 5 # magnification parameter for antialias
    width = (xlims[2] - xlims[1]) * unitlength[1]
    normalized_strain(u¹, u²) = E⁽⁰⁾₁₁(M, u¹, u²) / maximumstrain # bounded in -1 to 1

    mkpath(joinpath(dir, "bspline"))
    mkpath(joinpath(dir, "strain"))
    mkpath(joinpath(dir, "colorbar"))
    mkpath(joinpath(dir, "combined"))

    path_svg_bspline = joinpath(dir, "bspline", "bspline-$(index).svg")
    path_png_bspline = joinpath(dir, "bspline", "bspline-$(index).png")
    path_png_strain = joinpath(dir, "strain", "strain-$(index).png")
    path_png_colorbar = joinpath(dir, "colorbar", "colorbar-$(index).png")
    path_png_combined = joinpath(dir, "combined", "combined-$(index).png")

    colorfunc(u¹, u²) = normalized_strain(u¹, u²) * RGB(0.5, -0.5, -0.5) + RGB(0.5, 0.5, 0.5) # red to cyan

    save_svg(path_svg_bspline, M, xlims = xlims, ylims = ylims, mesh = mesh, unitlength = Int(unitlength[1]))
    save_png(path_png_bspline, M, xlims = xlims, ylims = ylims, mesh = mesh, unitlength = Int(unitlength[1]))
    save_png(path_png_strain, M, colorfunc, xlims = xlims, ylims = ylims, unitlength = Int(aa * unitlength[1]))
    _colorbar(max = maximumstrain, filename = path_png_colorbar, width = aa * colorbarsize * width)
    _changeunit(path_svg_bspline, "pt" => unitlength[2])

    img_bspline = load(path_png_bspline)
    img_strain = load(path_png_strain)
    img_colorbar = load(path_png_colorbar)

    img_bspline = convert(Array{RGBA{Float64},2}, img_bspline)
    img_strain = convert(Array{RGBA{Float64},2}, img_strain)
    img_colorbar = convert(Array{RGBA{Float64},2}, img_colorbar)

    size_bspline = size(img_bspline)
    size_strain = size(img_strain)
    size_colorbar = size(img_colorbar)

    img_bspline_white_background =
        ColorBlendModes.blend.(RGB(1, 1, 1), img_bspline, op = ColorBlendModes.CompositeSourceOver)
    img_strain_white_background =
        ColorBlendModes.blend.(RGB(1, 1, 1), img_strain, op = ColorBlendModes.CompositeSourceOver)
    Δ = size_strain .- size_colorbar

    img_offset_colorbar = OffsetArray(img_colorbar, Δ...)
    img_strain_with_colorbar = copy(img_strain_white_background)
    img_strain_with_colorbar[axes(img_offset_colorbar)...] =
        ColorBlendModes.blend.(
            img_strain_with_colorbar[axes(img_offset_colorbar)...],
            img_offset_colorbar,
            op = ColorBlendModes.CompositeSourceOver,
        )
    img_strain_with_colorbar =
        [RGB(mean(img_strain_with_colorbar[5i-4:5i, 5j-4:5j])) for i in 1:size_bspline[1], j in 1:size_bspline[2]]
    # img_strain_with_colorbar = imresize(img_strain_with_colorbar, (800,800)) # could be coded like this, but the previous one is better for anti-alias
    img_combined = hcat(img_bspline_white_background, img_strain_with_colorbar)

    save(path_png_combined, img_combined)
end

"""
    export_pinned_steps(; unitlength = (10, "mm"), cutout = (0.1, 5), mesh::Int = 60)

Export all pinned steps for final output
"""
function export_pinned_steps(
    dir::AbstractString,
    steptree::StepTree;
    xlims = (-2, 2),
    ylims = (-2, 2),
    mesh = (10, 1),
    unitlength::Tuple{<:Real,<:AbstractString} = (50, "mm"),
    # cutout=(0.1, 5),
)
    dir_pinned = joinpath(dir, "pinned")
    # Delete current pinned directory
    rm(dir_pinned, recursive = true, force = true)
    # Make path to pinned directory
    mkpath(dir_pinned)

    pinned_steps = findall(steptree.pinned)
    paths_output = Vector{String}(undef, length(pinned_steps))

    for (i, index) in enumerate(pinned_steps)
        M = loadM(steptree, index = index)
        path_svg = joinpath(dir_pinned, "pinned-$(index).svg")
        save_svg(path_svg, M, xlims = xlims, ylims = ylims, mesh = mesh, unitlength = unitlength[1], points = false)
        _changeunit(path_svg, "pt" => unitlength[2])
        paths_output[i] = path_svg
    end
    return paths_output
end
