"""
    pin(; parent::Int = 0)

Add a pin 📌 for the given index
"""
function pin(allsteps; index::Int=0)
    index = _realparent(index)
    allsteps.steps[index][1].pinned = true
    return allsteps
end

"""
    unpin(index::Integer)

Remeve the pin 📌 with the given index
"""
function unpin(allsteps; index::Int=0)
    index = _realparent(index)
    allsteps.steps[index][1].pinned = false
    return allsteps
end

function _find_all_pinned_states(allsteps)
    return findall(step -> step[1].pinned, allsteps.steps)
end

"""
    export_pinned_states(; unitlength = (10, "mm"), cutout = (0.1, 5), mesh::Int = 60)

Export all pinned states for final output
"""
function export_pinned_states(allsteps; unitlength::Tuple{<:Real,<:AbstractString}, cutout=(0.1, 5), mesh::Int=60)
    dir_pinned = joinpath(DIR, "pinned")
    # Delete current pinned directory
    rm(dir_pinned, recursive=true, force=true)
    # Make path to pinned directory
    mkpath(dir_pinned)

    pinned_states = _find_all_pinned_states(allsteps)

    for index in pinned_states
        M = loadM(index=index)
        filename = joinpath(DIR, "pinned", "$(index).svg")
        save_svg(filename, M, xlims=XLIMS, ylims=YLIMS, mesh=MESH, unitlength=unitlength[1], points=false)

        # P = bsplinespaces(M)
        # p₁, p₂ = degree.(P)
        # k₁, k₂ = knotvector.(P)
        # D₁, D₂ = k₁[1+p₁]..k₁[end-p₁], k₂[1+p₂]..k₂[end-p₂]

        # 𝒆⁽⁰⁾₁(u¹,u²) = normalize(𝒑₁₍ₜ₎(M,u¹,u²))
        # 𝒆⁽⁰⁾₂(u¹,u²) = [0.0 -1.0; 1.0 0.0] * 𝒆⁽⁰⁾₁(u¹,u²)
        # 𝒑a(i, t) = 𝒑₍ₜ₎(M, t, leftendpoint(D₂)) + 𝒆⁽⁰⁾₂(t, leftendpoint(D₂)) * i * cutout[1] / unitlength[1]
        # 𝒑b(i, t) = 𝒑₍ₜ₎(M, t, rightendpoint(D₂)) - 𝒆⁽⁰⁾₂(t, rightendpoint(D₂)) * i * cutout[1] / unitlength[1]
        # _svgcurve(
        #     [[t -> 𝒑a(i, t) for i in 0:cutout[2]]..., [t -> 𝒑b(i, t) for i in 0:cutout[2]]...],
        #     D₁,
        #     filename = joinpath(DIR, "pinned", "$(_get_tag(index))-cutout.svg"),
        #     up = UP,
        #     down = DOWN,
        #     right = RIGHT,
        #     left = LEFT,
        #     thickness = 0.1,
        #     mesh = mesh,
        #     unitlength = unitlength[1]
        # )
    end

    for name in readdir(dir_pinned)
        file = joinpath(dir_pinned, name)
        _changeunit(file, "pt"=>unitlength[2])
    end
end
