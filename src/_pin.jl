"""
    add_pin(; parent::Int = 0, tag::String = "")

Add a pin 📌 for the given index
"""
function add_pin(; parent::Int = 0, tag::String = "")
    if tag == ""
        tag = Dates.format(now(), "yyyy-mm-dd_H-M-S")
    end
    _check_filename(tag)
    if _tag_exists(tag)
        error("The tag $(tag) is already exists.")
    end
    dict = _loadresultdict()
    parent = _realparent(parent)
    M = loadM(index = parent)

    index = latest_index(dict) + 1
    dict["result"][string(index)] = Dict{String,Any}("parent" => string(parent))
    dict["result"][string(index)]["BSplineManifold"] = toJSON(M)

    comment = "📌 - tag: " * tag
    dict["result"][string(index)]["comment"] = comment

    # Save as json
    mkpath(DIR)
    write(joinpath(DIR, NAME*".json"), JSON.json(dict, 4))

    # Send messages
    allsteps = _tmp_allsteps_from_result(dict["result"])
    message = _tree_as_string(allsteps)
    println(message)
    _send_file_to_slack("", comment="```\n" * message * "```")
end

function _tag_exists(tag)
    pinned_states = _find_all_pinned_states()
    for index in pinned_states
        if _get_tag(index) == tag
            return true
        end
    end
    return false
end

function _get_tag(index)::String
    dict = _loadresultdict()

    comment = dict["result"][repr(index)]["comment"]
    if startswith(comment, "📌 ")
        return replace(comment, "📌 - tag: " => "")
    else
        error("The index $(index) is not pinned.")
    end
end

"""
    remove_pin(index::Integer)

Remeve a pin 💨 with the given index
"""
function remove_pin(index::Integer)
    if !(index in _find_all_pinned_states())
        @warn "There's no pin with the index $index"
        return
    end
    dict = _loadresultdict()
    comment = dict["result"][repr(index)]["comment"]
    comment = replace(comment, "📌" => "💨")
    dict["result"][repr(index)]["comment"] = comment

    # Save as json
    mkpath(DIR)
    write(joinpath(DIR, NAME*".json"), JSON.json(dict, 4))

    # Send messages
    allsteps = _tmp_allsteps_from_result(dict["result"])
    message = _tree_as_string(allsteps)
    println(message)
    _send_file_to_slack("", comment="```\n" * message * "```")
end

"""
    remove_pin(; tag::AbstractString)

Remeve a pin 💨 with the given tag
"""
function remove_pin(; tag::AbstractString)
    pinned_states = _find_all_pinned_states()
    for index in pinned_states
        if _get_tag(index) == tag
            remove_pin(index)
            return
        end
    end
    @warn "There's no pin with the tag $tag"
end

function _find_all_pinned_states()
    dict = _loadresultdict()
    pinned_states = Int[]
    for key in keys(dict["result"])
        comment = dict["result"][key]["comment"]
        if startswith(comment, "📌 ")
            push!(pinned_states, parse(Int, key))
        end
    end
    return pinned_states
end

"""
    export_all_pinned_states(; unitlength = (10, "mm"), cutout = (0.1, 5), mesh::Int = 60)

Export all pinned states for final output
"""
function export_all_pinned_states(; unitlength::Tuple{<:Real,<:AbstractString}, cutout=(0.1, 5), mesh::Int=60)
    dir_pinned = joinpath(DIR, "pinned")
    # Delete current pinned directory
    rm(dir_pinned, recursive=true, force=true)
    # Make path to pinned directory
    mkpath(dir_pinned)

    pinned_states = _find_all_pinned_states()

    for index in pinned_states
        M = loadM(index=index)
        filename = joinpath(DIR, "pinned", "$(_get_tag(index)).svg")
        save_svg(filename, M, xlims=XLIMS, ylims=YLIMS, mesh=MESH, unitlength=unitlength[1], points=false)

        P = bsplinespaces(M)
        p₁, p₂ = degree.(P)
        k₁, k₂ = knotvector.(P)
        D₁, D₂ = k₁[1+p₁]..k₁[end-p₁], k₂[1+p₂]..k₂[end-p₂]

        𝒆⁽⁰⁾₁(u¹,u²) = normalize(𝒑₁₍ₜ₎(M,u¹,u²))
        𝒆⁽⁰⁾₂(u¹,u²) = [0.0 -1.0; 1.0 0.0] * 𝒆⁽⁰⁾₁(u¹,u²)
        𝒑a(i, t) = 𝒑₍ₜ₎(M, t, leftendpoint(D₂)) + 𝒆⁽⁰⁾₂(t, leftendpoint(D₂)) * i * cutout[1] / unitlength[1]
        𝒑b(i, t) = 𝒑₍ₜ₎(M, t, rightendpoint(D₂)) - 𝒆⁽⁰⁾₂(t, rightendpoint(D₂)) * i * cutout[1] / unitlength[1]
        _svgcurve(
            [[t -> 𝒑a(i, t) for i in 0:cutout[2]]..., [t -> 𝒑b(i, t) for i in 0:cutout[2]]...],
            D₁,
            filename = joinpath(DIR, "pinned", "$(_get_tag(index))-cutout.svg"),
            up = UP,
            down = DOWN,
            right = RIGHT,
            left = LEFT,
            thickness = 0.1,
            mesh = mesh,
            unitlength = unitlength[1]
        )
    end

    for name in readdir(dir_pinned)
        file = joinpath(dir_pinned, name)
        _changeunit(file, "pt"=>unitlength[2])
    end
end
