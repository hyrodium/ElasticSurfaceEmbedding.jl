# BSpline
function N′(P₁::FastBSplineSpace, P₂::FastBSplineSpace, I₁, I₂, i, u)::Float64
    if i == 1
        return bsplinebasis′₊₀(I₁, P₁, u[1]) * bsplinebasis(I₂, P₂, u[2])
    else
        return bsplinebasis(I₁, P₁, u[1]) * bsplinebasis′₊₀(I₂, P₂, u[2])
    end
end

function N′_cont(P₁::FastBSplineSpace, P₂::FastBSplineSpace, I₁, I₂, i, u)::Float64
    if i == 1
        return bsplinebasis′(I₁, P₁, u[1]) * bsplinebasis(I₂, P₂, u[2])
    else
        return bsplinebasis(I₁, P₁, u[1]) * bsplinebasis′(I₂, P₂, u[2])
    end
end

"""
Affine transform of control points.
"""
function _affine(𝒂, A, b)
    # x'=Ax+b
    n₁, n₂, d = size(𝒂)
    return [(A*𝒂[I₁, I₂, :]+b)[i] for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:d]
end

function _rotate(𝒂)
    n₁, n₂, _ = size(𝒂)
    ind0 = [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2]
    ind1 = ind0 - [0, 1]
    v = 𝒂[ind1..., :] - 𝒂[ind0..., :]
    R = -[v[2] -v[1]; v[1] v[2]] / norm(v)
    return _affine(𝒂, R, [0.0, 0.0])
end

function _center(𝒂)
    x_min = minimum(𝒂[:,:,1])
    x_max = maximum(𝒂[:,:,1])
    y_min = minimum(𝒂[:,:,2])
    y_max = maximum(𝒂[:,:,2])
    x = (x_min+x_max)/2
    y = (y_min+y_max)/2
    return _affine(𝒂, I(2), -[x,y])
end

function _positioning(𝒂)
    return _center(_rotate(𝒂))
end

function _positioning(M::AbstractBSplineManifold)
    Ps = collect(bsplinespaces(M))
    𝒂 = controlpoints(M)
    if length(Ps) ≠ 2
        error("dimension does not match")
    end

    𝒂′ = _positioning(𝒂)
    return typeof(M)(Ps, 𝒂′)
end

"""
    spline_refinement(; p₊::Array{Int,1}=[0, 0], k₊::Array{Knots,1}=[Knots(), Knots()], parent::Int=0)

Compute a refinement of the B-spline manifold
"""
function spline_refinement(; p₊=(0,0), k₊=(Knots(),Knots()), parent::Int=0)
    parent = _realparent(parent)
    M = loadM(index=parent)

    P₁, P₂ = collect(bsplinespaces(M))
    k₁, k₂ = knots(P₁), knots(P₂)

    p₊₁, p₊₂ = p₊
    k₊₁, k₊₂ = k₊

    if (k₊₁ ≠ Knots()) && !(k₁[1] < k₊₁[1] && k₊₁[end] < k₁[end])
        error("given additional knots for refinement are out of range")
    end

    if (k₊₂ ≠ Knots()) && !(k₂[1] < k₊₂[1] && k₊₂[end] < k₂[end])
        error("given additional knots for refinement are out of range")
    end

    comment = "Refinement - p₊:$((p₊₁, p₊₂)), k₊:$((k₊₁.vector, k₊₂.vector))"
    comment = replace(comment, "Float64"=>"")
    M = refinement(M, p₊=[p₊₁, p₊₂], k₊=[k₊₁, k₊₂])
    _export(M, parent, comment=comment)
    return
end

"""
    show_knots(; index=0)

Show current knots and suggestions for knot insertions (with given index).
"""
function show_knots(; index=0)
    M = loadM(index = index)

    P = bsplinespaces(M)
    k₁, k₂ = knots.(P)
    k₁′ = unique(k₁)
    k₂′ = unique(k₂)
    msg = """
    Current knots (k₁, k₂) and suggestions for knot insertions (k₁₊, k₂₊)
    k₁: , $(k₁.vector)
    k₂: , $(k₂.vector)
    k₁₊: , $([(k₁′[i] + k₁′[i+1]) / 2 for i in 1:(length(k₁′)-1)])
    k₂₊: , $([(k₂′[i] + k₂′[i+1]) / 2 for i in 1:(length(k₂′)-1)])
    """
    @info msg
    return
end
