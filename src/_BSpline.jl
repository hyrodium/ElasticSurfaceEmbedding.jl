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
function affine(𝒂::Array{Float64,3}, A::Array{Float64,2}, b::Array{Float64,1})::Array{Float64,3}
    # x'=Ax+b
    n₁, n₂, d = size(𝒂)
    return [(A*𝒂[I₁, I₂, :]+b)[i] for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:d]
end

function Positioning(𝒂::Array{Float64,3})::Array{Float64,3} # 制御点の位置調整
    n₁, n₂, _ = size(𝒂)
    ind0 = [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2]
    ind1 = ind0 - [0, 1]
    v = 𝒂[ind1..., :] - 𝒂[ind0..., :]
    R = -[v[2] -v[1]; v[1] v[2]] / norm(v)
    return affine(𝒂, R, -R * 𝒂[ind0..., :])
end

function Positioning(M::AbstractBSplineManifold) # 制御点の位置調整
    Ps = collect(bsplinespaces(M))
    𝒂 = controlpoints(M)
    if length(Ps) ≠ 2
        error("dimension does not match")
    end

    𝒂′ = Positioning(𝒂)
    return typeof(M)(Ps, 𝒂′)
end

export spline_refinement
function spline_refinement(; p₊::Array{Int,1}=[0, 0], k₊::Array{Knots,1}=[Knots(), Knots()], parent::Int=0)
    parent = Parent(parent)
    M = loadM(index = parent)

    P₁, P₂ = collect(bsplinespaces(M))
    k₁, k₂ = knots(P₁), knots(P₂)

    k₊₁, k₊₂ = k₊

    if (k₊₁ ≠ Knots()) && !(k₁[1] < k₊₁[1] && k₊₁[end] < k₁[end])
        error("given additional knots for refinement are out of range")
    end

    if (k₊₂ ≠ Knots()) && !(k₂[1] < k₊₂[1] && k₊₂[end] < k₂[end])
        error("given additional knots for refinement are out of range")
    end

    comment = "Refinement - p₊:" * string(p₊) * ", k₊:" * string([k₊₁.vector, k₊₂.vector])
    M = refinement(M, p₊ = p₊, k₊ = k₊)
    Export(M, parent, comment = comment)
    return
end

export print_knots
function print_knots(; index = 0)
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
