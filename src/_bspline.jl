function _arrayofvector2array(a::AbstractArray{SVector{2,Float64},2})
    n1, n2 = size(a)
    a_2dim = [a[i1, i2][j] for i1 in 1:n1, i2 in 1:n2, j in 1:2]
    return a_2dim
end

function _array2arrayofvector(a::Array{<:Real,3})
    n1, n2 = size(a)
    a_vec = [SVector{2}(a[i1, i2, :]) for i1 in 1:n1, i2 in 1:n2]
    return a_vec
end

"""
Affine transform of control points.
"""
function _affine(𝒂::Matrix{<:SVector}, A::SMatrix{2,2}, b::SVector{2})
    # x'=Ax+b
    n₁, n₂ = size(𝒂)
    return [(A * 𝒂[I₁, I₂] + b) for I₁ in 1:n₁, I₂ in 1:n₂]
end

function _rotate(𝒂::Matrix{<:SVector})
    n₁, n₂ = size(𝒂)
    ind0 = [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2]
    ind1 = ind0 - [0, 1]
    v = 𝒂[ind1...] - 𝒂[ind0...]
    R = -(@SMatrix [v[2] -v[1]; v[1] v[2]]) / norm(v)
    return _affine(𝒂, R, SVector(0.0, 0.0))
end

function _center(𝒂::Matrix{<:SVector})
    xs = [p[1] for p in 𝒂]
    ys = [p[2] for p in 𝒂]
    x_min = minimum(xs)
    x_max = maximum(xs)
    y_min = minimum(ys)
    y_max = maximum(ys)
    x = (x_min + x_max) / 2
    y = (y_min + y_max) / 2
    return _affine(𝒂, one(SMatrix{2,2}), -SVector(x, y))
end

function _positioning(𝒂::Matrix{<:SVector})
    return _center(_rotate(𝒂))
end

function _positioning(M::BSplineManifold{2})
    Ps = bsplinespaces(M)
    𝒂 = controlpoints(M)
    𝒂′ = _positioning(𝒂)
    return BSplineManifold(𝒂′, Ps)
end

"""
    refinement!(steptree; p₊::Tuple{Int,Int}=(0, 0), k₊::Tuple{AbstractKnotVector,AbstractKnotVector}=(EmptyKnotVector(),EmptyKnotVector()), parent::Int=0)

Compute a refinement of the B-spline manifold
"""
function refinement!(steptree, parent::Int = 0; p₊ = (0, 0), k₊ = (EmptyKnotVector(), EmptyKnotVector()))
    parent = _validindex(steptree, parent)
    M = loadM(steptree, index = parent)

    P₁, P₂ = bsplinespaces(M)
    k₁, k₂ = knotvector(P₁), knotvector(P₂)

    p₊₁, p₊₂ = p₊
    k₊₁, k₊₂ = k₊

    if !iszero(k₊₁) && !(k₁[1] < k₊₁[1] && k₊₁[end] < k₁[end])
        error("given additional knots for refinement are out of range")
    end

    if !iszero(k₊₂) && !(k₂[1] < k₊₂[1] && k₊₂[end] < k₂[end])
        error("given additional knots for refinement are out of range")
    end

    comment = "Refinement - p₊:$((p₊₁, p₊₂)), k₊:$((BasicBSpline._vec(k₊₁), BasicBSpline._vec(k₊₂)))"
    comment = replace(comment, "Float64" => "")
    M = refinement_I(M, (Val(p₊₁), Val(p₊₂)), (k₊₁, k₊₂))
    info = Dict(["type" => "refinement"])
    step = Step(M, comment, info)
    addstep!(steptree, step, parent)
end

function suggest_knotvector(steptree; index=0)
    M = loadM(steptree, index = index)

    P = bsplinespaces(M)
    k₁, k₂ = knotvector.(P)
    k₁′ = unique(k₁)
    k₂′ = unique(k₂)
    k₁₊ = KnotVector([(k₁′[i] + k₁′[i+1]) / 2 for i in 1:(length(k₁′)-1)])
    k₂₊ = KnotVector([(k₂′[i] + k₂′[i+1]) / 2 for i in 1:(length(k₂′)-1)])
    return k₁₊, k₂₊
end

"""
    show_knotvector(::StepTree; index=0)

Show current knotvector and suggestions for knot insertions (with given index).
"""
function show_knotvector(steptree; index = 0)
    M = loadM(steptree, index = index)

    P = bsplinespaces(M)
    k₁, k₂ = knotvector.(P)
    k₁₊, k₂₊ = suggest_knotvector(steptree, index=index)
    msg = """
    Current knotvectors (k₁, k₂) and suggestions for knot insertions (k₁₊, k₂₊)
    k₁: $(BasicBSpline._vec(k₁))
    k₂: $(BasicBSpline._vec(k₂))
    k₁₊: $(BasicBSpline._vec(k₁₊))
    k₂₊: $(BasicBSpline._vec(k₂₊))
    """
    @info msg
    return
end

function integrate(C::BSplineManifold{1})
    a = controlpoints(C)
    P = bsplinespaces(C)[1]
    p = degree(P)
    k = knotvector(P)
    k′ = k + k[[begin, end]]
    p′ = p+1
    P′ = BSplineSpace{p′}(k′)
    A = [ifelse(i≤j, 0.0, (k′[p′+j+1]-k′[j+1])/(p′)) for i in 1:dim(P′), j in 1:dim(P)]
    return BSplineManifold(A*a, P′)
end

function _interpolate2(ts::AbstractVector{<:Real}, fs::AbstractVector{T}, f′0::T) where T
    # Quadric open B-spline space
    p = 2
    k = KnotVector(ts) + KnotVector([ts[1],ts[end]]) * p
    P = BSplineSpace{p}(k)

    # dimensions
    m = length(ts)
    n = dim(P)

    # The interpolant function has a f''=0 property at bounds.
    dP = BSplineDerivativeSpace{1}(P)
    d0 = [bsplinebasis(dP,j,ts[1]) for j in 1:n]

    # Compute the interpolant function (1-dim B-spline manifold)
    M = [bsplinebasis(P,j,ts[i]) for i in 1:m, j in 1:n]
    M = vcat(d0', M)
    y = vcat([f′0], fs)
    return BSplineManifold(inv(M)*y, P)
end

function _merge(manifolds::Vector{<:BSplineManifold{2, p}}) where p
    # Assume all B-spline manifolds have open knot vectors.
    p₁, p₂ = p

    k₁ = copy(knotvector(bsplinespaces(manifolds[1])[1]))
    k₂ = knotvector(bsplinespaces(manifolds[1])[2])
    for i in 2:length(manifolds)
        pop!(k₁.vector)
        k₁ += knotvector(bsplinespaces(manifolds[i])[1])[p₁+2:end]
    end
    P₁ = BSplineSpace{p₁}(k₁)
    P₂ = BSplineSpace{p₂}(k₂)

    𝒂 = controlpoints(manifolds[1])
    for i in 2:length(manifolds)
        _𝒂 = controlpoints(manifolds[i])
        v = 𝒂[end,:]
        _v = _𝒂[1,:]
        Δ = v[end] - v[1]
        _Δ = _v[end] - _v[1]
        a = dot(Δ, _Δ)
        b = cross(Δ, _Δ)
        r = (@SMatrix [a b;-b a]) / norm([a,b])
        _w = [r*p for p in _v]
        c = sum(v)/length(v)
        _c = sum(_w)/length(_w)
        _𝒂 = [r*p-_c+c for p in _𝒂]
        𝒂 = vcat(𝒂[1:end-1, :], (𝒂[end:end, :]+_𝒂[1:1, :])/2, _𝒂[2:end, :])
    end
    return BSplineManifold(𝒂, P₁, P₂)
end
