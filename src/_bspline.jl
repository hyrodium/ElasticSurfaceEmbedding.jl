function _arrayofvector2array(a::AbstractArray{SVector{2,Float64},2})
    n1,n2 = size(a)
    a_2dim = [a[i1,i2][j] for i1 in 1:n1, i2 in 1:n2, j in 1:2]
    return a_2dim
end

function _array2arrayofvector(a::Array{<:Real,3})
    n1,n2 = size(a)
    a_vec = [SVector{2}(a[i1,i2,:]) for i1 in 1:n1, i2 in 1:n2]
    return a_vec
end

"""
Affine transform of control points.
"""
function _affine(𝒂::Matrix{<:SVector}, A::SMatrix{2,2}, b::SVector{2})
    # x'=Ax+b
    n₁, n₂ = size(𝒂)
    return [(A*𝒂[I₁,I₂]+b) for I₁ in 1:n₁, I₂ in 1:n₂]
end

function _rotate(𝒂::Matrix{<:SVector})
    n₁, n₂ = size(𝒂)
    ind0 = [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2]
    ind1 = ind0 - [0, 1]
    v = 𝒂[ind1...] - 𝒂[ind0...]
    R = - (@SMatrix [v[2] -v[1]; v[1] v[2]]) / norm(v)
    return _affine(𝒂, R, SVector(0.0, 0.0))
end

function _center(𝒂::Matrix{<:SVector})
    xs = [p[1] for p in 𝒂]
    ys = [p[2] for p in 𝒂]
    x_min = minimum(xs)
    x_max = maximum(xs)
    y_min = minimum(ys)
    y_max = maximum(ys)
    x = (x_min+x_max)/2
    y = (y_min+y_max)/2
    return _affine(𝒂, one(SMatrix{2,2}), -SVector(x,y))
end

function _positioning(𝒂::Matrix{<:SVector})
    return _center(_rotate(𝒂))
end

function _positioning(M::BSplineManifold{2})
    Ps = bsplinespaces(M)
    𝒂 = controlpoints(M)
    𝒂′ = _positioning(𝒂)
    return BSplineManifold(𝒂′,Ps)
end

"""
    refinement!(allsteps; p₊::Tuple{Int,Int}=(0, 0), k₊::Tuple{AbstractKnotVector,AbstractKnotVector}=(EmptyKnotVector(),EmptyKnotVector()), parent::Int=0)

Compute a refinement of the B-spline manifold
"""
function refinement!(allsteps; p₊=(0,0), k₊=(EmptyKnotVector(),EmptyKnotVector()), parent::Int=0)
    parent = _validindex(allsteps, parent)
    M = loadM(allsteps, index=parent)

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
    comment = replace(comment, "Float64"=>"")
    M = refinement(M, (Val(p₊₁), Val(p₊₂)), (k₊₁, k₊₂))
    step = Step(M, comment)
    addstep!(allsteps, step, parent)
end

"""
    show_knotvector(; index=0)

Show current knotvector and suggestions for knot insertions (with given index).
"""
function show_knotvector(allsteps; index=0)
    M = loadM(allsteps, index=index)

    P = bsplinespaces(M)
    k₁, k₂ = knotvector.(P)
    k₁′ = unique(k₁)
    k₂′ = unique(k₂)
    msg = """
    Current knotvectors (k₁, k₂) and suggestions for knot insertions (k₁₊, k₂₊)
    k₁: , $BasicBSpline._vec((k₁))
    k₂: , $BasicBSpline._vec((k₂))
    k₁₊: , $([(k₁′[i] + k₁′[i+1]) / 2 for i in 1:(length(k₁′)-1)])
    k₂₊: , $([(k₂′[i] + k₂′[i+1]) / 2 for i in 1:(length(k₂′)-1)])
    """
    @info msg
    return
end
