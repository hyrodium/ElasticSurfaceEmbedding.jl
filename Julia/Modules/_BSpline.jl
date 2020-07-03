using Luxor
import ParametricDraw.ChangeUnit
import ParametricDraw.BézPts
import ParametricDraw.LxrPt

# BSpline
function FittingBSpline(f, P::FastBSplineSpace{p}; nip=NIP) where p # 1-dimensional
    k=knots(P)
    D=k[1+p]..k[end-p]
    function a(i,j)
        D′=(max(k[i],k[j])..min(k[i+p+1],k[j+p+1])) ∩ D
        if width(D′)==0
            return 0.0
        else
            return GaussianQuadrature(t->bsplinebasis(i,P,t)*bsplinebasis(j,P,t), D′)
        end
    end
    n=dim(P)
    A=[a(i,j) for i ∈ 1:n, j ∈ 1:n]
    b=[GaussianQuadrature(t->bsplinebasis(i,P,t)*f(t), ((k[i]..k[i+p+1]) ∩ D)) for i ∈ 1:n]
    return inv(A)*b
end

function N′(P₁::FastBSplineSpace, P₂::FastBSplineSpace, I₁, I₂, i, u)::Float64
    if i==1
        return bsplinebasis′(I₁,P₁,u[1])*bsplinebasis(I₂,P₂,u[2])
    else
        return bsplinebasis(I₁,P₁,u[1])*bsplinebasis′(I₂,P₂,u[2])
    end
end

"""
Affine transform of control points.
"""
function affine(𝒂::Array{Float64,3},A::Array{Float64,2},b::Array{Float64,1})::Array{Float64,3}
    #x'=Ax+b
    n₁, n₂, d = size(𝒂)
    return [(A*𝒂[I₁,I₂,:]+b)[i] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d]
end

function Positioning(𝒂::Array{Float64,3})::Array{Float64,3} # 制御点の位置調整
    n₁, n₂, _ = size(𝒂)
    ind0 = [(n₁+1)÷2,(n₂+1)÷2]
    ind1 = ind0-[0,1]
    v = 𝒂[ind1...,:]-𝒂[ind0...,:]
    R = -[v[2] -v[1];v[1] v[2]]/norm(v)
    return affine(𝒂,R,-R*𝒂[ind0...,:])
end

function Positioning(M::FastBSplineManifold)::FastBSplineManifold # 制御点の位置調整
    Ps = M.bsplinespaces
    𝒂 = M.controlpoints
    if length(Ps) ≠ d
        error("dimension does not match")
    end

    p¹, p² = p = [degree(M.bsplinespaces[i]) for i ∈ 1:d]
    k¹, k² = k = [knots(M.bsplinespaces[i]) for i ∈ 1:d]

    n₁, n₂, _ = size(𝒂)
    𝒂′ = Positioning(𝒂)
    return FastBSplineManifold(Ps,𝒂′)
end

export SplineRefinement
function SplineRefinement(; p₊::Array{Int,1}=[0,0], k₊::Array{Knots,1}=[Knots([]),Knots([])], parent::Int=0)
    parent=Parent(parent)
    M=loadM(index=parent)

    P₁,P₂=P=M.bsplinespaces
    p₁,p₂=p=degree.(P)
    k₁,k₂=k=knots.(P)
    D₁,D₂=D=k₁[1+p₁]..k₁[end-p₁],k₂[1+p₂]..k₂[end-p₂]
    n₁,n₂=n=dim.(P)

    k₊₁,k₊₂=k₊

    if (k₊₁ ≠ Knots([])) && !( k₁[1]<k₊₁[1] && k₊₁[end]<k₁[end] )
        error("given additional knots for refinement are out of range")
    end

    if (k₊₂ ≠ Knots([])) && !( k₂[1]<k₊₂[1] && k₊₂[end]<k₂[end] )
        error("given additional knots for refinement are out of range")
    end

    comment="Refinement - p₊:"*string(p₊)*", k₊:"*string([k₊₁.vector, k₊₂.vector])
    M=refinement(M,p₊=p₊,k₊=k₊)
    Export(M,parent,comment=comment)
    return nothing
end

export ShowKnots
function ShowKnots(; index=0)
    M=loadM(index=index)

    P₁,P₂=P=M.bsplinespaces
    p₁,p₂=degree.(P)
    k₁,k₂=knots.(P)
    println("k₁: ",k₁.vector)
    println("k₂: ",k₂.vector)
    println("Suggestion:")
    k₁′=unique(k₁)
    k₂′=unique(k₂)
    println("k₁₊: ",[(k₁′[i]+k₁′[i+1])/2 for i ∈ 1:(length(k₁′)-1)])
    println("k₂₊: ",[(k₂′[i]+k₂′[i+1])/2 for i ∈ 1:(length(k₂′)-1)])
    return nothing
end
