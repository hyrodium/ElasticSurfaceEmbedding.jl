# module Bspline
#
# using IntervalSets
# using Luxor
# using ElementaryCalculus
# import ParametricDraw.ChangeUnit
# import ParametricDraw.BézPts
# import ParametricDraw.LxrPt
#
# export Bs1mfd, Bs2mfd, Bs, Ḃs, Bsupp, BsCoef2, BsMapping, href, pref, BsDraw, BsWrite, BsRead

push!(LOAD_PATH, "/home/hyrodium/Git/I4SM-julia/")
using BenchmarkTools
using IntervalSets
using Luxor
import ParametricDraw.ChangeUnit
import ParametricDraw.BézPts
import ParametricDraw.LxrPt

# Knots
struct Knots
    vector :: Array{Float64,1}
    function Knots(vector)
        new(sort(vector))
    end
end

Base.:+(k₁::Knots, k₂::Knots) = Knots(sort([k₁.vector...,k₂.vector...]))
Base.:*(p₊::Int, k::Knots) = (
        if (p₊==0)
            Knots([])
        elseif (p₊>0)
            sum(k for _ ∈ 1:p₊)
        else
            error("p₊ must be non-negative")
        end
    )

Base.in(r::Real, k::Knots) = in(r,k.vector)
Base.getindex(k::Knots, i::Int) = k.vector[i]
Base.getindex(k::Knots, v::AbstractArray{Int64,1}) = Knots(k.vector[v])
Base.length(k::Knots) = length(k.vector)
♯(k::Knots) = length(k::Knots)
Base.firstindex(k) = 1
Base.lastindex(k) = length(k)
Base.unique(k::Knots) = Knots(unique(k.vector))

function Base.:⊆(k::Knots, k′::Knots)
    K′=copy(k′.vector)
    for kᵢ ∈ k.vector
        i=findfirst(x->x==kᵢ,K′)
        if (i isa Nothing)
            return false
        end
        deleteat!(K′,i)
    end
    return true
end


# B-Spline Space
struct BSplineSpace
    degree::Int
    knots::Knots
    function BSplineSpace(degree::Int, knots::Knots)
        if (degree < 0)
            error("degree of polynominal must be non-negative")
        end
        new(degree,knots)
    end
end

const 𝒫 = BSplineSpace

function dim(bsplinespace::BSplineSpace)
    p=bsplinespace.degree
    k=bsplinespace.knots
    return length(k)-p-1
end

function Base.:⊆(P::BSplineSpace, P′::BSplineSpace)
    p=P.degree
    k=P.knots
    p′=P′.degree
    k′=P′.knots
    p₊=p′-p

    return (k+p₊*unique(k) ⊆ k′) && p₊ ≥ 0
end

function Base.iszero(P::BSplineSpace)
    p=P.degree
    k=P.knots
    n=dim(P)
    return [k[i]==k[i+p+1] for i ∈ 1:n]
end

# B-Spline functions
function BSplineBasis₊₀(P::BSplineSpace, t)::Array{Float64,1}
    p=P.degree
    k=P.knots

    n=length(k)-p-1
    if (p==0)
        return [k[i] ≤ t < k[i+1] for i ∈ 1:n]
    end
    K=[ifelse(k[i+p]==k[i],0,(t-k[i])/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B=BSplineBasis₊₀(𝒫(p-1,k),t)
    return [K[i]*B[i]+(1-K[i+1])*B[i+1] for i ∈ 1:n]
end

function BSplineBasis₋₀(P::BSplineSpace, t)::Array{Float64,1}
    p=P.degree
    k=P.knots

    n=length(k)-p-1
    if (p==0)
        return [k[i] < t ≤ k[i+1] for i ∈ 1:n]
    end
    K=[ifelse(k[i+p]==k[i],0,(t-k[i])/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B=BSplineBasis₋₀(𝒫(p-1,k),t)
    return [K[i]*B[i]+(1-K[i+1])*B[i+1] for i ∈ 1:n]
end

function BSplineBasis(P::BSplineSpace, t)::Array{Float64,1}
    p=P.degree
    k=P.knots

    n=length(k)-p-1
    if (p==0)
        return [k[i] ≤ t < k[i+1] || (k[i] ≠ k[i+1] == k[end] == t) for i ∈ 1:n]
    end
    K=[ifelse(k[i+p]==k[i],0,(t-k[i])/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B=BSplineBasis(𝒫(p-1,k),t)
    return [K[i]*B[i]+(1-K[i+1])*B[i+1] for i ∈ 1:n]
end

# BSplineBasis₊₀ = BSplineBasis

function BSplineBasis′(P::BSplineSpace, t)::Array{Float64,1}
    p=P.degree
    k=P.knots

    n=length(k)-p-1
    if (p==0)
        return [0.0 for _ ∈ 1:n]
    end
    K=[ifelse(k[i+p]==k[i],0,p/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B=BSplineBasis(𝒫(p-1,k),t)
    return [K[i]*B[i]-K[i+1]*B[i+1] for i ∈ 1:n]
end

function BSplineBasis′₊₀(P::BSplineSpace, t)::Array{Float64,1}
    p=P.degree
    k=P.knots

    n=length(k)-p-1
    if (p==0)
        return [0.0 for _ ∈ 1:n]
    end
    K=[ifelse(k[i+p]==k[i],0,p/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B=BSplineBasis₊₀(𝒫(p-1,k),t)
    return [K[i]*B[i]-K[i+1]*B[i+1] for i ∈ 1:n]
end

function BSplineBasis′₋₀(P::BSplineSpace, t)::Array{Float64,1}
    p=P.degree
    k=P.knots

    n=length(k)-p-1
    if (p==0)
        return [0.0 for _ ∈ 1:n]
    end
    K=[ifelse(k[i+p]==k[i],0,p/(k[i+p]-k[i])) for i ∈ 1:n+1]
    B=BSplineBasis₋₀(𝒫(p-1,k),t)
    return [K[i]*B[i]-K[i+1]*B[i+1] for i ∈ 1:n]
end

# BSplineBasis′₊₀ = BSplineBasis′

function BSplineSupport(i::Int64, P::BSplineSpace)::ClosedInterval
    p=P.degree
    k=P.knots

    return k[i]..k[i+p+1]
end

function BSplineCoefficient(P::BSplineSpace, P′::BSplineSpace)::Array{Float64,2}
    p=P.degree
    k=P.knots
    p′=P′.degree
    k′=P′.knots
    p₊=p′-p
    if (P ⊈ P′)
        error("𝒫[p,k] ⊄ 𝒫[p′,k′]")
    end

    if (p == 0)
        n=length(k)-1
        n′=length(k′)-p₊-1
        A⁰=Float64[BSplineSupport(j,𝒫(p₊,k′)) ⊆ BSplineSupport(i,𝒫(0,k)) for i ∈ 1:n, j ∈ 1:n′]
        A⁰[:,findall(iszero(P′))].=NaN
        return A⁰
    end

    Aᵖ⁻¹=BSplineCoefficient(𝒫(p-1, k), 𝒫(p′-1, k′))
    n=dim(P)
    n′=dim(P′)
    Z=iszero(𝒫(p′-1,k′))
    W=findall(Z)
    K′=[k′[i+p′]-k′[i] for i ∈ 1:n′+1]
    K=[ifelse(k[i+p]≠k[i], 1/(k[i+p]-k[i]), 0.0) for i ∈ 1:n+1]
    Δ=(p/p′)*[K′[j]*(K[i]*Aᵖ⁻¹[i,j]-K[i+1]*Aᵖ⁻¹[i+1,j]) for i ∈ 1:n, j ∈ 1:n′+1]
    Aᵖ=zeros(n,n′)
    Aᵖ[:,1]=Δ[:,1]
    Aᵖ[:,n′]=-Δ[:,n′+1]

    if (length(W)==0)
        Q=[1:n′]
    else
        Q=[1:W[1]-1,[W[i]:W[i+1]-1 for i ∈ 1:length(W)-1]...,W[end]:n′]
    end
    l=length(Q)
    L=length.(Q)
    Ãᵖ=[Aᵖ[:,q] for q ∈ Q]

    for ȷ ∈ 2:l-1
        if (L[ȷ]==1)
            Ãᵖ[ȷ] .= NaN
        end
    end
    for ȷ ∈ 1:l-1
        if (L[ȷ] ≥ 2)
            t=k′[W[ȷ]]
            Ãᵖ[ȷ][:,end]=BSplineBasis₋₀(𝒫(p,k),t)
        end
    end
    for ȷ ∈ 2:l
        if (L[ȷ] ≥ 2)
            t=k′[W[ȷ-1]+p]
            Ãᵖ[ȷ][:,1]=BSplineBasis₊₀(𝒫(p,k),t)
        end
    end
    for ȷ ∈ 1:l
        if (L[ȷ] ≥ 3)
            r=Q[ȷ]
            A₊=copy(Ãᵖ[ȷ])
            A₋=copy(Ãᵖ[ȷ])
            for j ∈ 1:L[ȷ]-2
                A₊[:,j+1]=A₊[:,j]+Δ[:,j+r[1]]
                A₋[:,L[ȷ]-j]=A₋[:,L[ȷ]-j+1]-Δ[:,L[ȷ]-j+r[1]]
            end
            Ãᵖ[ȷ]=(A₊+A₋)/2
        end
    end
    Aᵖ=hcat(Ãᵖ...)
    return Aᵖ .* Float64[BSplineSupport(j,P′) ⊆ BSplineSupport(i,P) for i ∈ 1:n, j ∈ 1:n′]
end

function ⊗(X::Array{Float64},Y::Array{Float64})::Array{Float64}
    m=size(X)
    n=size(Y)
    reshape(reshape(X,length(X)) * reshape(Y,length(Y))', m..., n...)
end

function tensorprod(X::Array{T,1}) where T <: Array{Float64}
    n=length(X)
    # X[1] ⊗ … ⊗ X[n]
    @inbounds Y=X[1]
    for i ∈ 2:n
        @inbounds Y = Y ⊗ X[i]
    end
    return Y
end

struct BSplineManifold
    bsplinespaces::Array{BSplineSpace,1}
    controlpoints::Array{Float64}
    function BSplineManifold(bsplinespaces::Array{BSplineSpace,1}, controlpoints::Array{Float64})
        if (collect(size(controlpoints)[1:end-1]) ≠ dim.(bsplinespaces))
            error("dimension does not match")
        else
            new(bsplinespaces, controlpoints)
        end
    end
end


function Refinement(M::BSplineManifold, 𝒫s′::Array{BSplineSpace,1})
    𝒫s = M.bsplinespaces
    𝒂 = M.controlpoints
    d̂ = size(𝒂)[end]
    n = dim.(𝒫s)
    n′ = dim.(𝒫s′)
    if (prod(𝒫s .⊆ 𝒫s′))
        A = BSplineCoefficient.(𝒫s,𝒫s′)
        𝒂′ = [sum(A[1][I₁,J₁]*A[2][I₂,J₂]*𝒂[I₁,I₂,i] for I₁ ∈ 1:n[1], I₂ ∈ 1:n[2]) for J₁ ∈ 1:n′[1], J₂ ∈ 1:n′[2], i ∈ 1:d̂]
        return BSplineManifold(𝒫s′, 𝒂′)
    else
        error("𝒫[p,k] ⊄ 𝒫[p′,k′]")
    end
end

# function BSplineBasis(𝒫s::Array{BSplineSpace,1},t)
#     if (length(𝒫s)==length(t)==1)
#         return BSplineBasis(𝒫s[1],t[1])
#     elseif (length(𝒫s)==length(t)==2)
#         return BSplineBasis(𝒫s[1],t[1])*BSplineBasis(𝒫s[2],t[2])'
#     else
#         error("dimension does not match")
#     end
# end

function BSplineBasis(𝒫s::Array{BSplineSpace,1},t)
    d=length(t)
    Bs=[BSplineBasis(𝒫s[i],t[i]) for i ∈ 1:d]
    return tensorprod(Bs)
end

# function Mapping(M::BSplineManifold, t::Array{Float64,1})
#     𝒫s = M.bsplinespaces
#     𝒂 = M.controlpoints
#     d=length(𝒫s)
#     d̂=size(𝒂)[end]
#     return [sum(BSplineBasis(𝒫s,t).*𝒂[:,:,i]) for i ∈ 1:d̂]
# end

function Mapping(M::BSplineManifold, t::Array{Float64,1})
    𝒫s = M.bsplinespaces
    𝒂 = M.controlpoints
    d=length(𝒫s)
    d̂=size(𝒂)[end]
    return [sum(BSplineBasis(𝒫s,t).*𝒂[:,:,i]) for i ∈ 1:d̂]
end

function BSplineSvg(M::BSplineManifold; filename="Bspline.svg", up=5, down=-5, right=5, left=-5, zoom=1, mesh=(10,10), unitlength=(100,"pt"), points=true)
    step, unit=(unitlength[1]*zoom,unitlength[2])
    Drawing(step*(right-left),step*(up-down),filename)

    Luxor.origin(-step*left,step*up)
    setline(zoom)
    background("white")

    p¹,p²=p=[M.bsplinespaces[i].degree for i ∈ 1:2]
    k¹,k²=k=[M.bsplinespaces[i].knots for i ∈ 1:2]
    𝒂=M.controlpoints

    n¹,n²=n=length.(k)-p.-1
    𝒑(u)=Mapping(M,u)

    K¹,K²=K=[unique(k[i][1+p[i]:end-p[i]]) for i ∈ 1:2]
    N¹,N²=length.(K).-1
    m¹,m²=mesh

    sethue(1,.5,.5)
    drawbezierpath(BezierPath(vcat(
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,K²[1]]),K¹[i],K¹[i+1]))...) for i ∈ 1:N¹],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([K¹[end],u²]),K²[i],K²[i+1]))...) for i ∈ 1:N²],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,K²[end]]),K¹[end-i+1],K¹[end-i]))...) for i ∈ 1:N¹],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([K¹[1],u²]),K²[end-i+1],K²[end-i]))...) for i ∈ 1:N²]
    )),:fill,close=true)

    sethue("red")
    for u¹ ∈ range(K¹[1],stop=K¹[end],length=m¹+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([u¹,u²]),K²[i],K²[i+1]))...) for i ∈ 1:N²]),:stroke)
    end
    for u² ∈ range(K²[1],stop=K²[end],length=m²+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,u²]),K¹[i],K¹[i+1]))...) for i ∈ 1:N¹]),:stroke)
    end

    if (points)
        sethue("black")
        setline(zoom)
        CtrlPts=[LxrPt(𝒂[i,j,:],step) for i ∈ 1:size(𝒂)[1], j ∈ 1:size(𝒂)[2]]
        map(p->circle(p,3*zoom,:fill), CtrlPts)
        for i ∈ 1:n¹
            poly(CtrlPts[i,:], :stroke)
        end
        for j ∈ 1:n²
            poly(CtrlPts[:,j], :stroke)
        end
    end

    finish()
    ChangeUnit(filename,"pt",unit)
    return nothing
end







## Examples below

# k¹ = Knots([0,0,1,1])
# k² = Knots([0,0,1/2,1,1])
k¹ = Knots([0,0,1,1])
k² = Knots([0,0,1/2,1,1])
p¹ = 1
p² = 1
𝒂=Float64[ifelse(i==1, 2*I₁, 3*I₂/2) for I₁ ∈ 0:1,  I₂ ∈ 0:2, i ∈ 1:2]
𝒫s=[𝒫(p¹,k¹), 𝒫(p²,k²)]
M=BSplineManifold(𝒫s,𝒂)

[M.bsplinespaces[i].knots for i ∈ 1:2]

BSplineSvg(M)

p′¹=3
p′²=5
k′¹=k¹+(p′¹-p¹)*unique(k¹)+Knots(rand(3))
k′²=k²+(p′²-p²)*unique(k²)+Knots(rand(5))
M′=Refinement(M, [𝒫(p′¹,k′¹), 𝒫(p′²,k′²)])

BSplineSvg(M′,filename="hoge.svg")

M == M′

a=M.controlpoints
a′=M′.controlpoints

t=[rand(),2*rand()]
println(Mapping(M,t)-Mapping(M′,t))

@benchmark Mapping(M′,t)


@benchmark Mapping2(M′,t)


# exit()

# end

println(a[:,1,1])
println(a′[:,1,1])

@benchmark a′[:,1,1]
