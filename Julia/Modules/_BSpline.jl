using Luxor
import ParametricDraw.ChangeUnit
import ParametricDraw.BézPts
import ParametricDraw.LxrPt

# BSpline
function FittingBSpline(f, P::BSplineSpace; nip=NIP) # 1-dimensional
    p=P.degree
    k=P.knots
    D=k[1+p]..k[end-p]
    function a(i,j)
        D′=(max(k[i],k[j])..min(k[i+p+1],k[j+p+1])) ∩ D
        if width(D′)==0
            return 0
        else
            return GaussianQuadrature(t->BSplineBasis(i,P,t)*BSplineBasis(j,P,t), D′)
        end
    end
    n=dim(P)
    A=[a(i,j) for i ∈ 1:n, j ∈ 1:n]
    b=[GaussianQuadrature(t->BSplineBasis(i,P,t)*f(t), ((k[i]..k[i+p+1]) ∩ D)) for i ∈ 1:n]
    return inv(A)*b
end

function N′(P₁::BSplineSpace,P₂::BSplineSpace,I₁,I₂,i,u)::Float64
    if i==1
        return BSplineBasis′(I₁,P₁,u[1])*BSplineBasis(I₂,P₂,u[2])
    else
        return BSplineBasis(I₁,P₁,u[1])*BSplineBasis′(I₂,P₂,u[2])
    end
end


function BSplineSvg2(M::BSplineManifold; filename="Bspline.svg", up=5, down=-5, right=5, left=-5, zoom=1, mesh=(10,10), unitlength=(100,"pt"), points=true)
    step, unit = (unitlength[1]*zoom,unitlength[2])
    Drawing(step*(right-left),step*(up-down),filename)

    Luxor.origin(-step*left,step*up)
    setline(zoom)
    background("white")

    p¹,p² = p = [M.bsplinespaces[i].degree for i ∈ 1:2]
    k¹,k² = k = [M.bsplinespaces[i].knots for i ∈ 1:2]
    𝒂 = M.controlpoints

    n¹,n² = n = length.(k)-p.-1
    𝒑(u) = Mapping(M,u)

    K¹,K² = K = [unique(k[i][1+p[i]:end-p[i]]) for i ∈ 1:2]
    N¹,N² = length.(K).-1
    m¹,m² = mesh

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

    if points
        sethue("black")
        setline(zoom)
        CtrlPts = [LxrPt(𝒂[i,j,:],step) for i ∈ 1:size(𝒂)[1], j ∈ 1:size(𝒂)[2]]
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

function DrawBSpline(M::BSplineManifold; filename="Bspline.svg", up=5, down=-5, right=5, left=-5, zoom=1, mesh=(10,10), unitlength=100, points=true)
    step = unitlength
    p¹,p² = p = [M.bsplinespaces[i].degree for i ∈ 1:2]
    k¹,k² = k = [M.bsplinespaces[i].knots for i ∈ 1:2]
    𝒂 = M.controlpoints
    n¹,n² = n = length.(k)-p.-1
    𝒑(u) = Mapping(M,u)

    K¹,K² = K = [unique(k[i][1+p[i]:end-p[i]]) for i ∈ 1:2]
    N¹,N² = length.(K).-1
    m¹,m² = mesh

    Drawing(step*(right-left),step*(up-down),filename)
    Luxor.origin(-step*left,step*up)
    setline(zoom)
    background("white")

    sethue(1,.5,.5) # Pale Red
    drawbezierpath(BezierPath(vcat(
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,K²[1]]),K¹[i],K¹[i+1]))...) for i ∈ 1:N¹],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([K¹[end],u²]),K²[i],K²[i+1]))...) for i ∈ 1:N²],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,K²[end]]),K¹[end-i+1],K¹[end-i]))...) for i ∈ 1:N¹],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([K¹[1],u²]),K²[end-i+1],K²[end-i]))...) for i ∈ 1:N²]
    )),:fill,close=true)

    sethue("red") # Red
    for u¹ ∈ range(K¹[1],stop=K¹[end],length=m¹+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([u¹,u²]),K²[i],K²[i+1]))...) for i ∈ 1:N²]),:stroke)
    end
    for u² ∈ range(K²[1],stop=K²[end],length=m²+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,u²]),K¹[i],K¹[i+1]))...) for i ∈ 1:N¹]),:stroke)
    end

    if points
        sethue(.1,.1,.1) # Dark Gray
        setline(zoom)
        CtrlPts = [LxrPt(𝒂[i,j,:],step) for i ∈ 1:size(𝒂)[1], j ∈ 1:size(𝒂)[2]]
        map(p->circle(p,3*zoom,:fill), CtrlPts)

        sethue(.3,.3,.3) # Light Gray
        for i ∈ 1:n¹
            poly(CtrlPts[i,:], :stroke)
        end
        for j ∈ 1:n²
            poly(CtrlPts[:,j], :stroke)
        end
    end
    finish()

    return nothing
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

function Positioning(M::BSplineManifold)::BSplineManifold # 制御点の位置調整
    𝒫s = M.bsplinespaces
    𝒂 = M.controlpoints
    if length(𝒫s) ≠ d
        error("dimension does not match")
    end

    p¹, p² = p = [M.bsplinespaces[i].degree for i ∈ 1:d]
    k¹, k² = k = [M.bsplinespaces[i].knots for i ∈ 1:d]

    n₁, n₂, _ = size(𝒂)
    𝒂′ = Positioning(𝒂)
    return BSplineManifold(𝒫s,𝒂′)
end

export SplineRefinement
function SplineRefinement( ;p₊::Array{Int,1}=[0,0], k₊::Array{Knots,1}=[Knots([]),Knots([])], parent::Int=0)
    parent=Parent(parent)
    M=loadM(index=parent)

    P₁,P₂=P=M.bsplinespaces
    p₁,p₂=p=P₁.degree,P₂.degree
    k₁,k₂=k=P₁.knots,P₂.knots
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
    M=BSpline.Refinement(M,p₊=p₊,k₊=k₊)
    Export(M,parent,comment=comment)
    return nothing
end

export ShowKnots
function ShowKnots( ;index=0)
    M=loadM(index=index)

    P₁,P₂=M.bsplinespaces
    p₁,p₂=P₁.degree,P₂.degree
    k₁,k₂=P₁.knots,P₂.knots
    println("k₁: ",k₁.vector)
    println("k₂: ",k₂.vector)
    println("Suggestion:")
    k₁′=unique(k₁)
    k₂′=unique(k₂)
    println("k₁₊: ",[(k₁′[i]+k₁′[i+1])/2 for i ∈ 1:(length(k₁′)-1)])
    println("k₂₊: ",[(k₂′[i]+k₂′[i+1])/2 for i ∈ 1:(length(k₂′)-1)])
    return nothing
end
