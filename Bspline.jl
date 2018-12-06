module Bspline

using IntervalSets
using Luxor
using ElementaryCalculus
import ParametricDraw.ChangeUnit
import ParametricDraw.BézPts
import ParametricDraw.LxrPt

export Bs1mfd, Bs2mfd, Bs, Ḃs, Bsupp, BsCoef2, BsMapping, href, pref, BsDraw, BsWrite, BsRead

mutable struct Bs1mfd
    p::Int64
    k::Array{Float64,1}
    a::Array{Float64,2}
    function Bs1mfd(p,k,a)
        if (size(a)≠((length(k)-p-1),2))
            error("dim-error")
        elseif (!issorted(k))
            error("knots not sorted")
        else
            new(p,k,a)
        end
    end
end

mutable struct Bs2mfd
    p::Array{Int64,1}
    k::Array{Array{Float64,1},1}
    a::Array{Float64,3}
    function Bs2mfd(p,k,a)
        if (length(p)≠2)
            error("p-error")
        elseif (length(k)≠2)
            error("k-error")
        elseif (size(a)≠((length.(k)-p.-1)...,2))
            error("dim-error")
        elseif (!*(issorted.(k)...))
            error("knots not sorted")
        else
            new(p,k,a)
        end
    end
end

function Bs(i::Int64,p::Int64,k,t)::Float64
    if(p==0)
        return k[i]≤t<k[i+1]||(k[i]≠k[i+1]==k[end]==t)
    else
        return (((k[i+p]-k[i]≠0) ? Bs(i,p-1,k,t)*(t-k[i])/(k[i+p]-k[i]) : 0)
        +((k[i+p+1]-k[i+1]≠0) ? Bs(i+1,p-1,k,t)*(k[i+p+1]-t)/(k[i+p+1]-k[i+1]) : 0))
    end
end

function Ḃs(i::Int64,p::Int64,k,t)::Float64
    return p*(((k[i+p]-k[i]≠0) ? Bs(i,p-1,k,t)/(k[i+p]-k[i]) : 0)
    -((k[i+p+1]-k[i+1]≠0) ? Bs(i+1,p-1,k,t)/(k[i+p+1]-k[i+1]) : 0))
end

function Bsupp(i,p,k)::ClosedInterval
    return k[i]..k[i+p+1]
end

function BsCoef2(f,p::Int64,k::Array{Float64,1};nip=25)
    n=length(k)-p-1
    A=[INT₊(t->Bs(i,p,k,t)*Bs(j,p,k,t),Bsupp(i,p,k)∩Bsupp(j,p,k),nip=nip) for i ∈ 1:n, j ∈ 1:n]
    B=[INT(t->Bs(i,p,k,t)*f(t),Bsupp(i,p,k)) for i ∈ 1:n]
    C=inv(A)*B
    return [C[I][i] for I ∈ 1:n, i ∈ 1:2]
end

function BsCoef2(f,p::Array{Int64,1},k::Array{Array{Float64,1},1};nip=25)
    n₁,n₂=n=length.(k)-p.-1
    p₁,p₂=p
    k₁,k₂=k
    A=reshape([INT2₊(u->Bs(i₁,p₁,k₁,u[1])*Bs(i₂,p₂,k₂,u[2])*Bs(j₁,p₁,k₁,u[1])*Bs(j₂,p₂,k₂,u[2]),(Bsupp(i₁,p₁,k₁)∩Bsupp(j₁,p₁,k₁),Bsupp(i₂,p₂,k₂)∩Bsupp(j₂,p₂,k₂)),nip=nip) for i₁ ∈ 1:n₁, i₂ ∈ 1:n₂, j₁ ∈ 1:n₁, j₂ ∈ 1:n₂],n₁*n₂,n₁*n₂)
    B=reshape([INT2(u->Bs(i₁,p₁,k₁,u[1])*Bs(i₂,p₂,k₂,u[2])*f(u),(Bsupp(i₁,p₁,k₁),Bsupp(i₂,p₂,k₂))) for i₁ ∈ 1:n₁, i₂ ∈ 1:n₂],n₁*n₂)
    C=reshape(inv(A)*B,n₁,n₂)
    return [C[I₁,I₂][i] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:2]
end

# function BsCoef(f,p::Int64,k::Array{Float64,1})
#     n=length(k)-p-1
#     κ=[((n+1-i)*k[i]+i*k[i+p+1])/(n+1) for i ∈ 1:n]
#     [Bs(j,p,k,κ[i]) for i ∈ 1:n, j ∈ 1:n]\map(f,κ)
# end

# function BsCoef(f,p::Array{Int64,1},k::Array{Array{Float64,1},1})
#     n=length.(k).-p-1
#     d=length(p)
#     κ=[[((n[j]+1-i)*k[j][i]+i*k[j][i+p[j]+1])/(n[j]+1) for i ∈ 1:n[j]] for j ∈ 1:d]
#     N=vcat([1],[prod(n[1:i]) for i ∈ 1:d])
#     K=Array{Array{Float64,1},1}()
#     I=Array{Array{Int64,1},1}()
#     for i ∈ 0:(N[end]-1)
#         ind=Array{Int64,1}()
#         k=0
#         for j ∈ 1:d
#             k=mod((i-k)÷N[j],n[j])
#             push!(ind,k+1)
#         end
#         push!(K,[κ[i][ind[i]] for i ∈ 1:d])
#         push!(I,[ind[i] for i ∈ 1:d])
#     end
#
#     [prod(Bs(I[j][l],p[l],k[l],K[i]) for l ∈ 1:d) for i ∈ 1:N[end], j ∈ 1:N[end]]\map(f,K)
# end

function BsMapping(B1::Bs1mfd,t)
    p=B1.p
    k=B1.k
    a=B1.a
    n=length(k)-p-1
    return sum(Bs(I,p,k,t)*a[I,:] for I ∈ 1:n)
end

function BsMapping(B2::Bs2mfd,u)
    p₁,p₂=p=B2.p
    k₁,k₂=k=B2.k
    a=B2.a
    n₁,n₂=n=length.(k)-p.-1
    return sum(Bs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])*a[I₁,I₂,:] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂)
end

# function href(B2::Bs2mfd,k₊::Array{Array{Float64,1},1})
#     p,k,a = B2.p,B2.k,B2.a
#     d=2
#     n=length.(k)-p.-1
#     pᵣ=p
#     kᵣ=[sort(convert(Array{Float64,1},vcat(k[l],k₊[l]))) for l ∈ 1:d]
#     nᵣ=length.(kᵣ)-pᵣ.-1
#     κ=[[((nᵣ[l]+1-i)/(nᵣ[l]+1))*kᵣ[l][i]+((i)/(nᵣ[l]+1))*kᵣ[l][i+pᵣ[l]+1] for i ∈ 1:nᵣ[l]] for l ∈ 1:d]
#     C=[[Bs(i,pᵣ[l],kᵣ[l],κ[l][j]) for j ∈ 1:nᵣ[l], i ∈ 1:nᵣ[l]] for l ∈ 1:d]
#     D=[C[l]\[Bs(j,p[l],k[l],κ[l][i]) for i ∈ 1:nᵣ[l], j ∈ 1:n[l]] for l ∈ 1:d]
#     aᵣ=[sum(D[1][l₁,i₁]*D[2][l₂,i₂]*a[i₁,i₂,j] for i₁ ∈ 1:n[1], i₂ ∈ 1:n[2]) for l₁ ∈ 1:nᵣ[1], l₂ ∈ 1:nᵣ[2], j ∈ 1:d]
#     return Bs2mfd(pᵣ,kᵣ,aᵣ)
# end
#
# function pref(B2::Bs2mfd,p₊::Array{Int64,1})
#     p,k,a=B2.p,B2.k,B2.a
#     d=2
#     n=length.(k)-p.-1
#     pᵣ=p+p₊
#     k₊=[repeat(DelDpl(k[l]),inner=p₊[l]) for l ∈ 1:d]
#     kᵣ=[sort(convert(Array{Float64,1},vcat(k[l],k₊[l]))) for l ∈ 1:d]
#     nᵣ=length.(kᵣ)-pᵣ.-1
#     κ=[[((nᵣ[l]+1-i)/(nᵣ[l]+1))*kᵣ[l][i]+((i)/(nᵣ[l]+1))*kᵣ[l][i+pᵣ[l]+1] for i ∈ 1:nᵣ[l]] for l ∈ 1:d]
#     C=[[Bs(i,pᵣ[l],kᵣ[l],κ[l][j]) for j ∈ 1:nᵣ[l], i ∈ 1:nᵣ[l]] for l ∈ 1:d]
#     D=[C[l]\[Bs(j,p[l],k[l],κ[l][i]) for i ∈ 1:nᵣ[l], j ∈ 1:n[l]] for l ∈ 1:d]
#     aᵣ=[sum(D[1][l₁,i₁]*D[2][l₂,i₂]*a[i₁,i₂,j] for i₁ ∈ 1:n[1], i₂ ∈ 1:n[2]) for l₁ ∈ 1:nᵣ[1], l₂ ∈ 1:nᵣ[2], j ∈ 1:d]
#     return Bs2mfd(pᵣ,kᵣ,aᵣ)
# end

function href(B2::Bs2mfd,k₊::Array{Array{Float64,1},1};nip=25)
    p,k,a = B2.p,B2.k,B2.a
    d=2
    n=length.(k)-p.-1
    pᵣ=p
    kᵣ=[sort(convert(Array{Float64,1},vcat(k[l],k₊[l]))) for l ∈ 1:d]
    aᵣ=BsCoef2(u->BsMapping(B2,u),pᵣ,kᵣ,nip=nip)
    return Bs2mfd(pᵣ,kᵣ,aᵣ)
end

function pref(B2::Bs2mfd,p₊::Array{Int64,1};nip=25)
    p,k,a=B2.p,B2.k,B2.a
    d=2
    n=length.(k)-p.-1
    pᵣ=p+p₊
    k₊=[repeat(DelDpl(k[l]),inner=p₊[l]) for l ∈ 1:d]
    kᵣ=[sort(convert(Array{Float64,1},vcat(k[l],k₊[l]))) for l ∈ 1:d]
    aᵣ=BsCoef2(u->BsMapping(B2,u),pᵣ,kᵣ,nip=nip)
    return Bs2mfd(pᵣ,kᵣ,aᵣ)
end

function BsDraw(B1::Bs1mfd;filename="BsplineCurve.svg",up=5,down=-5,right=5,left=-5,zoom=1,unitlength=(100,"pt"),points=true)
    step, unit=(unitlength[1]*zoom,unitlength[2])
    Drawing(step*(right-left),step*(up-down),filename)

    Luxor.origin(-step*left,step*up)
    setline(zoom)
    background("white")

    p=B1.p
    k=B1.k
    a=B1.a
    n=length(k)-p-1
    𝒑(t)=BsMapping(B1,t)

    K=DelDpl(k[1+p:end-p])
    N=length(K)-1

    sethue("red")
    drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(𝒑,K[i],K[i+1]))...) for i ∈ 1:N]),:stroke)

    if (points)
        sethue("black")
        setline(zoom)
        CtrlPts=[LxrPt(a[i,:],step) for i ∈ 1:size(a)[1]]
        map(p->circle(p,3*zoom,:fill), CtrlPts)
        poly(CtrlPts[:], :stroke)
    end

    finish()
    ChangeUnit(filename,"pt",unit)
    return nothing
end

function BsDraw(B2::Bs2mfd;filename="BsplineSurface.svg",up=5,down=-5,right=5,left=-5,zoom=1,mesh=(10,10),unitlength=(100,"pt"),points=true)
    step, unit=(unitlength[1]*zoom,unitlength[2])
    Drawing(step*(right-left),step*(up-down),filename)

    Luxor.origin(-step*left,step*up)
    setline(zoom)
    background("white")

    p₁,p₂=p=B2.p
    k₁,k₂=k=B2.k
    a=B2.a
    n₁,n₂=n=length.(k)-p.-1
    𝒑(u)=BsMapping(B2,u)

    K₁,K₂=K=[DelDpl(k[i][1+p[i]:end-p[i]]) for i ∈ 1:2]
    N₁,N₂=length.(K).-1
    m₁,m₂=mesh

    sethue(1,.5,.5)
    drawbezierpath(BezierPath(vcat(
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₁->𝒑([u₁,K₂[1]]),K₁[i],K₁[i+1]))...) for i ∈ 1:N₁],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₂->𝒑([K₁[end],u₂]),K₂[i],K₂[i+1]))...) for i ∈ 1:N₂],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₁->𝒑([u₁,K₂[end]]),K₁[end-i+1],K₁[end-i]))...) for i ∈ 1:N₁],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₂->𝒑([K₁[1],u₂]),K₂[end-i+1],K₂[end-i]))...) for i ∈ 1:N₂]
    )),:fill,close=true)

    sethue("red")
    for u₁ ∈ range(K₁[1],stop=K₁[end],length=m₁+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₂->𝒑([u₁,u₂]),K₂[i],K₂[i+1]))...) for i ∈ 1:N₂]),:stroke)
    end
    for u₂ ∈ range(K₂[1],stop=K₂[end],length=m₂+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u₁->𝒑([u₁,u₂]),K₁[i],K₁[i+1]))...) for i ∈ 1:N₁]),:stroke)
    end

    if (points)
        sethue("black")
        setline(zoom)
        CtrlPts=[LxrPt(a[i,j,:],step) for i ∈ 1:size(a)[1], j ∈ 1:size(a)[2]]
        map(p->circle(p,3*zoom,:fill), CtrlPts)
        for i ∈ 1:n₁
            poly(CtrlPts[i,:], :stroke)
        end
        for j ∈ 1:n₂
            poly(CtrlPts[:,j], :stroke)
        end
    end

    finish()
    ChangeUnit(filename,"pt",unit)
    return nothing
end

# function BsWrite(B2::Bs2mfd,i::Int64;filename="BsplineData.h5")
#     #ShapeName
#     p=B2.p
#     k₁,k₂=B2.k
#     a=B2.a
#     h5write(filename,"p-"*string(i),p)
#     h5write(filename,"k₁-"*string(i),k₁)
#     h5write(filename,"k₂-"*string(i),k₂)
#     h5write(filename,"a-"*string(i),a)
#     return nothing
# end
#
# function BsRead(i::Int64;filename="BsplineData.h5")
#     p=h5read(filename,"p-"*string(i))
#     k₁=h5read(filename,"k₁-"*string(i))
#     k₂=h5read(filename,"k₂-"*string(i))
#     a=h5read(filename,"a-"*string(i))
#     return Bs2mfd(p,[k₁,k₂],a)
# end

end
