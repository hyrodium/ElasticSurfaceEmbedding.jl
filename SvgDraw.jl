module SvgDraw

using IntervalSets
using Luxor

export SvgCurve

function BézPts(𝒑,a,b) # Bézier曲線の制御点
    𝒑(a),
    3*(𝒑(2*a/3+b/3)-(8*𝒑(a)+𝒑(b))/27)-3*(𝒑(a/3+2*b/3)-(𝒑(a)+8*𝒑(b))/27)/2,
    -3*(𝒑(2*a/3+b/3)-(8*𝒑(a)+𝒑(b))/27)/2+3*(𝒑(a/3+2*b/3)-(𝒑(a)+8*𝒑(b))/27),
    𝒑(b)
end

function LxrPt(p::Array{T,1},step) where T<:Real
    Point(step*[1,-1].*p...)
end

function ChangeUnit(name,before,after)
    ss=open(name) do file
        strn=read(file, String)
        replace(strn,Regex("(\\d+)"*before)=>SubstitutionString("\\1"*after))
    end
    open(name,"w") do file
        write(file,ss)
    end
    return ss
end

function SvgCurve(𝒑,k::Array{T,1};name="BCA.svg",up=5,down=-5,right=5,left=-5,zoom=1,unitlength=(100,"pt")) where T<:Real
    n=length(k)-1
    step, unit=(unitlength[1]*zoom,unitlength[2])
    Drawing((right-left)*step,(up-down)*step,name)

    Luxor.origin(-left*step,up*step)
    background("white")

    BézPth=BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(𝒑,k[i],k[i+1]))...) for i in 1:n])

    setline(zoom)
    sethue("red")
    drawbezierpath(BézPth, :stroke)
    finish()
    ChangeUnit(name,"pt",unit)
    return true
end

function SvgCurve(𝒑,I::ClosedInterval;name="BCA.svg",up=5,down=-5,right=5,left=-5,zoom=1,unitlength=(100,"pt")) where T<:Real
    k=collect(range(endpoints(I)...,length=60))
    n=length(k)-1
    step, unit=(unitlength[1]*zoom,unitlength[2])
    Drawing((right-left)*step,(up-down)*step,name)

    Luxor.origin(-left*step,up*step)
    background("white")

    BézPth=BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(𝒑,k[i],k[i+1]))...) for i in 1:n])

    setline(zoom)
    sethue("red")
    drawbezierpath(BézPth, :stroke)
    finish()
    ChangeUnit(name,"pt",unit)
    return true
end


function SvgSurface(𝒑,k,n;name="BSA.svg",up=5,down=-5,right=5,left=-5,step=50)
    k₁,k₂=k
    n₁,n₂=n
    K₁,K₂=range(minimum(k₁),maximum(k₁),length=n₁+1), range(minimum(k₂),maximum(k₂),length=n₂+1)
    N₁,N₂=length(k₁),length(k₂)

    Drawing(step*(right-left),step*(up-down),name)
    Luxor.origin(-step*left,step*up)
    background("white")
    sethue("blue")

    for j in 1:(n₂+1)
        move(step*Point([1,-1].*𝒑([k₁[1],K₂[j]])...))
        for i in 1:(N₁-1)
            a₁,a₂,a₃,a₄=BézPts(t->𝒑([t,K₂[j]]),k₁[i],k₁[i+1])
            curve(step*Point([1,-1].*a₂...),step*Point([1,-1].*a₃...),step*Point([1,-1].*a₄...))
        end
    end
    for j in 1:(n₁+1)
        move(step*Point([1,-1].*𝒑([K₁[j],k₂[1]])...))
        for i in 1:(N₂-1)
            a₁,a₂,a₃,a₄=BézPts(t->𝒑([K₁[j],t]),k₂[i],k₂[i+1])
            curve(step*Point([1,-1].*a₂...),step*Point([1,-1].*a₃...),step*Point([1,-1].*a₄...))
        end
    end
    strokepath()
    finish()
    ChangeUnit!(name,"pt","mm")
end

end
