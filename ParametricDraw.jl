module ParametricDraw

using IntervalSets
using Colors
using Luxor
using Printf

export SvgCurve, ParametricColor, ColorBar

function BézPts(𝒑,a,b) # Bézier曲線の制御点
    𝒑(a),
    3*(𝒑(2*a/3+b/3)-(8*𝒑(a)+𝒑(b))/27)-3*(𝒑(a/3+2*b/3)-(𝒑(a)+8*𝒑(b))/27)/2,
    -3*(𝒑(2*a/3+b/3)-(8*𝒑(a)+𝒑(b))/27)/2+3*(𝒑(a/3+2*b/3)-(𝒑(a)+8*𝒑(b))/27),
    𝒑(b)
end

function LxrPt(p::Array{T,1},step) where T<:Real
    Point(step*[1,-1].*p...)
end

function ChangeUnit(filename,before,after)
    ss=open(filename) do file
        strn=read(file, String)
        replace(strn,Regex("(\\d+)"*before)=>SubstitutionString("\\1"*after))
    end
    open(filename,"w") do file
        write(file,ss)
    end
    return ss
end

function SvgCurve(𝒑,k::Array{T,1};filename="BCA.svg",up=5,down=-5,right=5,left=-5,thickness=1,unitlength=(100,"pt")) where T<:Real
    n=length(k)-1
    step, unit=(unitlength[1],unitlength[2])
    Drawing((right-left)*step,(up-down)*step,filename)

    Luxor.origin(-left*step,up*step)
    background("white")

    BézPth=BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(𝒑,k[i],k[i+1]))...) for i ∈ 1:n])

    setline(thickness)
    sethue("red")
    drawbezierpath(BézPth, :stroke)
    finish()
    ChangeUnit(filename,"pt",unit)
    return nothing
end

function SvgCurve(𝒑,I::ClosedInterval;filename="BCA.svg",up=5,down=-5,right=5,left=-5,thickness=1,mesh=50,unitlength=(100,"pt"))
    k=collect(range(endpoints(I)...,length=mesh))
    n=length(k)-1
    step, unit=(unitlength[1],unitlength[2])
    Drawing((right-left)*step,(up-down)*step,filename)

    Luxor.origin(-left*step,up*step)
    background("white")

    setline(thickness)
    sethue("red")

    BézPth=BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(𝒑,k[i],k[i+1]))...) for i ∈ 1:n])
    drawbezierpath(BézPth, :stroke)

    finish()
    ChangeUnit(filename,"pt",unit)
    return nothing
end

function SvgCurve(𝒑s::Array{T,1},I::ClosedInterval;filename="BCA.svg",up=5,down=-5,right=5,left=-5,thickness=1,mesh=50,unitlength=(100,"pt")) where T<:Any
    k=collect(range(endpoints(I)...,length=mesh))
    n=length(k)-1
    step, unit=(unitlength[1],unitlength[2])
    Drawing((right-left)*step,(up-down)*step,filename)

    Luxor.origin(-left*step,up*step)
    background("white")

    setline(thickness)
    sethue("cyan")

    for 𝒑 ∈ 𝒑s
        BézPth=BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(𝒑,k[i],k[i+1]))...) for i ∈ 1:n])
        drawbezierpath(BézPth, :stroke)
    end

    finish()
    ChangeUnit(filename,"pt",unit)
    return nothing
end

function SvgSurface(𝒑,k,n;filename="BSA.svg",up=5,down=-5,right=5,left=-5,step=50)
    k₁,k₂=k
    n₁,n₂=n
    K₁,K₂=range(minimum(k₁),maximum(k₁),length=n₁+1), range(minimum(k₂),maximum(k₂),length=n₂+1)
    N₁,N₂=length(k₁),length(k₂)

    Drawing(step*(right-left),step*(up-down),filename)
    Luxor.origin(-step*left,step*up)
    background("white")
    sethue("blue")

    for j ∈ 1:(n₂+1)
        move(step*Point([1,-1].*𝒑([k₁[1],K₂[j]])...))
        for i ∈ 1:(N₁-1)
            a₁,a₂,a₃,a₄=BézPts(t->𝒑([t,K₂[j]]),k₁[i],k₁[i+1])
            curve(step*Point([1,-1].*a₂...),step*Point([1,-1].*a₃...),step*Point([1,-1].*a₄...))
        end
    end
    for j ∈ 1:(n₁+1)
        move(step*Point([1,-1].*𝒑([K₁[j],k₂[1]])...))
        for i ∈ 1:(N₂-1)
            a₁,a₂,a₃,a₄=BézPts(t->𝒑([K₁[j],t]),k₂[i],k₂[i+1])
            curve(step*Point([1,-1].*a₂...),step*Point([1,-1].*a₃...),step*Point([1,-1].*a₄...))
        end
    end
    strokepath()
    finish()
    ChangeUnit!(filename,"pt","mm")
end

function ParametricColor(𝒑,D;rgb=(u->[0.5,0.5,0.5]),filename="ParametricColor.png",up=5,down=-5,right=5,left=-5,mesh=(10,10),unit=100)
    Drawing((right-left)*unit,(up-down)*unit,filename)
    Luxor.origin(-left*unit,up*unit)

    k₁=range(leftendpoint(D[1]),stop=rightendpoint(D[1]),length=mesh[1]+1)
    k₂=range(leftendpoint(D[2]),stop=rightendpoint(D[2]),length=mesh[2]+1)

    for I₁ ∈ 1:mesh[1], I₂ ∈ 1:mesh[2]
        BézPth=BezierPath([
                BezierPathSegment(map(p->LxrPt(p,unit),BézPts(t->𝒑([t,k₂[I₂]]),k₁[I₁],k₁[I₁+1]))...),
                BezierPathSegment(map(p->LxrPt(p,unit),BézPts(t->𝒑([k₁[I₁+1],t]),k₂[I₂],k₂[I₂+1]))...),
                BezierPathSegment(map(p->LxrPt(p,unit),BézPts(t->𝒑([t,k₂[I₂+1]]),k₁[I₁+1],k₁[I₁]))...),
                BezierPathSegment(map(p->LxrPt(p,unit),BézPts(t->𝒑([k₁[I₁],t]),k₂[I₂+1],k₂[I₂]))...)])
        mesh1 = Luxor.mesh(BézPth, [
            Colors.RGB(rgb([k₁[I₁], k₂[I₂]])...), # (k₁[I₁], k₂[I₂])
            Colors.RGB(rgb([k₁[I₁+1], k₂[I₂]])...), # (k₁[I₁+1], k₂[I₂])
            Colors.RGB(rgb([k₁[I₁+1], k₂[I₂+1]])...), # (k₁[I₁+1], k₂[I₂+1])
            Colors.RGB(rgb([k₁[I₁], k₂[I₂+1]])...)  # (k₁[I₁], k₂[I₂+1])
            ])
        setmesh(mesh1)
        box(LxrPt([right+left,up+down]/2,unit), (right-left)*unit,(up-down)*unit,:fill)
    end
    finish()
    return nothing
end

function ColorBar(;max=1.234,filename="ColorBar.png",width=100)
    up=4
    down=-4
    right=4.6
    left=-2
    Length=3.5
    FontSize=1
    unit=width/(right-left)
    Thickness=unit/10

    Drawing(round(width),round((up-down)*unit),filename)
    Luxor.origin(-left*unit,up*unit)
    setblend(blend(Point(0, -Length*unit), Point(0, Length*unit), "red", "cyan"))
    box(LxrPt([-0.9,0],unit), 1.8*unit, 7*unit, :fill)
    sethue("Black")
    fontface("Cica")
    fontsize(unit*FontSize)
    setline(Thickness)
    setlinecap("round")
    text(" "*@sprintf("%.3f",max),LxrPt([1.5,Length-0.32*FontSize],unit))
    text(" "*@sprintf("%.3f",0),LxrPt([1.5,-0.32*FontSize],unit))
    text("-"*@sprintf("%.3f",max),LxrPt([1.5,-Length-0.32*FontSize],unit))
    line(LxrPt([0.5,0],unit),LxrPt([1.2,0],unit),:stroke)
    line(LxrPt([0.5,-Length],unit),LxrPt([1.2,-Length],unit),:stroke)
    line(LxrPt([0.5,Length],unit),LxrPt([1.2,Length],unit),:stroke)

    finish()
    return nothing
end

end
