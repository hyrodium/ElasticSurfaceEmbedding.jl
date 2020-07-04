module ParametricDraw

using IntervalSets
using Colors
using Luxor
using Printf
using ExportNURBS

export SvgCurve, ColorBar

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

function SvgCurve(𝒑s::Array{T,1},I::ClosedInterval;filename="BCA.svg",up=5,down=-5,right=5,left=-5,thickness=1,mesh=50,unitlength=(100,"pt")) where T<:Any
    k=collect(range(endpoints(I)...,length=mesh+1))
    n=length(k)-1
    step, unit=(unitlength[1],unitlength[2])
    Drawing((right-left)*step,(up-down)*step,filename)

    Luxor.origin(-left*step,up*step)
    background("white")

    setline(thickness)
    sethue("cyan")

    for 𝒑 in 𝒑s
        BézPth=BezierPath([BezierPathSegment(map(p->ExportNURBS.LxrPt(p,step),ExportNURBS.BézPts(𝒑,k[i],k[i+1]))...) for i ∈ 1:n])
        drawbezierpath(BézPth, :stroke)
    end

    finish()
    ChangeUnit(filename,"pt",unit)
    return nothing
end

function ColorBar(; max=1.234,filename="ColorBar.png",width=100)
    up=4
    down=-4
    right=4.6
    right=6.2
    left=-2
    Length=3.5
    FontSize=1
    unit=width/(right-left)
    Thickness=unit/10

    Drawing(round(width),round((up-down)*unit),filename)
    Luxor.origin(-left*unit,up*unit)
    setblend(blend(Point(0, -Length*unit), Point(0, Length*unit), "red", "cyan"))
    box(ExportNURBS.LxrPt([-0.9,0],unit), 1.8*unit, 7*unit, :fill)
    sethue("Black")
    fontface("Cica")
    fontsize(unit*FontSize)
    setline(Thickness)
    setlinecap("round")
    text(" "*@sprintf("%.6f",max),ExportNURBS.LxrPt([1.5,Length-0.32*FontSize],unit))
    text(" "*@sprintf("%.6f",0),ExportNURBS.LxrPt([1.5,-0.32*FontSize],unit))
    text("-"*@sprintf("%.6f",max),ExportNURBS.LxrPt([1.5,-Length-0.32*FontSize],unit))
    line(ExportNURBS.LxrPt([0.5,0],unit),ExportNURBS.LxrPt([1.2,0],unit),:stroke)
    line(ExportNURBS.LxrPt([0.5,-Length],unit),ExportNURBS.LxrPt([1.2,-Length],unit),:stroke)
    line(ExportNURBS.LxrPt([0.5,Length],unit),ExportNURBS.LxrPt([1.2,Length],unit),:stroke)

    finish()
    return nothing
end

end
