# The Package ColorBlendModes might work better than the code below, I guess.
# https://github.com/kimikage/ColorBlendModes.jl

import Images
import Images.load
import Images.save
import Images.RGB
import Images.RGBA
import Images.AbstractRGB
import Images.TransparentRGB
import Images.weighted_color_mean
import Images.OffsetArray
import Statistics.mean
using Colors
using Luxor
using Printf
using ExportNURBS

# Images related
function Base.:/(c1::AbstractRGB, c2::Union{AbstractRGB,TransparentRGB})
    return c1
end
function Base.:/(c1::TransparentRGB, c2::TransparentRGB)
    rgb1, rgb2 = RGB(c1), RGB(c2)
    alpha1, alpha2 = c1.alpha, c2.alpha
    # rgb = (alpha1*rgb1 + (1-alpha1)*alpha2*rgb2) / (alpha1 + alpha2 - alpha1*alpha2)
    rgb = weighted_color_mean(alpha1 / (alpha1 + alpha2 - alpha1 * alpha2), rgb1, rgb2)
    alpha = 1 - (1 - alpha1) * (1 - alpha2)
    return typeof(c1)(RGBA(rgb, alpha))
end
function Base.:/(c1::TransparentRGB, c2::AbstractRGB)
    rgb1, rgb2 = RGB(c1), RGB(c2)
    alpha1 = c1.alpha
    # rgb = alpha1*rgb1 + (1-alpha1)*rgb2
    rgb = weighted_color_mean(alpha1, rgb1, rgb2)
    return typeof(c1)(rgb)
end


# Luxor related
function ChangeUnit(filename, before, after)
    ss = open(filename) do file
        strn = read(file, String)
        replace(strn, Regex("(\\d+)" * before) => SubstitutionString("\\1" * after))
    end
    open(filename, "w") do file
        write(file, ss)
    end
    return ss
end
function SvgCurve(
    𝒑s::Array{T,1},
    I::ClosedInterval;
    filename = "BCA.svg",
    up = 5,
    down = -5,
    right = 5,
    left = -5,
    thickness = 1,
    mesh = 50,
    unitlength = (100, "pt"),
) where {T<:Any}
    k = collect(range(endpoints(I)..., length = mesh + 1))
    n = length(k) - 1
    step, unit = (unitlength[1], unitlength[2])
    Drawing((right - left) * step, (up - down) * step, filename)

    Luxor.origin(-left * step, up * step)
    background("white")

    setline(thickness)
    sethue("cyan")

    for 𝒑 in 𝒑s
        BézPth = BezierPath([BezierPathSegment(map(p -> ExportNURBS.LxrPt(p, step), ExportNURBS.BézPts(𝒑, k[i], k[i+1]))...) for i in 1:n])
        drawbezierpath(BézPth, :stroke)
    end

    finish()
    ChangeUnit(filename, "pt", unit)
    return nothing
end
function ColorBar(; max = 1.234, filename = "ColorBar.png", width = 100)
    up = 4
    down = -4
    right = 4.6
    right = 6.2
    left = -2
    Length = 3.5
    FontSize = 1
    unit = width / (right - left)
    Thickness = unit / 10

    Drawing(round(width), round((up - down) * unit), filename)
    Luxor.origin(-left * unit, up * unit)
    setblend(blend(Point(0, -Length * unit), Point(0, Length * unit), "red", "cyan"))
    box(ExportNURBS.LxrPt([-0.9, 0], unit), 1.8 * unit, 7 * unit, :fill)
    sethue("Black")
    fontface("Cica")
    fontsize(unit * FontSize)
    setline(Thickness)
    setlinecap("round")
    text(" " * @sprintf("%.6f", max), ExportNURBS.LxrPt([1.5, Length - 0.32 * FontSize], unit))
    text(" " * @sprintf("%.6f", 0), ExportNURBS.LxrPt([1.5, -0.32 * FontSize], unit))
    text("-" * @sprintf("%.6f", max), ExportNURBS.LxrPt([1.5, -Length - 0.32 * FontSize], unit))
    line(ExportNURBS.LxrPt([0.5, 0], unit), ExportNURBS.LxrPt([1.2, 0], unit), :stroke)
    line(ExportNURBS.LxrPt([0.5, -Length], unit), ExportNURBS.LxrPt([1.2, -Length], unit), :stroke)
    line(ExportNURBS.LxrPt([0.5, Length], unit), ExportNURBS.LxrPt([1.2, Length], unit), :stroke)

    finish()
    return nothing
end
