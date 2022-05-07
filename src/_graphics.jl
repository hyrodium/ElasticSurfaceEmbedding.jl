# The Package ColorBlendModes might work better than the code below, I guess.
# https://github.com/kimikage/ColorBlendModes.jl

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
function _changeunit(path_svg, units::Pair{String,String})
    old_unit, new_unit = units
    acceptable_units = ["px", "in", "pt", "pc", "cm", "mm"]
    if !(new_unit in acceptable_units)
        error("The unit $(new_unit) is not support in SVG format.")
    end
    script = read(path_svg, String)
    lines = split(script, "\n")
    lines[2] = replace(lines[2],"$(old_unit)\""=>"$(new_unit)\"")
    write(path_svg, join(lines,"\n"))
end

function _svgcurve(𝒑s::Array{T,1}, I::ClosedInterval; filename, up=5, down=-5, right=5, left=-5, thickness=1, mesh=50, unitlength=100) where {T<:Any}
    k = collect(range(endpoints(I)..., length = mesh + 1))
    n = length(k) - 1
    step = unitlength
    Drawing((right - left) * step, (up - down) * step, filename)

    Luxor.origin(-left * step, up * step)
    background("white")

    setline(thickness)
    sethue("cyan")

    for 𝒑 in 𝒑s
        BézPth = BezierPath([BezierPathSegment(map(p -> BasicBSplineExporter._luxor_pt(p, step), BasicBSplineExporter._bezier(𝒑, k[i], k[i+1]))...) for i in 1:n])
        drawbezierpath(BézPth, :stroke)
    end

    finish()
    return
end
function _colorbar(; max=1.000, filename="ColorBar.png", width=100)
    up = 4
    down = -4
    right = 4.6
    right = 6.2
    left = -2
    Length = 3.5
    FontSize = 0.85
    unit = width / (right - left)
    Thickness = unit / 10

    Drawing(round(width), round((up - down) * unit), filename)
    Luxor.origin(-left * unit, up * unit)
    setblend(blend(Point(0, -Length * unit), Point(0, Length * unit), "red", "cyan"))
    box(BasicBSplineExporter._luxor_pt([-0.9, 0], unit), 1.8 * unit, 7 * unit, :fill)
    sethue("Black")
    fontface("JuliaMono")
    fontsize(unit * FontSize)
    setline(Thickness)
    setlinecap("round")
    text(" " * @sprintf("%.6f", max), BasicBSplineExporter._luxor_pt([1.4, Length - 0.28 * FontSize], unit))
    text(" " * @sprintf("%.6f", 0), BasicBSplineExporter._luxor_pt([1.4, -0.28 * FontSize], unit))
    text("-" * @sprintf("%.6f", max), BasicBSplineExporter._luxor_pt([1.4, -Length - 0.28 * FontSize], unit))
    line(BasicBSplineExporter._luxor_pt([0.5, 0], unit), BasicBSplineExporter._luxor_pt([1.2, 0], unit), :stroke)
    line(BasicBSplineExporter._luxor_pt([0.5, -Length], unit), BasicBSplineExporter._luxor_pt([1.2, -Length], unit), :stroke)
    line(BasicBSplineExporter._luxor_pt([0.5, Length], unit), BasicBSplineExporter._luxor_pt([1.2, Length], unit), :stroke)

    finish()
    return
end
