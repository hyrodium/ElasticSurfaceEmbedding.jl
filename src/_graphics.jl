# Luxor related
function _changeunit(path_svg, units::Pair{String,String})
    old_unit, new_unit = units
    acceptable_units = ("px", "in", "pt", "pc", "cm", "mm")
    if !(new_unit in acceptable_units)
        error("The unit $(new_unit) is not supported in SVG format.")
    end
    script = read(path_svg, String)
    lines = split(script, "\n")
    lines[2] = replace(lines[2], "$(old_unit)\"" => "$(new_unit)\"")
    widthequal = match(r"(width=\"\d+)\"", lines[2])
    if !isnothing(widthequal)
        lines[2] = replace(lines[2], widthequal.captures[1] => widthequal.captures[1]*"mm")
    end
    heightequal = match(r"(height=\"\d+)\"", lines[2])
    if !isnothing(heightequal)
        lines[2] = replace(lines[2], heightequal.captures[1] => heightequal.captures[1]*"mm")
    end
    write(path_svg, join(lines, "\n"))
end

function _colorbar(; max = 1.000, filename = "ColorBar.png", width = 100)
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
    setblend(Luxor.blend(Point(0, -Length * unit), Point(0, Length * unit), "red", "cyan"))
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
    line(
        BasicBSplineExporter._luxor_pt([0.5, -Length], unit),
        BasicBSplineExporter._luxor_pt([1.2, -Length], unit),
        :stroke,
    )
    line(
        BasicBSplineExporter._luxor_pt([0.5, Length], unit),
        BasicBSplineExporter._luxor_pt([1.2, Length], unit),
        :stroke,
    )

    finish()
end
