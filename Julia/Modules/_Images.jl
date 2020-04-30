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

function Base.:/(c1::AbstractRGB, c2::Union{AbstractRGB, TransparentRGB})
    return c1
end
function Base.:/(c1::TransparentRGB, c2::TransparentRGB)
    rgb1, rgb2 = RGB(c1), RGB(c2)
    alpha1, alpha2 = c1.alpha, c2.alpha
    # rgb = (alpha1*rgb1 + (1-alpha1)*alpha2*rgb2) / (alpha1 + alpha2 - alpha1*alpha2)
    rgb = weighted_color_mean(alpha1/(alpha1 + alpha2 - alpha1*alpha2), rgb1, rgb2)
    alpha = 1 - (1-alpha1) * (1-alpha2)
    return typeof(c1)(RGBA(rgb, alpha))
end
function Base.:/(c1::TransparentRGB, c2::AbstractRGB)
    rgb1, rgb2 = RGB(c1), RGB(c2)
    alpha1 = c1.alpha
    # rgb = alpha1*rgb1 + (1-alpha1)*rgb2
    rgb = weighted_color_mean(alpha1, rgb1, rgb2)
    return typeof(c1)(rgb)
end
