using FastGaussQuadrature

"""
Numerical integral in 2-dimension.
"""
function GaussianQuadrature(f::Function, D₁::ClosedInterval, D₂::ClosedInterval; nip = NIP)
    nodes, weights = gausslegendre(nip)
    return sum(
               (weights * weights') .*
               [f([x, y]) for x in (width(D₁) * nodes .+ sum(extrema(D₁))) / 2, y in (width(D₂) * nodes .+ sum(extrema(D₂))) / 2],
           ) *
           width(D₁) *
           width(D₂) / 4
end

"""
Numerical integral in 1-dimension.
"""
function GaussianQuadrature(f::Function, D::ClosedInterval; nip = NIP)
    nodes, weights = gausslegendre(nip)
    return sum(weights .* [f(x) for x in (width(D) * nodes .+ sum(extrema(D))) / 2]) * width(D)
end
