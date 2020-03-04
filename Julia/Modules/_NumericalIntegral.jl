using FastGaussQuadrature

function GaussianQuadrature(f,D₁,D₂;nip=NIP)
    nodes, weights = gausslegendre(nip)
    return sum(
    (weights*weights').*
    [f([x,y]) for
            x ∈ (width(D₁)*nodes.+sum(extrema(D₁)))/2,
            y ∈ (width(D₂)*nodes.+sum(extrema(D₂)))/2
    ])*width(D₁)*width(D₂)/4
end

function GaussianQuadrature(f,D;nip=NIP)
    nodes, weights = gausslegendre(nip)
    return sum(
    weights.*
    [f(x) for x ∈ (width(D)*nodes.+sum(extrema(D)))/2
    ])*width(D)
end
