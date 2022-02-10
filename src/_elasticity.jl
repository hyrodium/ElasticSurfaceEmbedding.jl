# Strain related functions
E(M,u¹,u²) = (g₍ₜ₎(M,u¹,u²) - g₍₀₎(u¹,u²)) / 2
E₁₁(M::BSplineManifold{2},u¹,u²) = (g₍ₜ₎₁₁(M,u¹,u²) - g₍₀₎₁₁(u¹,u²)) / 2
E⁽⁰⁾₁₁(M::BSplineManifold{2},u¹,u²) = E₁₁(M,u¹,u²) / g₍₀₎₁₁(u¹,u²)

function Ẽ⁽⁰⁾₁₁(D₂::ClosedInterval,u¹,u²)
    # Breadth of the strip-like shape
    b = width(D₂) / 2
    # Center coordinate of u²
    c = sum(extrema(D₂)) / 2
    # Normalized coordinate of u²
    r = (u² - c) / b
    # Compute the predicted strain with the Strain Approximation Theorem
    return (1/2) * K₍₀₎(u¹,D₂) * B̃(u¹, D₂)^2 * (r^2 - 1 / 3)
end

function Ẽ⁽⁰⁾₁₁(M::BSplineManifold{2},u¹,u²)
    _, P₂ = bsplinespaces(M)
    p₂ = degree(P₂)
    k₂ = knotvector(P₂)
    D₂ = k₂[1+p₂]..k₂[end-p₂]
    return Ẽ⁽⁰⁾₁₁(D₂,u¹,u²)
end

function _compute_minmax_strain(; index=0, mesh=tuple(20 * [MESH...]...))
    M = loadM(index = index)
    P = bsplinespaces(M)
    p₁, p₂ = degree.(P)
    k₁, k₂ = knotvector.(P)
    D₁, D₂ = k₁[1+p₁]..k₁[end-p₁], k₂[1+p₂]..k₂[end-p₂]

    κ₁ = range(leftendpoint(D₁) + 0.0001, stop = rightendpoint(D₁) - 0.0001, length = mesh[1] + 1)
    κ₂ = range(leftendpoint(D₂) + 0.0001, stop = rightendpoint(D₂) - 0.0001, length = mesh[2] + 1)

    E = [E⁽⁰⁾₁₁(M, u₁, u₂) for u₁ in κ₁, u₂ in κ₂]

    return (minimum(E), maximum(E))
end

function _predict_minmax_strain(D; mesh = tuple(20 * [MESH...]...))
    D₁, D₂ = D

    κ₁ = range(leftendpoint(D₁), stop = rightendpoint(D₁), length = mesh[1] + 1)
    κ₂ = range(leftendpoint(D₂), stop = rightendpoint(D₂), length = mesh[2] + 1)

    E = [Ẽ⁽⁰⁾₁₁(D₂, u₁, u₂) for u₁ in κ₁, u₂ in κ₂]

    return (minimum(E), maximum(E))
end

"""
    show_strain(D; index=0)

Show the predicted maximum strain and, if possible, also the computed strain with the given index.
"""
function show_strain(D; index=0)
    minE, maxE = _predict_minmax_strain(D)

    D₁, D₂ = D
    msg = "Strain - domain: " * repr([endpoints(D₁)...]) * "×" * repr([endpoints(D₂)...]) * "\n"
    msg *= "Predicted: (min: $(minE), max: $(maxE))\n"

    if isTheShapeComputed()
        minE, maxE = _compute_minmax_strain(index = index)
        msg *= "Computed: (min: $(minE), max: $(maxE))\n"
    end
    
    @info msg

    return
end


# Elastic Modulus
function C(i, j, k, l, g⁻)
    𝝀 * g⁻[i, j] * g⁻[k, l] + 𝝁 * (g⁻[i, k] * g⁻[j, l] + g⁻[i, l] * g⁻[j, k])
end
