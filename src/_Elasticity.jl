# Strain related functions
E(M, u) = (g₍ₜ₎(M, u) - g₍₀₎(u)) / 2
E₁₁(M::AbstractBSplineManifold, u) = (g₍ₜ₎₁₁(M, u) - g₍₀₎₁₁(u)) / 2
E⁽⁰⁾₁₁(M::AbstractBSplineManifold, u) = E₁₁(M, u) / g₍₀₎₁₁(u)
E₁₁_cont(M::AbstractBSplineManifold, u) = (g₍ₜ₎₁₁_cont(M, u) - g₍₀₎₁₁(u)) / 2
E⁽⁰⁾₁₁_cont(M::AbstractBSplineManifold, u) = E₁₁_cont(M, u) / g₍₀₎₁₁(u)

function Ẽ⁽⁰⁾₁₁(D₂::ClosedInterval, u)
    # Breadth of the strip-like shape
    b = width(D₂) / 2
    # Center coordinate of u²
    c = sum(extrema(D₂)) / 2
    # Normalized coordinate of u²
    r = (u[2] - c) / b
    # Compute the predicted strain with the Strain Approximation Theorem
    return (1 / 2) * K₍₀₎(D₂, u[1]) * B̃(D₂, u[1])^2 * (r^2 - 1 / 3)
end

function Ẽ⁽⁰⁾₁₁(M::AbstractBSplineManifold, u)
    _, P₂ = bsplinespaces(M)
    p₂ = degree(P₂)
    k₂ = P₂.knots
    D₂ = k₂[1+p₂]..k₂[end-p₂]
    return Ẽ⁽⁰⁾₁₁(D₂, u)
end

function ComputeMaximumStrain(; index = 0, mesh = tuple(20 * [MESH...]...))
    M = loadM(index = index)
    P = bsplinespaces(M)
    p₁, p₂ = degree.(P)
    k₁, k₂ = knots.(P)
    D₁, D₂ = k₁[1+p₁]..k₁[end-p₁], k₂[1+p₂]..k₂[end-p₂]

    κ₁ = range(leftendpoint(D₁) + 0.0001, stop = rightendpoint(D₁) - 0.0001, length = mesh[1] + 1)
    κ₂ = range(leftendpoint(D₂) + 0.0001, stop = rightendpoint(D₂) - 0.0001, length = mesh[2] + 1)

    E = [E⁽⁰⁾₁₁(M, [u₁, u₂]) for u₁ in κ₁, u₂ in κ₂]

    return (minimum(E), maximum(E))
end

function PredictMaximumStrain(D; mesh = tuple(20 * [MESH...]...))
    D₁, D₂ = D

    κ₁ = range(leftendpoint(D₁), stop = rightendpoint(D₁), length = mesh[1] + 1)
    κ₂ = range(leftendpoint(D₂), stop = rightendpoint(D₂), length = mesh[2] + 1)

    E = [Ẽ⁽⁰⁾₁₁(D₂, [u₁, u₂]) for u₁ in κ₁, u₂ in κ₂]

    return (minimum(E), maximum(E))
end

"""
    print_strain(D; index = 0)

Show the predicted maximum strain and, if possible, also the computed strain with the given index.
"""
function print_strain(D; index = 0)
    minE, maxE = PredictMaximumStrain(D)

    D₁, D₂ = D
    msg = "Strain - domain: " * repr([endpoints(D₁)...]) * "×" * repr([endpoints(D₂)...]) * "\n"
    msg *= "Predicted: (min: $(minE), max: $(maxE))\n"

    if isTheShapeComputed()
        minE, maxE = ComputeMaximumStrain(index = index)
        msg *= "Computed: (min: $(minE), max: $(maxE))\n"
    end
    
    @info msg

    return
end


# Elastic Modulus
function C(i, j, k, l, g⁻)
    𝝀 * g⁻[i, j] * g⁻[k, l] + 𝝁 * (g⁻[i, k] * g⁻[j, l] + g⁻[i, l] * g⁻[j, k])
end
