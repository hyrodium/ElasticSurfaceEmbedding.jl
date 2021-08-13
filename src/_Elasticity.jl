using ForwardDiff

# Strain
E(M, u) = (g₍ₜ₎(M, u) - g₍₀₎(u)) / 2
E₁₁(M::AbstractBSplineManifold, u) = (g₍ₜ₎₁₁(M, u) - g₍₀₎₁₁(u)) / 2
E⁽⁰⁾₁₁(M::AbstractBSplineManifold, u) = E₁₁(M, u) / g₍₀₎₁₁(u)
E₁₁_cont(M::AbstractBSplineManifold, u) = (g₍ₜ₎₁₁_cont(M, u) - g₍₀₎₁₁(u)) / 2
E⁽⁰⁾₁₁_cont(M::AbstractBSplineManifold, u) = E₁₁_cont(M, u) / g₍₀₎₁₁(u)

function Ẽ⁽⁰⁾₁₁(D₂::ClosedInterval, u)
    b = width(D₂) / 2
    c = sum(extrema(D₂)) / 2
    r = (u[2] - c) / b
    return (1 / 2) * K₍₀₎(D₂, u[1]) * B̃(D₂, u[1])^2 * (r^2 - 1 / 3)
end

function Ẽ⁽⁰⁾₁₁(M::AbstractBSplineManifold, u)
    P₁, P₂ = P = bsplinespaces(M)
    p₁, p₂ = degree.(P)
    k₁, k₂ = P₁.knots, P₂.knots
    D₂ = k₂[1+p₂]..k₂[end-p₂]
    return Ẽ⁽⁰⁾₁₁(D₂, u)
end

function ComputeMaximumStrain(; index = 0, mesh = tuple(20 * [MESH...]...))
    M = loadM(index = index)
    𝒂 = controlpoints(M)
    P₁, P₂ = P = bsplinespaces(M)
    p₁, p₂ = p = degree.(P)
    k₁, k₂ = k = knots.(P)
    D₁, D₂ = D = k₁[1+p₁]..k₁[end-p₁], k₂[1+p₂]..k₂[end-p₂]

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

export ShowMaximumStrain
function ShowMaximumStrain(D; index = 0)
    minE, maxE = PredictMaximumStrain(D)

    println("Predicted: (min: ", minE, ", max: ", maxE, ")")

    if isTheShapeComputed()
        minE, maxE = ComputeMaximumStrain(index = index)
        println("Computed: (min: ", minE, ", max: ", maxE, ")")
    end

    return
end


# Elastic Modulus
function C(i, j, k, l, g⁻)
    𝝀 * g⁻[i, j] * g⁻[k, l] + 𝝁 * (g⁻[i, k] * g⁻[j, l] + g⁻[i, l] * g⁻[j, k])
end
