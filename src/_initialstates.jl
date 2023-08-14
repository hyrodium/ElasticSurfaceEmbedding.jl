"""
    initial_state(D; n₁ = 15)

Compute the initial state, by solving a ODE of center curve.
"""
function initial_state(D::Tuple{ClosedInterval{<:Real}, ClosedInterval{<:Real}}; n₁ = 15)
    D₁, D₂ = D
    M = _initialize(D, n₁)
    comment = "Initial state - domain: " * repr([endpoints(D₁)...]) * "×" * repr([endpoints(D₂)...])
    info = Dict(["type" => "initial"])

    step = Step(M, comment, info)
    steptree = StepTree()
    addstep!(steptree, step, 0)
end

"""
    initial_state!(steptree, D; n₁ = 15)

Compute the initial state, by solving a ODE of center curve.
"""
function initial_state!(steptree, D::Tuple{ClosedInterval{<:Real}, ClosedInterval{<:Real}}; n₁ = 15)
    D₁, D₂ = D
    M = _initialize(D, n₁)
    comment = "Initial state - domain: " * repr([endpoints(D₁)...]) * "×" * repr([endpoints(D₂)...])
    info = Dict(["type" => "initial"])

    step = Step(M, comment, info)
    addstep!(steptree, step, 0)
end

# Coefficient matrix of the center-curve ODE
A(t, D₂) = @SMatrix [
    ṡ₍₀₎(t, D₂)/s₍₀₎(t, D₂) -𝜅₍₀₎(t, D₂)*s₍₀₎(t, D₂)
    𝜅₍₀₎(t, D₂)*s₍₀₎(t, D₂) ṡ₍₀₎(t, D₂)/s₍₀₎(t, D₂)
]

function _initialize(D::Tuple{ClosedInterval{<:Real}, ClosedInterval{<:Real}}, n₁)
    D₁, D₂ = D

    # Definitions for the center curve
    # 1e-14 is ad-hoc number to avoid non-smooth singularity on the boundary.
    t₋ = minimum(D₁) + 1e-14
    t₊ = maximum(D₁) - 1e-14

    # Number of divisions for ODE
    N = 100

    # Initial condition for ODE
    𝒄̇₀ = SVector(s₍₀₎(t₋, D₂), 0.0)

    # Solve ODE 𝒄̈₍ₛ₎(t) = A(t)𝒄̇₍ₛ₎(t) with Runge-Kutta method (and interpolation)
    Δt = (t₊ - t₋) / N
    ts = range(t₋, stop = t₊, length = N + 1)
    𝒄̇₍ₛ₎s = zeros(SVector{2,Float64}, N + 1)
    𝒄̇₍ₛ₎s[1] = 𝒄̇₀
    for i in 1:N
        t = ts[i]
        𝒄̇ = 𝒄̇₍ₛ₎s[i]

        k1 = A(t, D₂) * 𝒄̇
        k2 = A(t + Δt / 2, D₂) * (𝒄̇ + k1 * Δt / 2)
        k3 = A(t + Δt / 2, D₂) * (𝒄̇ + k2 * Δt / 2)
        k4 = A(t + Δt, D₂) * (𝒄̇ + k3 * Δt)

        Δ𝒄̇₀ = Δt * (k1 + 2k2 + 2k3 + k4) / 6
        𝒄̇₍ₛ₎s[i+1] = 𝒄̇ + Δ𝒄̇₀
    end
    𝒄̇₍ₛ₎ = _interpolate2(ts, 𝒄̇₍ₛ₎s, A(t₋, D₂)*𝒄̇₀)

    # Integrate 𝒄̇₍ₛ₎ and obtain the center-curve 𝒄₍ₛ₎
    𝒄₍ₛ₎(t) = unbounded_mapping(integrate(𝒄̇₍ₛ₎), t)

    # Construct initial state M₍ₛ₎
    𝒒₍ₛ₎₁(t) = unbounded_mapping(𝒄̇₍ₛ₎, t)
    𝒒₍ₛ₎₂(t) = (@SMatrix [g₍₀₎₁₂(t, D₂) -𝝊₍₀₎(t, D₂); 𝝊₍₀₎(t, D₂) g₍₀₎₁₂(t, D₂)]) * 𝒒₍ₛ₎₁(t) / g₍₀₎₁₁(t, D₂)
    c = (t₋+t₊)/2
    𝒑₍ₛ₎(u¹, u²) = 𝒄₍ₛ₎(u¹) + (u²-c)*𝒒₍ₛ₎₂(u¹)

    p₁ = 3
    p₂ = 1
    k₁ = KnotVector(range(extrema(D₁)..., length = n₁ - p₁ + 1)) + p₁ * KnotVector([extrema(D₁)...])
    k₂ = KnotVector(repeat(collect(extrema(D₂)), inner = 2))
    P₁ = BSplineSpace{p₁}(k₁)
    P₂ = BSplineSpace{p₂}(k₂)

    𝒂 = fittingcontrolpoints(𝒑₍ₛ₎, P₁, P₂)
    M = BSplineManifold(𝒂, (P₁, P₂))
    M′ = refinement(M, (Val(0), Val(1)))
    return _positioning(M′)
end
