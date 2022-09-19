"""
    initial_state(D; n₁ = 15)

Compute the initial state, by solving a ODE of center curve.
"""
function initial_state(D; n₁ = 15)
    D₁, D₂ = D
    M = _initialize(D, n₁)
    comment = "Initial state - domain: " * repr([endpoints(D₁)...]) * "×" * repr([endpoints(D₂)...])
    info = Dict(["type"=>"initial"])

    step = Step(M, comment, info)
    allsteps = AllSteps()
    addstep!(allsteps, step, 0)
end

"""
    initial_state!(allsteps, D; n₁ = 15)

Compute the initial state, by solving a ODE of center curve.
"""
function initial_state!(allsteps, D; n₁ = 15)
    D₁, D₂ = D
    M = _initialize(D, n₁)
    comment = "Initial state - domain: " * repr([endpoints(D₁)...]) * "×" * repr([endpoints(D₂)...])
    info = Dict(["type"=>"initial"])

    step = Step(M, comment, info)
    addstep!(allsteps, step, 0)
end

function _initialize(D, n₁)
    D₁, D₂ = D

    # Definition for center curve
    t₋, t₊ = extrema(D₁)
    p₁ = 3
    p₂ = 1
    k₁ = KnotVector(range(t₋, t₊, length=n₁-p₁+1)) + p₁ * KnotVector([t₋, t₊])
    k₂ = KnotVector(repeat(collect(extrema(D₂)), inner = 2))
    P₁ = BSplineSpace{p₁}(k₁)
    P₂ = BSplineSpace{p₂}(k₂)

    # Number of divisions for ODE
    N = 6n₁

    # Solve 𝒄̈(t) = A(t)𝒄̇(t) with Runge-Kutta method
    A(t) = @SMatrix[
        ṡ₍₀₎(t, D₂)/s₍₀₎(t, D₂)  -𝜅₍₀₎(t, D₂)*s₍₀₎(t, D₂)
        𝜅₍₀₎(t, D₂)*s₍₀₎(t, D₂)  ṡ₍₀₎(t, D₂)/s₍₀₎(t, D₂)
    ]

    # Initial condition
    𝒄̇₀ = SVector(1.0, 0.0) * s₍₀₎(minimum(D₁), D₂)

    Δt = (t₊-t₋)/N
    ts = range(t₋, stop=t₊, length=N+1)
    𝒄̇s = zeros(SVector{2,Float64}, N+1)
    𝒄̇s[1] = 𝒄̇₀
    for i in 1:N
        t = ts[i]
        𝒄̇ = 𝒄̇s[i]

        k1 = A(t)*𝒄̇
        k2 = A(t+Δt/2)*(𝒄̇+k1*Δt/2)
        k3 = A(t+Δt/2)*(𝒄̇+k2*Δt/2)
        k4 = A(t+Δt)*(𝒄̇+k3*Δt)

        Δ𝒄̇₀ = Δt * (k1+2k2+2k3+k4)/6
        𝒄̇s[i+1] = 𝒄̇ + Δ𝒄̇₀
    end

    # Approximate 𝒄̇=𝒒₁ with B-spline curve
    _p₁ = p₁-1
    _k₁ = KnotVector(range(t₋, t₊, length = n₁-_p₁)) + _p₁ * KnotVector([t₋, t₊])
    _P₁ = BSplineSpace{_p₁}(_k₁)
    _n₁ = dim(_P₁)
    _B = [bsplinebasis(_P₁,i,t) for i in 1:_n₁, t in ts]
    _BB = _B * _B'
    _b = _B * 𝒄̇s
    𝒎̇ = inv(_BB)*_b  # control points of 𝒒̃₁

    # Approximate 𝒄 with B-spline curve
    Δk = (t₊-t₋)/(n₁-p₁)
    𝒎 = zeros(SVector{2,Float64},dim(P₁))  # control points of 𝒄̃
    𝒎[1] = zero(SVector{2,Float64})
    𝒎[2] = 𝒎[1] + 𝒎̇[1]*Δk*1/3
    𝒎[3] = 𝒎[2] + 𝒎̇[2]*Δk*2/3
    for i in 3:dim(P₁)-1
        𝒎[i+1] = 𝒎[i] + 𝒎̇[i]*Δk
    end
    𝒎[n₁-1] = 𝒎[n₁-2] + 𝒎̇[n₁-2]*Δk*2/3
    𝒎[n₁] = 𝒎[n₁-1] + 𝒎̇[n₁-1]*Δk*1/3

    # Approximate 𝒒₂ with B-spline curve
    𝒒₂s = [
        (@SMatrix [g₍₀₎₁₂(ts[i],D₂) -𝝊₍₀₎(ts[i],D₂); 𝝊₍₀₎(ts[i],D₂) g₍₀₎₁₂(ts[i],D₂)]) * 𝒄̇s[i] / g₍₀₎₁₁(ts[i],D₂)
        for i in 1:N+1
    ]

    _B = [bsplinebasis(P₁,i,t) for i in 1:n₁, t in ts]
    _BB = _B * _B'
    _b = _B * 𝒒₂s
    𝒓 = inv(_BB)*_b  # control points of 𝒄̃₂

    a1 = 𝒎 - width(D₂) * 𝒓/2
    a2 = 𝒎 + width(D₂) * 𝒓/2
    𝒂 = hcat(a1,a2)

    M = BSplineManifold(𝒂, (P₁, P₂))
    M′ = refinement(M, (Val(0), Val(1)))
    return _positioning(M′)
end
