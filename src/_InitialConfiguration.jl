"""
    initial_configulation(D; n₁ = 15)

Compute the initial configulation, by solving a ODE.
"""
function initial_configulation(D; n₁ = 15)
    parent = 0

    D₁, D₂ = D
    M = _initialize(D, n₁)
    comment = "Initial Configuration - domain: " * repr([endpoints(D₁)...]) * "×" * repr([endpoints(D₂)...])

    _export(M, parent, comment = comment)
end

function _initialize(D, n₁)
    D₁, D₂ = D

    # Definition for center curve
    t₋, t₊ = extrema(D₁)
    p₁ = 3
    p₂ = 1
    k₁ = Knots(range(t₋, t₊, length = n₁-p₁+1)) + p₁ * Knots(t₋, t₊)
    k₂ = Knots(repeat(collect(extrema(D₂)), inner = 2))
    P₁ = FastBSplineSpace(p₁, k₁)
    P₂ = FastBSplineSpace(p₂, k₂)
    n₂ = dim(P₂)


    # Number of divisions for ODE
    N = 128

    # Solve 𝒄̈(t) = A(t)𝒄̇(t) with Runge-Kutta method
    A(t) = [s̈₍₀₎(D₂, t) / ṡ₍₀₎(D₂, t) -𝜅₍₀₎(D₂, t) * ṡ₍₀₎(D₂, t)
    𝜅₍₀₎(D₂, t) * ṡ₍₀₎(D₂, t) s̈₍₀₎(D₂, t) / ṡ₍₀₎(D₂, t)]

    # Initial condition
    𝒄̇₀ = [1.0, 0.0] * ṡ₍₀₎(D₂, minimum(D₁))

    Δt = (t₊-t₋)/N
    ts = range(t₋, stop=t₊, length=N+1)
    𝒄̇s = zeros(N+1,2)
    𝒄̇s[1,:] = 𝒄̇₀
    for i in 1:N
        t = ts[i]
        𝒄̇ = 𝒄̇s[i,:]

        k1 = A(t)*𝒄̇
        k2 = A(t+Δt/2)*(𝒄̇+k1*Δt/2)
        k3 = A(t+Δt/2)*(𝒄̇+k2*Δt/2)
        k4 = A(t+Δt)*(𝒄̇+k3*Δt)

        Δ𝒄̇₀ = Δt * (k1+2k2+2k3+k4)/6
        𝒄̇s[i+1,:] = 𝒄̇ + Δ𝒄̇₀
    end

    # Approximate 𝒄̇=𝒄₁ with B-spline curve
    _p₁ = p₁-1
    _k₁ = Knots(range(t₋, t₊, length = n₁-_p₁)) + _p₁ * Knots(t₋, t₊)
    _P₁ = FastBSplineSpace(_p₁, _k₁)
    _n₁ =  dim(_P₁)
    _B = [bsplinebasis(i, _P₁, t) for i in 1:_n₁, t in ts]
    _BB = _B * _B'
    _b = _B * 𝒄̇s
    𝒎̇ = _BB\_b  # control points of 𝒄̃₁

    # Approximate 𝒄 with B-spline curve
    Δk = (t₊-t₋)/(n₁-p₁)
    𝒎 = zeros(dim(P₁),2)  # control points of 𝒄̃
    𝒎[1,:] = zeros(2)
    𝒎[2,:] = 𝒎[1,:] + 𝒎̇[1,:]*Δk*1/3
    𝒎[3,:] = 𝒎[2,:] + 𝒎̇[2,:]*Δk*2/3
    for i in 3:dim(P₁)-1
        𝒎[i+1,:] = 𝒎[i,:] + 𝒎̇[i,:]*Δk
    end
    𝒎[n₁-1,:] = 𝒎[n₁-2,:] + 𝒎̇[n₁-2,:]*Δk*2/3
    𝒎[n₁,:] = 𝒎[n₁-1,:] + 𝒎̇[n₁-1,:]*Δk*1/3

    # Approximate 𝒄₂ with B-spline curve
    𝒄₂s = [[g₍₀₎₁₂(c(D₂, ts[i])) -𝝊₍₀₎(c(D₂, ts[i])); 𝝊₍₀₎(c(D₂, ts[i])) g₍₀₎₁₂(c(D₂, ts[i]))] * 𝒄̇s[i,:] / g₍₀₎₁₁(c(D₂, ts[i])) for i in 1:N+1]
    𝒄₂s = hcat(𝒄₂s...)'

    _B = [bsplinebasis(i, P₁, t) for i in 1:n₁, t in ts]
    _BB = _B * _B'
    _b = _B * 𝒄₂s
    𝒓 = _BB\_b  # control points of 𝒄̃₂

    a1 = 𝒎 - width(D₂) * 𝒓/2
    a2 = 𝒎 + width(D₂) * 𝒓/2
    𝒂 = [(a1[I₁,i],a2[I₁,i])[I₂] for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:2]

    M = BSplineSurface([P₁, P₂], 𝒂)
    M′ = refinement(M, p₊ = [0, 1])
    return Positioning(M′)
end
