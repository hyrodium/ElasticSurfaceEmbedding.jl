using Dates

function DefaultOrientation(n₁, n₂)
    return ([(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2, 1], [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2, 2], [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2 - 1, 1])
end

function FixThreePoints(n₁, n₂)
    return (
        [1, (n₂ + 1) ÷ 2, 1],
        [1, (n₂ + 1) ÷ 2, 2],
        [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2, 1],
        [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2, 2],
        [n₁, (n₂ + 1) ÷ 2, 1],
        [n₁, (n₂ + 1) ÷ 2, 2],
    )
end

AbbStr(t::Week) = string(t.value) * "w "
AbbStr(t::Day) = string(t.value) * "d "
AbbStr(t::Hour) = string(t.value) * "h "
AbbStr(t::Minute) = string(t.value) * "m "
AbbStr(t::Second) = string(t.value) * "s "
AbbStr(t::Millisecond) = string(t.value) * "ms "
AbbStr(t::Vector{Period}) = *(AbbStr.(t)...)[1:end-1]

function SecondsToString(Δt::Float64)
    prds = Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(floor(1000Δt)))).periods
    return AbbStr(prds)
end

export NewtonMethodIteration
function NewtonMethodIteration(; fixingmethod = :DefaultOrientation, parent::Int = 0, nip = NIP)
    if fixingmethod == :DefaultOrientation
        fixed = DefaultOrientation
    elseif fixingmethod == :FixThreePoints
        fixed = FixThreePoints
    else
        error("No method for $(fixingmethod)")
    end
    parent = Parent(parent)
    M = loadM(index = parent)

    n₁, n₂ = n = dim.(bsplinespaces(M))
    if !isodd(n₁ * n₂)
        error("n₁ and n₂ should be odd numbers")
    end
    M = Positioning(M)
    M, F, Ǧ, Δt = NewtonIteration(M, fixed, nip = nip)
    comment =
        "Newton Iteration - residual norm: " *
        (@sprintf("%.4e", norm(F))) *
        ", Δa norm: " *
        (@sprintf("%.4e", norm(Ǧ))) *
        ", computation time: " *
        SecondsToString(Δt)
    Export(M, parent, comment = comment)
end

function NewtonIteration(M::AbstractBSplineManifold, fixed; nip = NIP)
    𝒂 = controlpoints(M)
    P₁, P₂ = P = bsplinespaces(M)
    p₁, p₂ = p = degree.(P)
    k₁, k₂ = k = knots.(P)
    D₁, D₂ = D = k₁[1+p₁]..k₁[end-p₁], k₂[1+p₂]..k₂[end-p₂]
    n₁, n₂ = n = dim.(P)
    function lineup(I₁::Int, I₂::Int, i::Int)::Int
        return (i - 1) * n₁ * n₂ + (I₂ - 1) * n₁ + (I₁ - 1) + 1
    end

    t₀ = time()
    if distributed
        f = Array{Union{Future,Nothing}}(nothing, n₁, n₂, d)
        for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:d
            f[I₁, I₂, i] = @spawn elm_F(M, I₁, I₂, i, nip = nip)
        end
        h = Array{Union{Future,Nothing}}(nothing, n₁, n₂, d, n₁, n₂, d)
        for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:d, R₁ in 1:n₁, R₂ in 1:n₂, r in 1:d
            if lineup(I₁, I₂, i) ≤ lineup(R₁, R₂, r)
                h[I₁, I₂, i, R₁, R₂, r] = h[R₁, R₂, r, I₁, I₂, i] = @spawn elm_H(M, I₁, I₂, i, R₁, R₂, r, nip = nip)
            end
        end
        F = fetch.(f)
        H = fetch.(h)
    else
        H = [elm_H(M, I₁, I₂, i, R₁, R₂, r, nip = nip) for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:d, R₁ in 1:n₁, R₂ in 1:n₂, r in 1:d]
        F = [elm_F(M, I₁, I₂, i, nip = nip) for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:d]
    end
    t₁ = time()

    𝕟 = n₁ * n₂ * d
    Fixed = sort(collect((i -> lineup(i...)).(fixed(n₁, n₂))))
    Unfixed = deleteat!(collect(1:𝕟), Fixed)

    F = reshape(F, 𝕟)
    H = reshape(H, 𝕟, 𝕟)
    a = aₒ = reshape(𝒂, 𝕟)
    Ȟ = H[Unfixed, Unfixed]
    ǎ = a[Unfixed]
    F̌ = F[Unfixed]
    Ǧ = Ȟ \ F̌
    ǎ = ǎ - Ǧ
    for i in Fixed
        insert!(ǎ, i, aₒ[i])
    end
    𝒂 = reshape(ǎ, n₁, n₂, d)
    M = typeof(M)(P, 𝒂)
    return (M, F, Ǧ, t₁ - t₀)
end

function elm_H(M::AbstractBSplineManifold, I₁, I₂, i, R₁, R₂, r; nip = NIP)
    𝒂 = controlpoints(M)
    P₁, P₂ = P = bsplinespaces(M)
    p₁, p₂ = p = degree.(P)
    k₁, k₂ = k = knots.(P)
    n₁, n₂ = n = dim.(P)

    P₁, P₂ = 𝒫(p₁, k₁), 𝒫(p₂, k₂)

    𝜹 = [1.0 0.0; 0.0 1.0]
    Σ₁ = (maximum([I₁, R₁]):minimum([I₁, R₁])+p₁)
    Σ₂ = (maximum([I₂, R₂]):minimum([I₂, R₂])+p₂)

    if length(Σ₁) == 0 || length(Σ₂) == 0
        return 0.0
    else
        return sum(
            GaussianQuadrature(
                u ->
                    (
                        g = g₍₀₎(u);
                        g⁻ = inv(g);
                        𝝊 = sqrt(det(g));
                        𝑁 = [N′(P₁, P₂, I₁, I₂, i, u) for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:d];
                        Q = [sum(𝒂[I₁, I₂, i] * 𝑁[I₁, I₂, j] for I₁ in 1:n₁, I₂ in 1:n₂) for i in 1:d, j in 1:d];
                        sum(
                            C(p, q, m, n, g⁻) *
                            𝑁[I₁, I₂, p] *
                            (𝜹[i, r] * 𝑁[R₁, R₂, q] * (sum(Q[o, m] * Q[o, n] for o in 1:d) - g[m, n]) + 2 * 𝑁[R₁, R₂, n] * Q[i, q] * Q[r, m])
                            for p in 1:d, q in 1:d, m in 1:d, n in 1:d
                        )
                    ) * 𝝊,
                k₁[ι₁]..k₁[ι₁+1],
                k₂[ι₂]..k₂[ι₂+1],
                nip = nip,
            ) for ι₁ in Σ₁, ι₂ in Σ₂
        )
    end
end

function elm_F(M::AbstractBSplineManifold, I₁, I₂, i; nip = NIP)
    𝒂 = controlpoints(M)
    P₁, P₂ = P = bsplinespaces(M)
    p₁, p₂ = p = degree.(P)
    k₁, k₂ = k = knots.(P)
    n₁, n₂ = n = dim.(P)

    D̂₁ = bsplinesupport(I₁, P₁)
    D̂₂ = bsplinesupport(I₂, P₂)
    Σ₁ = I₁:I₁+p₁
    Σ₂ = I₂:I₂+p₂

    return sum(
        GaussianQuadrature(
            u ->
                (
                    g = g₍₀₎(u);
                    g⁻ = inv(g);
                    𝝊 = sqrt(det(g));
                    𝑁 = [N′(P₁, P₂, I₁, I₂, i, u) for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:d];
                    Q = [sum(𝒂[I₁, I₂, i] * 𝑁[I₁, I₂, j] for I₁ in 1:n₁, I₂ in 1:n₂) for i in 1:d, j in 1:d];
                    sum(
                        sum(C(p, q, m, n, g⁻) * 𝑁[I₁, I₂, p] * Q[i, q] for p in 1:d, q in 1:d) * (sum(Q[o, m] * Q[o, n] for o in 1:d) - g[m, n])
                        for m in 1:d, n in 1:d
                    )
                ) * 𝝊,
            k₁[ι₁]..k₁[ι₁+1],
            k₂[ι₂]..k₂[ι₂+1],
            nip = nip,
        ) for ι₁ in Σ₁, ι₂ in Σ₂
    )
end
