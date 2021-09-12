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

"""
    newton_onestep(; fixingmethod=:default, parent::Int=0, nip=NIP)

Compute one step of Newton-Raphson method
"""
function newton_onestep(; fixingmethod=:default, parent::Int=0, nip=NIP)
    if fixingmethod == :default
        fixed = DefaultOrientation
    elseif fixingmethod == :fix3points
        fixed = FixThreePoints
    else
        error("No method for $(fixingmethod)")
    end
    parent = _realparent(parent)
    M = loadM(index=parent)

    n₁, n₂ = dim.(bsplinespaces(M))

    iseven(n₁) && error("n₁ should be odd numbers")
    iseven(n₂) && error("n₂ should be odd numbers")

    M = _positioning(M)
    M, F, Ǧ, Δt = _newton(M, fixed, nip=nip)
    comment =
        "Newton onestep - residual norm: " *
        (@sprintf("%.4e", norm(F))) *
        ", Δa norm: " *
        (@sprintf("%.4e", norm(Ǧ))) *
        ", computation time: " *
        SecondsToString(Δt)
    _export(M, parent, comment=comment)
end

function _newton(M::AbstractBSplineManifold, fix_method; nip=NIP)
    𝒂 = controlpoints(M)
    P = collect(bsplinespaces(M))
    n₁, n₂ = dim.(P)
    lineup(I₁, I₂, i) = (i-1)*n₁*n₂ + (I₂-1)*n₁ + (I₁-1) + 1

    t₀ = time()

    H = zeros(n₁,n₂,2,n₁,n₂,2)
    F = zeros(n₁,n₂,2)
    Threads.@threads for I₁ in 1:n₁
        for I₂ in 1:n₂, i in 1:2, R₁ in 1:n₁, R₂ in 1:n₂, r in 1:2
            if lineup(I₁, I₂, i) ≤ lineup(R₁, R₂, r)
                H[I₁, I₂, i, R₁, R₂, r] = H[R₁, R₂, r, I₁, I₂, i] = elm_H(M, I₁, I₂, i, R₁, R₂, r, nip=nip)
            end
        end
    end
    Threads.@threads for I₁ in 1:n₁
        for I₂ in 1:n₂, i in 1:2
            F[I₁, I₂, i] = elm_F(M, I₁, I₂, i, nip=nip)
        end
    end

    t₁ = time()

    N = 2n₁*n₂
    _fixed = sort(collect((i -> lineup(i...)).(fix_method(n₁, n₂))))
    _unfixed = deleteat!(collect(1:N), _fixed)

    F = reshape(F, N)
    H = reshape(H, N, N)
    𝒂 = 𝒂ₒ = reshape(𝒂, N)
    Ȟ = H[_unfixed, _unfixed]
    𝒂̌ = 𝒂[_unfixed]
    F̌ = F[_unfixed]
    Ǧ = Ȟ \ F̌
    𝒂̌ = 𝒂̌ - Ǧ
    for i in _fixed
        insert!(𝒂̌, i, 𝒂ₒ[i])
    end
    𝒂 = reshape(𝒂̌, n₁, n₂, 2)
    M = typeof(M)(P, 𝒂)
    return M, F, Ǧ, t₁ - t₀
end

function elm_H(M::AbstractBSplineManifold, I₁, I₂, i, R₁, R₂, r; nip=NIP)
    𝒂 = controlpoints(M)
    P₁, P₂ = P = bsplinespaces(M)
    p₁, p₂ = degree.(P)
    k₁, k₂ = knots.(P)
    n₁, n₂ = dim.(P)

    Σ₁ = max(I₁, R₁):min(I₁, R₁)+p₁
    Σ₂ = max(I₂, R₂):min(I₂, R₂)+p₂

    if length(Σ₁) == 0 || length(Σ₂) == 0
        return 0.0
    elseif i == r
        return sum(
            GaussianQuadrature(
                (u¹,u²) ->
                    (
                        g = g₍₀₎(u¹,u²);
                        g⁻ = inv(g);
                        𝝊 = sqrt(det(g));
                        𝑁 = [N′(P₁, P₂, I₁, I₂, i, u¹,u²) for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:2];
                        Q = @SMatrix [sum(𝒂[I₁, I₂, i]*𝑁[I₁, I₂, j] for I₁ in 1:n₁, I₂ in 1:n₂) for i in 1:2, j in 1:2];
                        QQ = @SMatrix [Q[1, m]*Q[1, n]+Q[2, m]*Q[2, n] for m in 1:2, n in 1:2];
                        sum(
                            C(p, q, m, n, g⁻)*
                            𝑁[I₁, I₂, p]*
                            (𝑁[R₁, R₂, q] * (QQ[m, n] - g[m, n]) + 2𝑁[R₁, R₂, n]*Q[i, q]*Q[r, m])
                            for p in 1:2, q in 1:2, m in 1:2, n in 1:2
                        )*𝝊
                    ),
                k₁[s₁]..k₁[s₁+1],
                k₂[s₂]..k₂[s₂+1],
                nip=nip,
            ) for s₁ in Σ₁, s₂ in Σ₂
        )
    else
        return sum(
            GaussianQuadrature(
                (u¹,u²) ->
                    (
                        g = g₍₀₎(u¹,u²);
                        g⁻ = inv(g);
                        𝝊 = sqrt(det(g));
                        𝑁 = [N′(P₁, P₂, I₁, I₂, i, u¹,u²) for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:2];
                        Q = @SMatrix [sum(𝒂[I₁, I₂, i]*𝑁[I₁, I₂, j] for I₁ in 1:n₁, I₂ in 1:n₂) for i in 1:2, j in 1:2];
                        2*sum(
                            C(p, q, m, n, g⁻)*𝑁[I₁, I₂, p]*𝑁[R₁, R₂, n]*Q[i, q]*Q[r, m]
                            for p in 1:2, q in 1:2, m in 1:2, n in 1:2
                        )*𝝊
                    ),
                k₁[s₁]..k₁[s₁+1],
                k₂[s₂]..k₂[s₂+1],
                nip=nip,
            ) for s₁ in Σ₁, s₂ in Σ₂
        )
    end
end

function elm_F(M::AbstractBSplineManifold, I₁, I₂, i; nip = NIP)
    𝒂 = controlpoints(M)
    P₁, P₂ = P = bsplinespaces(M)
    p₁, p₂ = degree.(P)
    k₁, k₂ = knots.(P)
    n₁, n₂ = dim.(P)

    Σ₁ = I₁:I₁+p₁
    Σ₂ = I₂:I₂+p₂

    return sum(
        GaussianQuadrature(
            (u¹,u²) ->
                (
                    g = g₍₀₎(u¹,u²);
                    g⁻ = inv(g);
                    𝝊 = sqrt(det(g));
                    𝑁 = [N′(P₁, P₂, I₁, I₂, i, u¹,u²) for I₁ in 1:n₁, I₂ in 1:n₂, i in 1:2];
                    Q = @SMatrix [sum(𝒂[I₁, I₂, i] * 𝑁[I₁, I₂, j] for I₁ in 1:n₁, I₂ in 1:n₂) for i in 1:2, j in 1:2];
                    QQ = @SMatrix [Q[1, m]*Q[1, n]+Q[2, m]*Q[2, n] for m in 1:2, n in 1:2];
                    sum(
                        sum(C(p, q, m, n, g⁻) * 𝑁[I₁, I₂, p] * Q[i, q] for p in 1:2, q in 1:2) * (QQ[m, n] - g[m, n])
                        for m in 1:2, n in 1:2
                    )
                ) * 𝝊,
            k₁[s₁]..k₁[s₁+1],
            k₂[s₂]..k₂[s₂+1],
            nip = nip,
        ) for s₁ in Σ₁, s₂ in Σ₂
    )
end
