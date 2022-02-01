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

function _newton(M::BSplineManifold{2}, fix_method; nip=NIP)
    𝒂 = controlpoints(M)
    P = bsplinespaces(M)
    n₁, n₂ = dim.(P)
    lineup(I₁, I₂, i) = (i-1)*n₁*n₂ + (I₂-1)*n₁ + (I₁-1) + 1

    t₀ = time()

    H = _matrix_H(M)
    F = _vector_F(M)

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
    M = BSplineManifold(𝒂, P)
    return M, F, Ǧ, t₁ - t₀
end

function _matrix_H(M::BSplineManifold{2,p}) where p
    rrr = StaticArrays.SUnitRange{1,10}()
    𝒂 = controlpoints(M)
    P₁, P₂ = P = bsplinespaces(M)
    p₁, p₂ = p
    k₁, k₂ = k = knotvector.(P)
    l₁, l₂ = length.(k)
    n₁, n₂ = dim.(P)

    H = zeros(n₁,n₂,2,n₁,n₂,2)
    _nodes, _weights = gausslegendre(10)
    nodes = SVector{10,Float64}(_nodes)
    weights = SVector{10,Float64}(_weights)
    nodes₁ = nodes
    nodes₂ = nodes
    weights₁ = weights
    weights₂ = weights
    for s₁ in 1:l₁-1, s₂ in 1:l₂-1
        a₁ = k₁[s₁]
        b₁ = k₁[s₁+1]
        a₂ = k₂[s₂]
        b₂ = k₂[s₂+1]
        w₁ = b₁-a₁
        w₂ = b₂-a₂
        iszero(w₁) && continue
        iszero(w₂) && continue
        dnodes₁ = (w₁ * nodes₁ .+ (a₁+b₁)) / 2
        dnodes₂ = (w₂ * nodes₂ .+ (a₂+b₂)) / 2
        for ii1 in rrr, ii2 in rrr
            u¹,u² = dnodes₁[ii1],dnodes₂[ii2]
            g₁₁ = g₍₀₎₁₁(u¹,u²)
            g₁₂ = g₂₁ = g₍₀₎₁₂(u¹,u²)
            g₂₂ = g₍₀₎₂₂(u¹,u²)
            g = @SMatrix [g₁₁ g₁₂;g₂₁ g₂₂]
            g⁻ = inv(g)
            𝝊 = sqrt(det(g))

            B₁ = bsplinebasisall(P₁,s₁-p₁,u¹)
            B₂ = bsplinebasisall(P₂,s₂-p₂,u²)
            Ḃ₁ = bsplinebasisall(BSplineDerivativeSpace{1}(P₁),s₁-p₁,u¹)
            Ḃ₂ = bsplinebasisall(BSplineDerivativeSpace{1}(P₂),s₂-p₂,u²)

            Q1 = @SVector [sum(𝒂[J₁+(s₁-p₁)-1,J₂+(s₂-p₂)-1,i] * Ḃ₁[J₁]*B₂[J₂] for J₁ in 1:p₁+1, J₂ in 1:p₂+1) for i in 1:2]
            Q2 = @SVector [sum(𝒂[J₁+(s₁-p₁)-1,J₂+(s₂-p₂)-1,i] * B₁[J₁]*Ḃ₂[J₂] for J₁ in 1:p₁+1, J₂ in 1:p₂+1) for i in 1:2]
            Q = hcat(Q1,Q2)
            QQ = @SMatrix [Q[1,m]*Q[1,n] + Q[2,m]*Q[2,n] for m in 1:2, n in 1:2]
            weight1 = weights₁[ii1]
            weight2 = weights₂[ii2]
            C¹¹¹¹ = C(1,1,1,1,g⁻)
            C¹¹¹² = C¹¹²¹ = C¹²¹¹ = C²¹¹¹ = C(1,1,1,2,g⁻)
            C¹¹²² = C²²¹¹ = C(1,1,2,2,g⁻)
            C¹²¹² = C¹²²¹ = C²¹¹² = C²¹²¹ = C(1,2,1,2,g⁻)
            C¹²²² = C²¹²² = C²²¹² = C²²²¹ = C(1,2,2,2,g⁻)
            C²²²² = C(2,2,2,2,g⁻)
            for i₁ in 1:p₁+1, i₂ in 1:p₂+1, i in 1:2, r₁ in 1:p₁+1, r₂ in 1:p₂+1, r in 1:2
                I₁ = i₁+(s₁-p₁)-1
                R₁ = r₁+(s₁-p₁)-1
                I₂ = i₂+(s₂-p₂)-1
                R₂ = r₂+(s₂-p₂)-1

                NI1 = Ḃ₁[i₁]*B₂[i₂]
                NI2 = B₁[i₁]*Ḃ₂[i₂]
                NR1 = Ḃ₁[r₁]*B₂[r₂]
                NR2 = B₁[r₁]*Ḃ₂[r₂]
                s = 0.0
                s += C¹¹¹¹ * NI1 * NR1*Q1[i]*Q1[r]
                s += C¹¹¹² * NI1 * NR2*Q1[i]*Q1[r]
                s += C¹¹²¹ * NI1 * NR1*Q1[i]*Q2[r]
                s += C¹¹²² * NI1 * NR2*Q1[i]*Q2[r]
                s += C¹²¹¹ * NI1 * NR1*Q2[i]*Q1[r]
                s += C¹²¹² * NI1 * NR2*Q2[i]*Q1[r]
                s += C¹²²¹ * NI1 * NR1*Q2[i]*Q2[r]
                s += C¹²²² * NI1 * NR2*Q2[i]*Q2[r]
                s += C²¹¹¹ * NI2 * NR1*Q1[i]*Q1[r]
                s += C²¹¹² * NI2 * NR2*Q1[i]*Q1[r]
                s += C²¹²¹ * NI2 * NR1*Q1[i]*Q2[r]
                s += C²¹²² * NI2 * NR2*Q1[i]*Q2[r]
                s += C²²¹¹ * NI2 * NR1*Q2[i]*Q1[r]
                s += C²²¹² * NI2 * NR2*Q2[i]*Q1[r]
                s += C²²²¹ * NI2 * NR1*Q2[i]*Q2[r]
                s += C²²²² * NI2 * NR2*Q2[i]*Q2[r]
                if i == r
                    s += C¹¹¹¹ * NI1 * NR1 * (QQ[1,1]-g₁₁)/2
                    s += C¹¹¹² * NI1 * NR1 * (QQ[1,2]-g₁₂)/2
                    s += C¹¹²¹ * NI1 * NR1 * (QQ[2,1]-g₂₁)/2
                    s += C¹¹²² * NI1 * NR1 * (QQ[2,2]-g₂₂)/2
                    s += C¹²¹¹ * NI1 * NR2 * (QQ[1,1]-g₁₁)/2
                    s += C¹²¹² * NI1 * NR2 * (QQ[1,2]-g₁₂)/2
                    s += C¹²²¹ * NI1 * NR2 * (QQ[2,1]-g₂₁)/2
                    s += C¹²²² * NI1 * NR2 * (QQ[2,2]-g₂₂)/2
                    s += C²¹¹¹ * NI2 * NR1 * (QQ[1,1]-g₁₁)/2
                    s += C²¹¹² * NI2 * NR1 * (QQ[1,2]-g₁₂)/2
                    s += C²¹²¹ * NI2 * NR1 * (QQ[2,1]-g₂₁)/2
                    s += C²¹²² * NI2 * NR1 * (QQ[2,2]-g₂₂)/2
                    s += C²²¹¹ * NI2 * NR2 * (QQ[1,1]-g₁₁)/2
                    s += C²²¹² * NI2 * NR2 * (QQ[1,2]-g₁₂)/2
                    s += C²²²¹ * NI2 * NR2 * (QQ[2,1]-g₂₁)/2
                    s += C²²²² * NI2 * NR2 * (QQ[2,2]-g₂₂)/2
                end
                s *= 𝝊*weight1*weight2*w₁*w₂/2
                H[I₁, I₂, i, R₁, R₂, r] += s
            end
        end
    end
    return H
end


function _vector_F(M::BSplineManifold{2,p}) where p
    rrr = StaticArrays.SUnitRange{1,10}()
    𝒂 = controlpoints(M)
    P₁, P₂ = P = bsplinespaces(M)
    p₁, p₂ = p
    k₁, k₂ = k = knotvector.(P)
    l₁, l₂ = length.(k)
    n₁, n₂ = dim.(P)

    F = zeros(n₁,n₂,2)
    _nodes, _weights = gausslegendre(10)
    nodes = SVector{10,Float64}(_nodes)
    weights = SVector{10,Float64}(_weights)
    nodes₁ = nodes
    nodes₂ = nodes
    weights₁ = weights
    weights₂ = weights
    for s₁ in 1:l₁-1, s₂ in 1:l₂-1
        a₁ = k₁[s₁]
        b₁ = k₁[s₁+1]
        a₂ = k₂[s₂]
        b₂ = k₂[s₂+1]
        w₁ = b₁-a₁
        w₂ = b₂-a₂
        iszero(w₁) && continue
        iszero(w₂) && continue
        dnodes₁ = (w₁ * nodes₁ .+ (a₁+b₁)) / 2
        dnodes₂ = (w₂ * nodes₂ .+ (a₂+b₂)) / 2
        for ii1 in rrr, ii2 in rrr
            u¹,u² = dnodes₁[ii1],dnodes₂[ii2]
            g₁₁ = g₍₀₎₁₁(u¹,u²)
            g₁₂ = g₂₁ = g₍₀₎₁₂(u¹,u²)
            g₂₂ = g₍₀₎₂₂(u¹,u²)
            g = @SMatrix [g₁₁ g₁₂;g₂₁ g₂₂]
            g⁻ = inv(g)
            𝝊 = sqrt(det(g))

            B₁ = bsplinebasisall(P₁,s₁-p₁,u¹)
            B₂ = bsplinebasisall(P₂,s₂-p₂,u²)
            Ḃ₁ = bsplinebasisall(BSplineDerivativeSpace{1}(P₁),s₁-p₁,u¹)
            Ḃ₂ = bsplinebasisall(BSplineDerivativeSpace{1}(P₂),s₂-p₂,u²)

            Q1 = @SVector [sum(𝒂[J₁+(s₁-p₁)-1,J₂+(s₂-p₂)-1,i] * Ḃ₁[J₁]*B₂[J₂] for J₁ in 1:p₁+1, J₂ in 1:p₂+1) for i in 1:2]
            Q2 = @SVector [sum(𝒂[J₁+(s₁-p₁)-1,J₂+(s₂-p₂)-1,i] * B₁[J₁]*Ḃ₂[J₂] for J₁ in 1:p₁+1, J₂ in 1:p₂+1) for i in 1:2]
            Q = hcat(Q1,Q2)
            QQ = @SMatrix [Q[1,m]*Q[1,n] + Q[2,m]*Q[2,n] for m in 1:2, n in 1:2]
            weight1 = weights₁[ii1]
            weight2 = weights₂[ii2]
            C¹¹¹¹ = C(1,1,1,1,g⁻)
            C¹¹¹² = C¹¹²¹ = C¹²¹¹ = C²¹¹¹ = C(1,1,1,2,g⁻)
            C¹¹²² = C²²¹¹ = C(1,1,2,2,g⁻)
            C¹²¹² = C¹²²¹ = C²¹¹² = C²¹²¹ = C(1,2,1,2,g⁻)
            C¹²²² = C²¹²² = C²²¹² = C²²²¹ = C(1,2,2,2,g⁻)
            C²²²² = C(2,2,2,2,g⁻)
            for i₁ in 1:p₁+1, i₂ in 1:p₂+1, i in 1:2
                I₁ = i₁+(s₁-p₁)-1
                I₂ = i₂+(s₂-p₂)-1

                NI1 = Ḃ₁[i₁]*B₂[i₂]
                NI2 = B₁[i₁]*Ḃ₂[i₂]
                s = 0.0
                s += C¹¹¹¹ * NI1 * Q[i,1] * (QQ[1,1]-g₁₁)/2
                s += C¹¹¹² * NI1 * Q[i,1] * (QQ[1,2]-g₁₂)/2
                s += C¹¹²¹ * NI1 * Q[i,1] * (QQ[2,1]-g₂₁)/2
                s += C¹¹²² * NI1 * Q[i,1] * (QQ[2,2]-g₂₂)/2
                s += C¹²¹¹ * NI1 * Q[i,2] * (QQ[1,1]-g₁₁)/2
                s += C¹²¹² * NI1 * Q[i,2] * (QQ[1,2]-g₁₂)/2
                s += C¹²²¹ * NI1 * Q[i,2] * (QQ[2,1]-g₂₁)/2
                s += C¹²²² * NI1 * Q[i,2] * (QQ[2,2]-g₂₂)/2
                s += C²¹¹¹ * NI2 * Q[i,1] * (QQ[1,1]-g₁₁)/2
                s += C²¹¹² * NI2 * Q[i,1] * (QQ[1,2]-g₁₂)/2
                s += C²¹²¹ * NI2 * Q[i,1] * (QQ[2,1]-g₂₁)/2
                s += C²¹²² * NI2 * Q[i,1] * (QQ[2,2]-g₂₂)/2
                s += C²²¹¹ * NI2 * Q[i,2] * (QQ[1,1]-g₁₁)/2
                s += C²²¹² * NI2 * Q[i,2] * (QQ[1,2]-g₁₂)/2
                s += C²²²¹ * NI2 * Q[i,2] * (QQ[2,1]-g₂₁)/2
                s += C²²²² * NI2 * Q[i,2] * (QQ[2,2]-g₂₂)/2
                s *= 𝝊*weight1*weight2*w₁*w₂/2
                F[I₁, I₂, i] += s
            end
        end
    end
    return F
end
