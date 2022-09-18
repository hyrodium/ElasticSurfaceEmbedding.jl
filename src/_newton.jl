function _defaultorientation(n₁, n₂)
    return ([(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2, 1], [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2, 2], [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2 - 1, 1])
end

function _fixthreepoints(n₁, n₂)
    return (
        [1, (n₂ + 1) ÷ 2, 1],
        [1, (n₂ + 1) ÷ 2, 2],
        [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2, 1],
        [(n₁ + 1) ÷ 2, (n₂ + 1) ÷ 2, 2],
        [n₁, (n₂ + 1) ÷ 2, 1],
        [n₁, (n₂ + 1) ÷ 2, 2],
    )
end

_abbstr(t::Week) = string(t.value) * "w "
_abbstr(t::Day) = string(t.value) * "d "
_abbstr(t::Hour) = string(t.value) * "h "
_abbstr(t::Minute) = string(t.value) * "m "
_abbstr(t::Second) = string(t.value) * "s "
_abbstr(t::Millisecond) = string(t.value) * "ms "
_abbstr(t::Vector{Period}) = *(_abbstr.(t)...)[1:end-1]

function _seconds2string(Δt::Float64)
    prds = Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(floor(1000Δt)))).periods
    return _abbstr(prds)
end

"""
    newton_onestep(allsteps; fixingmethod=:default, parent::Int=0)

Compute one step of Newton-Raphson method
"""
function newton_onestep!(allsteps; fixingmethod=:default, parent::Int=0)
    if fixingmethod == :default
        fixed = _defaultorientation
    elseif fixingmethod == :fix3points
        fixed = _fixthreepoints
    else
        error("No method for $(fixingmethod). Use :default or :fix3points.")
    end
    parent = _validindex(allsteps, parent)
    M = loadM(allsteps, index=parent)

    n₁, n₂ = dim.(bsplinespaces(M))

    iseven(n₁) && error("n₁ should be odd numbers")
    iseven(n₂) && error("n₂ should be odd numbers")

    M = _positioning(M)
    M, F, Ǧ, Δt = _newton(M, fixed)
    comment =
        "Newton onestep - residual norm: " *
        (@sprintf("%.4e", norm(F))) *
        ", Δa norm: " *
        (@sprintf("%.4e", norm(Ǧ))) *
        ", computation time: " *
        _seconds2string(Δt)
    info = Dict(["type"=>"newton", "fixingmethod"=>string(fixingmethod)])
    step = Step(M, comment, info)
    addstep!(allsteps, step, parent)
end

function _newton(M::BSplineManifold{2, p, <:SVector}, fix_method) where p
    𝒂 = _arrayofvector2array(controlpoints(M))
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
    M = BSplineManifold(_array2arrayofvector(𝒂), P)
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

            Q₁ = sum(𝒂[J₁+(s₁-p₁)-1,J₂+(s₂-p₂)-1] * Ḃ₁[J₁]*B₂[J₂] for J₁ in 1:p₁+1, J₂ in 1:p₂+1)
            Q₂ = sum(𝒂[J₁+(s₁-p₁)-1,J₂+(s₂-p₂)-1] * B₁[J₁]*Ḃ₂[J₂] for J₁ in 1:p₁+1, J₂ in 1:p₂+1)
            Q = hcat(Q₁,Q₂)
            QQ = @SMatrix [Q[1,m]*Q[1,n] + Q[2,m]*Q[2,n] for m in 1:2, n in 1:2]
            weight1 = weights₁[ii1]
            weight2 = weights₂[ii2]
            C¹¹¹¹ = C(1,1,1,1,g⁻)
            C¹¹¹² = C(1,1,1,2,g⁻)
            C¹¹²² = C(1,1,2,2,g⁻)
            C¹²¹² = C(1,2,1,2,g⁻)
            C¹²²² = C(1,2,2,2,g⁻)
            C²²²² = C(2,2,2,2,g⁻)
            C¹¹²¹ = C¹²¹¹ = C²¹¹¹ = C¹¹¹²
            C²²¹¹ = C¹¹²²
            C¹²²¹ = C²¹¹² = C²¹²¹ = C¹²¹²
            C²¹²² = C²²¹² = C²²²¹ = C¹²²²

            for i₁ in 1:p₁+1, i₂ in 1:p₂+1, i in 1:2, r₁ in 1:p₁+1, r₂ in 1:p₂+1, r in 1:2
                I₁ = i₁+(s₁-p₁)-1
                R₁ = r₁+(s₁-p₁)-1
                I₂ = i₂+(s₂-p₂)-1
                R₂ = r₂+(s₂-p₂)-1

                Ni₁ = Ḃ₁[i₁]*B₂[i₂]
                Ni₂ = B₁[i₁]*Ḃ₂[i₂]
                Nr₁ = Ḃ₁[r₁]*B₂[r₂]
                Nr₂ = B₁[r₁]*Ḃ₂[r₂]
                s =  C¹¹¹¹ * Ni₁ * Nr₁ * Q₁[i] * Q₁[r]
                s += C¹¹¹² * Ni₁ * Nr₂ * Q₁[i] * Q₁[r]
                s += C¹¹²¹ * Ni₁ * Nr₁ * Q₁[i] * Q₂[r]
                s += C¹¹²² * Ni₁ * Nr₂ * Q₁[i] * Q₂[r]
                s += C¹²¹¹ * Ni₁ * Nr₁ * Q₂[i] * Q₁[r]
                s += C¹²¹² * Ni₁ * Nr₂ * Q₂[i] * Q₁[r]
                s += C¹²²¹ * Ni₁ * Nr₁ * Q₂[i] * Q₂[r]
                s += C¹²²² * Ni₁ * Nr₂ * Q₂[i] * Q₂[r]
                s += C²¹¹¹ * Ni₂ * Nr₁ * Q₁[i] * Q₁[r]
                s += C²¹¹² * Ni₂ * Nr₂ * Q₁[i] * Q₁[r]
                s += C²¹²¹ * Ni₂ * Nr₁ * Q₁[i] * Q₂[r]
                s += C²¹²² * Ni₂ * Nr₂ * Q₁[i] * Q₂[r]
                s += C²²¹¹ * Ni₂ * Nr₁ * Q₂[i] * Q₁[r]
                s += C²²¹² * Ni₂ * Nr₂ * Q₂[i] * Q₁[r]
                s += C²²²¹ * Ni₂ * Nr₁ * Q₂[i] * Q₂[r]
                s += C²²²² * Ni₂ * Nr₂ * Q₂[i] * Q₂[r]
                if i == r
                    s += C¹¹¹¹ * Ni₁ * Nr₁ * (QQ[1,1]-g₁₁)/2
                    s += C¹¹¹² * Ni₁ * Nr₁ * (QQ[1,2]-g₁₂)/2
                    s += C¹¹²¹ * Ni₁ * Nr₁ * (QQ[2,1]-g₂₁)/2
                    s += C¹¹²² * Ni₁ * Nr₁ * (QQ[2,2]-g₂₂)/2
                    s += C¹²¹¹ * Ni₁ * Nr₂ * (QQ[1,1]-g₁₁)/2
                    s += C¹²¹² * Ni₁ * Nr₂ * (QQ[1,2]-g₁₂)/2
                    s += C¹²²¹ * Ni₁ * Nr₂ * (QQ[2,1]-g₂₁)/2
                    s += C¹²²² * Ni₁ * Nr₂ * (QQ[2,2]-g₂₂)/2
                    s += C²¹¹¹ * Ni₂ * Nr₁ * (QQ[1,1]-g₁₁)/2
                    s += C²¹¹² * Ni₂ * Nr₁ * (QQ[1,2]-g₁₂)/2
                    s += C²¹²¹ * Ni₂ * Nr₁ * (QQ[2,1]-g₂₁)/2
                    s += C²¹²² * Ni₂ * Nr₁ * (QQ[2,2]-g₂₂)/2
                    s += C²²¹¹ * Ni₂ * Nr₂ * (QQ[1,1]-g₁₁)/2
                    s += C²²¹² * Ni₂ * Nr₂ * (QQ[1,2]-g₁₂)/2
                    s += C²²²¹ * Ni₂ * Nr₂ * (QQ[2,1]-g₂₁)/2
                    s += C²²²² * Ni₂ * Nr₂ * (QQ[2,2]-g₂₂)/2
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

            Q₁ = sum(𝒂[J₁+(s₁-p₁)-1,J₂+(s₂-p₂)-1] * Ḃ₁[J₁]*B₂[J₂] for J₁ in 1:p₁+1, J₂ in 1:p₂+1)
            Q₂ = sum(𝒂[J₁+(s₁-p₁)-1,J₂+(s₂-p₂)-1] * B₁[J₁]*Ḃ₂[J₂] for J₁ in 1:p₁+1, J₂ in 1:p₂+1)
            Q = hcat(Q₁,Q₂)
            QQ = @SMatrix [Q[1,m]*Q[1,n] + Q[2,m]*Q[2,n] for m in 1:2, n in 1:2]
            weight1 = weights₁[ii1]
            weight2 = weights₂[ii2]
            C¹¹¹¹ = C(1,1,1,1,g⁻)
            C¹¹¹² = C(1,1,1,2,g⁻)
            C¹¹²² = C(1,1,2,2,g⁻)
            C¹²¹² = C(1,2,1,2,g⁻)
            C¹²²² = C(1,2,2,2,g⁻)
            C²²²² = C(2,2,2,2,g⁻)
            C¹¹²¹ = C¹²¹¹ = C²¹¹¹ = C¹¹¹²
            C²²¹¹ = C¹¹²²
            C¹²²¹ = C²¹¹² = C²¹²¹ = C¹²¹²
            C²¹²² = C²²¹² = C²²²¹ = C¹²²²
            for i₁ in 1:p₁+1, i₂ in 1:p₂+1, i in 1:2
                I₁ = i₁+(s₁-p₁)-1
                I₂ = i₂+(s₂-p₂)-1

                Ni₁ = Ḃ₁[i₁]*B₂[i₂]
                Ni₂ = B₁[i₁]*Ḃ₂[i₂]
                s =  C¹¹¹¹ * Ni₁ * Q₁[i] * (QQ[1,1]-g₁₁)/2
                s += C¹¹¹² * Ni₁ * Q₁[i] * (QQ[1,2]-g₁₂)/2
                s += C¹¹²¹ * Ni₁ * Q₁[i] * (QQ[2,1]-g₂₁)/2
                s += C¹¹²² * Ni₁ * Q₁[i] * (QQ[2,2]-g₂₂)/2
                s += C¹²¹¹ * Ni₁ * Q₂[i] * (QQ[1,1]-g₁₁)/2
                s += C¹²¹² * Ni₁ * Q₂[i] * (QQ[1,2]-g₁₂)/2
                s += C¹²²¹ * Ni₁ * Q₂[i] * (QQ[2,1]-g₂₁)/2
                s += C¹²²² * Ni₁ * Q₂[i] * (QQ[2,2]-g₂₂)/2
                s += C²¹¹¹ * Ni₂ * Q₁[i] * (QQ[1,1]-g₁₁)/2
                s += C²¹¹² * Ni₂ * Q₁[i] * (QQ[1,2]-g₁₂)/2
                s += C²¹²¹ * Ni₂ * Q₁[i] * (QQ[2,1]-g₂₁)/2
                s += C²¹²² * Ni₂ * Q₁[i] * (QQ[2,2]-g₂₂)/2
                s += C²²¹¹ * Ni₂ * Q₂[i] * (QQ[1,1]-g₁₁)/2
                s += C²²¹² * Ni₂ * Q₂[i] * (QQ[1,2]-g₁₂)/2
                s += C²²²¹ * Ni₂ * Q₂[i] * (QQ[2,1]-g₂₁)/2
                s += C²²²² * Ni₂ * Q₂[i] * (QQ[2,2]-g₂₂)/2
                s *= 𝝊*weight1*weight2*w₁*w₂/2
                F[I₁, I₂, i] += s
            end
        end
    end
    return F
end
