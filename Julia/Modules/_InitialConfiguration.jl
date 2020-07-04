using DifferentialEquations

export InitialConfiguration
function InitialConfiguration(D;n₁=15,nip=NIP)
    parent = 0

    D₁,D₂ = D
    M = InitBs(D,n₁,nip=nip)
    comment = "Initial Configuration - domain: "*repr([endpoints(D₁)...])*"×"*repr([endpoints(D₂)...])

    Export(M,parent,comment=comment)
end

function InitBs(D,n₁;nip=NIP)::FastBSplineManifold
    D₁,D₂ = D

    function ode(𝒄̇𝒄̈,𝒄𝒄̇,par,t)
        𝒄̇𝒄̈[1] = 𝒄𝒄̇[3]
        𝒄̇𝒄̈[2] = 𝒄𝒄̇[4]
        𝒄̇𝒄̈[3] = dot([s̈₍₀₎(D₂,t)/ṡ₍₀₎(D₂,t),-𝜅₍₀₎(D₂,t)*ṡ₍₀₎(D₂,t)],𝒄𝒄̇[3:4])
        𝒄̇𝒄̈[4] = dot([𝜅₍₀₎(D₂,t)*ṡ₍₀₎(D₂,t),s̈₍₀₎(D₂,t)/ṡ₍₀₎(D₂,t)],𝒄𝒄̇[3:4])
    end
    𝒄𝒄̇₀ = vcat([0.0,0.0],[1.,0.]*ṡ₍₀₎(D₂,minimum(D₁)))
    curve = solve(ODEProblem(ode,𝒄𝒄̇₀,extrema(D₁)))
    𝒄(t) = curve(t)[1:d] # center curve of the solution
    𝒄₁(t) = curve(t)[(1:d).+d] # its derivative
    𝒄₂(t) = [g₍₀₎₁₂(c(D₂,t)) -𝝊₍₀₎(c(D₂,t));𝝊₍₀₎(c(D₂,t)) g₍₀₎₁₂(c(D₂,t))]*𝒄₁(t)/g₍₀₎₁₁(c(D₂,t)) # 中心曲線上の幅方向のベクトル場

    D₁₋, D₁₊ = extrema(D₁)
    p₁ = 3
    k₁ = Knots(range(D₁₋, D₁₊, length=n₁-p₁+1)) + p₁ * Knots(D₁₋, D₁₊)
    P₁ = FastBSplineSpace(p₁,k₁)

    𝒎 = FittingControlPoints(t->𝒄(t[1]),[P₁])
    𝒓 = FittingControlPoints(t->𝒄₂(t[1]),[P₁])
    a1 = 𝒎-width(D₂)*𝒓/2
    a2 = 𝒎+width(D₂)*𝒓/2
    p₂ = 1
    k₂ = Knots(repeat(collect(extrema(D₂)),inner=2))
    n₂ = length(k₂)-p₂-1

    P₂ = FastBSplineSpace(p₂,k₂)
    𝒂 = hcat(a1,a2)

    M = FastBSplineManifold([P₁,P₂],𝒂)
    M′ = refinement(M,p₊=[0,1])
    return Positioning(M′)
end
