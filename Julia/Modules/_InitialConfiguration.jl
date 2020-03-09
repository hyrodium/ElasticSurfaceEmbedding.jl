using DifferentialEquations

export InitialConfiguration
function InitialConfiguration(D;n₁=15,nip=NIP)
    parent=0
    mkpath(DIR)
    mkpath(DIR*"/nurbs")
    mkpath(DIR*"/strain")
    mkpath(DIR*"/colorbar")
    mkpath(DIR*"/slack")

    D₁,D₂=D
    M=InitBs(D,n₁,nip=nip)
    comment="Initial Configuration with a domain "*repr([endpoints(D₁)...])*"×"*repr([endpoints(D₂)...])

    Export(M,parent,comment=comment)
end

function InitBs(D,n₁;nip=NIP)::BSplineManifold
    D₁,D₂=D

    function ode(𝒄̇𝒄̈,𝒄𝒄̇,par,t)
        𝒄̇𝒄̈[1]=𝒄𝒄̇[3]
        𝒄̇𝒄̈[2]=𝒄𝒄̇[4]
        𝒄̇𝒄̈[3]=dot([s̈₍₀₎(D₂,t)/ṡ₍₀₎(D₂,t),-𝜅₍₀₎(D₂,t)*ṡ₍₀₎(D₂,t)],𝒄𝒄̇[3:4])
        𝒄̇𝒄̈[4]=dot([𝜅₍₀₎(D₂,t)*ṡ₍₀₎(D₂,t),s̈₍₀₎(D₂,t)/ṡ₍₀₎(D₂,t)],𝒄𝒄̇[3:4])
    end
    𝒄𝒄̇₀=vcat([0.0,0.0],[1.,0.]*ṡ₍₀₎(D₂,minimum(D₁)))
    curve=solve(ODEProblem(ode,𝒄𝒄̇₀,extrema(D₁)))
    𝒄(t)=curve(t)[1:d] # center curve of the solution
    𝒄₁(t)=curve(t)[(1:d).+d] # its derivative
    𝒄₂(t)=[g₍₀₎₁₂(c(D₂,t)) -𝝊₍₀₎(c(D₂,t));𝝊₍₀₎(c(D₂,t)) g₍₀₎₁₂(c(D₂,t))]*𝒄₁(t)/g₍₀₎₁₁(c(D₂,t)) # 中心曲線上の幅方向のベクトル場

    p₁=3
    k₁=Knots(sort(vcat(repeat(collect(extrema(D₁)),inner=p₁),collect(range(leftendpoint(D₁),stop=rightendpoint(D₁),length=n₁-2)))))
    P₁=BSplineSpace(p₁,k₁)

    𝒎=FittingBSpline(𝒄,P₁,nip=nip)
    𝒓=FittingBSpline(𝒄₂,P₁,nip=nip)
    a1=𝒎-width(D₂)*𝒓/2
    a2=𝒎+width(D₂)*𝒓/2
    p₂=1
    k₂=Knots(repeat(collect(extrema(D₂)),inner=2))
    n₂=length(k₂)-p₂-1

    P₂=BSplineSpace(p₂,k₂)
    𝒂=[[a1[I₁][i],a2[I₁][i]][I₂] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d]
    M=BSplineManifold([P₁,P₂],𝒂)
    M′=BSpline.Refinement(M,p₊=[0,1])
    return Positioning(M′)
end
