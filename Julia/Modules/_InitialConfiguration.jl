export InitialConfiguration
function InitialConfiguration(D;n₁=15,nip=NIP)
    if (isfile(DIR*"/"*NAME*".jld"))
        error("jld file already exists")
    end
    mkpath(DIR)
    mkpath(DIR*"/nurbs")
    mkpath(DIR*"/strain")
    mkpath(DIR*"/colorbar")
    mkpath(DIR*"/slack")
    BsJLD=Dict{String,Any}("Expr"=>EXPR)

    M=InitBs(D,n₁,nip=nip)
    comment="Initial Configuration"
    BsTree=Tree()

    Export(M,BsTree,BsJLD,comment=comment)
end

function InitBs(D,n₁;nip=NIP)::BSplineManifold
    D₁,D₂=D

    c(t)=[t,sum(extrema(D[2]))/2] # 中心線に沿った座標
    ṡ₍₀₎(t)=sqrt(g₍₀₎(c(t))[1,1])
    s̈₍₀₎(t)=ForwardDiff.derivative(ṡ₍₀₎,t)
    𝜅₍₀₎(t)=𝛤₍₀₎²₁₁(c(t))*𝝊₍₀₎(c(t))/ṡ₍₀₎(t)^3 # 測地的曲率

    function ode(𝒄̇𝒄̈,𝒄𝒄̇,par,t)
        𝒄̇𝒄̈[1]=𝒄𝒄̇[3]
        𝒄̇𝒄̈[2]=𝒄𝒄̇[4]
        𝒄̇𝒄̈[3]=dot([s̈₍₀₎(t)/ṡ₍₀₎(t),-𝜅₍₀₎(t)*ṡ₍₀₎(t)],𝒄𝒄̇[3:4])
        𝒄̇𝒄̈[4]=dot([𝜅₍₀₎(t)*ṡ₍₀₎(t),s̈₍₀₎(t)/ṡ₍₀₎(t)],𝒄𝒄̇[3:4])
    end
    𝒄𝒄̇₀=vcat([0.0,0.0],[1.,0.]*ṡ₍₀₎(minimum(D₁)))
    sol=solve(ODEProblem(ode,𝒄𝒄̇₀,extrema(D₁)))
    𝒄(t)=sol(t)[1:d] # 解となる中心曲線
    𝒄₁(t)=sol(t)[(d+1):(2d)] # その導関数
    𝒄₂(t)=[g₍₀₎₁₂(c(t)) -𝝊₍₀₎(c(t));𝝊₍₀₎(c(t)) g₍₀₎₁₂(c(t))]*𝒄₁(t)/g₍₀₎₁₁(c(t)) # 中心曲線上の幅方向のベクトル場

    p₁=3
    k₁=Knots(sort(vcat(repeat(collect(extrema(D₁)),inner=p₁),collect(range(leftendpoint(D₁),stop=rightendpoint(D₁),length=n₁-2)))))
    P₁=BSplineSpace(p₁,k₁)

    global 𝒎=FittingBSpline(𝒄,P₁,nip=nip)
    global 𝒓=FittingBSpline(𝒄₂,P₁,nip=nip)
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