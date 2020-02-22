module ElasticSurfaceEmbedding

using LinearAlgebra
using Printf
using Distributed
using IntervalSets
using ForwardDiff
using FastGaussQuadrature
using Dates
using DifferentialEquations
using JLD

using BSpline
using Slack
using ParametricDraw

export @DefineShape, InitialConfiguration, Refinement, NewtonMethodIteration, FinalOutput, ShowKnots, ShowMaximumStrain, Settings, Restoration

include("../Config.jl")
const d=2 # Dimension
const Y=1.0 # Young率Y
const 𝝀=𝝂*Y/((1+𝝂)*(1-(d-1)*𝝂)) # Lamé constant
const 𝝁=1/2(1+𝝂) # Lamé constant

macro DefineShape(ex)
    global EXPR=ex
    return :(@everywhere $ex)
end

𝒑′₍₀₎(u)=ForwardDiff.jacobian(Main.𝒑₍₀₎,u) # 接ベクトル
𝒑₁₍₀₎(u)=ForwardDiff.derivative(u₁->Main.𝒑₍₀₎([u₁,u[2]]),u[1])
𝒑₂₍₀₎(u)=ForwardDiff.derivative(u₂->Main.𝒑₍₀₎([u[1],u₂]),u[2])
𝒑₁₁₍₀₎(u)=ForwardDiff.derivative(u₁->Main.𝒑₁₍₀₎([u₁,u[2]]),u[1])
𝒑₁₂₍₀₎(u)=ForwardDiff.derivative(u₂->Main.𝒑₁₍₀₎([u[1],u₂]),u[2])
𝒑₂₁₍₀₎(u)=ForwardDiff.derivative(u₁->Main.𝒑₂₍₀₎([u₁,u[2]]),u[1])
𝒑₂₂₍₀₎(u)=ForwardDiff.derivative(u₂->Main.𝒑₂₍₀₎([u[1],u₂]),u[2])
𝒆₍₀₎(u)=normalize(cross(𝒑₁₍₀₎(u),𝒑₂₍₀₎(u)))
g₍₀₎(u)=𝒑′₍₀₎(u)'𝒑′₍₀₎(u) # 第1基本量
g₍₀₎₁₁(u)=𝒑₁₍₀₎(u)'𝒑₁₍₀₎(u) # 第1基本量
g₍₀₎₁₂(u)=𝒑₁₍₀₎(u)'𝒑₂₍₀₎(u) # 第1基本量
g₍₀₎₂₁(u)=𝒑₂₍₀₎(u)'𝒑₁₍₀₎(u) # 第1基本量
g₍₀₎₂₂(u)=𝒑₂₍₀₎(u)'𝒑₂₍₀₎(u) # 第1基本量
h₍₀₎(u)=[(𝒆₍₀₎(u)'*𝒑₁₁₍₀₎(u)) (𝒆₍₀₎(u)'*𝒑₁₂₍₀₎(u)) ; (𝒆₍₀₎(u)'*𝒑₂₁₍₀₎(u)) (𝒆₍₀₎(u)'*𝒑₂₂₍₀₎(u))] # 第2基本量
K₍₀₎(u)=det(h₍₀₎(u))/det(g₍₀₎(u))
𝝊₍₀₎(u)=norm(cross(𝒑₁₍₀₎(u),𝒑₂₍₀₎(u))) # 体積要素υ
g⁻₍₀₎(u)=inv(g₍₀₎(u)) # 第一基本量の逆


function GaussianQuadrature(f,D₁,D₂;nip=NIP)
    nodes, weights = gausslegendre(nip)
    return sum(
    (weights*weights').*
    [f([x,y]) for
            x ∈ (width(D₁)*nodes.+sum(extrema(D₁)))/2,
            y ∈ (width(D₂)*nodes.+sum(extrema(D₂)))/2
    ])*width(D₁)*width(D₂)/4
end

function GaussianQuadrature(f,D;nip=NIP)
    nodes, weights = gausslegendre(nip)
    return sum(
    weights.*
    [f(x) for x ∈ (width(D)*nodes.+sum(extrema(D)))/2
    ])*width(D)
end

function FittingBSpline(f, P::BSplineSpace; nip=NIP)
    p=P.degree
    k=P.knots
    D=k[1+p]..k[end-p]
    function a(i,j)
        D′=(max(k[i],k[j])..min(k[i+p+1],k[j+p+1])) ∩ D
        if width(D′)==0
            return 0
        else
            return GaussianQuadrature(t->BSplineBasis(i,P,t)*BSplineBasis(j,P,t), D′)
        end
    end
    n=dim(P)
    A=[a(i,j) for i ∈ 1:n, j ∈ 1:n]
    b=[GaussianQuadrature(t->BSplineBasis(i,P,t)*f(t), ((k[i]..k[i+p+1]) ∩ D)) for i ∈ 1:n]
    return inv(A)*b
end

function aff(𝒂::Array{Float64,3},A::Array{Float64,2},b::Array{Float64,1})
    #x'=Ax+b
    n₁,n₂=size(𝒂)[1:d]
    return [sum(A[i,j]*𝒂[I₁,I₂,j] for j ∈ 1:d)+b[i] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d]
end

function Positioning(𝒂::Array{Float64,3}) # 制御点の位置調整
    n₁,n₂=size(𝒂)[1:d]
    ind0=[(n₁+1)÷2,(n₂+1)÷2]
    ind1=ind0-[0,1]
    v=𝒂[ind1...,:]-𝒂[ind0...,:]
    R=-[v[2] -v[1];v[1] v[2]]/norm(v)
    return 𝒂=aff(𝒂,R,-R*𝒂[ind0...,:])
end

function Positioning(M::BSplineManifold) # 制御点の位置調整
    𝒫s = M.bsplinespaces
    𝒂 = M.controlpoints
    if (length(𝒫s) ≠ d)
        error("dimension does not match")
    end

    p¹,p²=p=[M.bsplinespaces[i].degree for i ∈ 1:2]
    k¹,k²=k=[M.bsplinespaces[i].knots for i ∈ 1:2]

    n₁,n₂=size(𝒂)[1:d]
    𝒂′=Positioning(𝒂)
    return BSplineManifold(𝒫s,𝒂′)
end

function InitBs(D,n₁;nip=NIP)::BSplineManifold
    D₁,D₂=D

    g′₍₀₎(u)=reshape(ForwardDiff.jacobian(g₍₀₎,u),d,d,d) # 第一基本量の微分
    c(t)=[t,sum(extrema(D[2]))/2] # 中心線に沿った座標
    ṡ₍₀₎(t)=sqrt(g₍₀₎(c(t))[1,1])
    s̈₍₀₎(t)=ForwardDiff.derivative(ṡ₍₀₎,t)
    𝛤₍₀₎²₁₁(u)=(g⁻₍₀₎(u)[2,1]*g′₍₀₎(u)[1,1,1]+g⁻₍₀₎(u)[2,2]*(2g′₍₀₎(u)[2,1,1]-g′₍₀₎(u)[1,1,2]))/2 # Christoffel記号
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

function N′(M::BSplineManifold,I₁,I₂,i,u)
    P₁,P₂=P=M.bsplinespaces

    if (i==1)
        return BSplineBasis′(I₁,P₁,u[1])*BSplineBasis(I₂,P₂,u[2])
    else
        return BSplineBasis(I₁,P₁,u[1])*BSplineBasis′(I₂,P₂,u[2])
    end
end

function C(i,j,k,l,g⁻)
    return 𝝀*g⁻[i,j]*g⁻[k,l]+𝝁*(g⁻[i,k]*g⁻[j,l]+g⁻[i,l]*g⁻[j,k])
end

function elm_H(g₍₀₎,M::BSplineManifold,I₁,I₂,i,R₁,R₂,r;nip=NIP)
    𝒂=M.controlpoints
    P₁,P₂=P=M.bsplinespaces
    p₁,p₂=p=P₁.degree,P₂.degree
    k₁,k₂=k=P₁.knots,P₂.knots
    D₁,D₂=D=k₁[1+p₁]..k₁[end-p₁],k₂[1+p₂]..k₂[end-p₂]
    n₁,n₂=n=dim.(P)

    𝜹=[1.0 0.0;0.0 1.0]
    Σ₁=(maximum([I₁,R₁]):minimum([I₁,R₁])+p₁)
    Σ₂=(maximum([I₂,R₂]):minimum([I₂,R₂])+p₂)

    if (length(Σ₁)==0 || length(Σ₂)==0)
        return 0.0
    else
        return sum(GaussianQuadrature(
            u->(
                g=g₍₀₎(u);
                g⁻=inv(g);
                𝝊=sqrt(det(g));
                𝑁=[N′(M,I₁,I₂,i,u) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d];
                Q=[sum(𝒂[I₁,I₂,i]*𝑁[I₁,I₂,j] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂) for i ∈ 1:d, j ∈ 1:d];
                sum(
                    C(p,q,m,n,g⁻)*𝑁[I₁,I₂,p]*(𝜹[i,r]*𝑁[R₁,R₂,q]*(sum(Q[o,m]*Q[o,n] for o ∈ 1:d)-g[m,n])+2*𝑁[R₁,R₂,n]*Q[i,q]*Q[r,m])
                for p ∈ 1:d, q ∈ 1:d, m ∈ 1:d, n ∈ 1:d)
            )*𝝊, k₁[ι₁]..k₁[ι₁+1], k₂[ι₂]..k₂[ι₂+1], nip=nip
        ) for ι₁ ∈ Σ₁, ι₂ ∈ Σ₂)
    end
end


function elm_F(g₍₀₎,M::BSplineManifold,I₁,I₂,i;nip=NIP)
    𝒂=M.controlpoints
    P₁,P₂=P=M.bsplinespaces
    p₁,p₂=p=P₁.degree,P₂.degree
    k₁,k₂=k=P₁.knots,P₂.knots
    D₁,D₂=D=k₁[1+p₁]..k₁[end-p₁],k₂[1+p₂]..k₂[end-p₂]
    n₁,n₂=n=dim.(P)

    D̂₁=BSplineSupport(I₁,P₁)
    D̂₂=BSplineSupport(I₂,P₂)
    Σ₁=(I₁:I₁+p₁)
    Σ₂=(I₂:I₂+p₂)

    return sum(GaussianQuadrature(
        u->(
            g=g₍₀₎(u);
            g⁻=inv(g);
            𝝊=sqrt(det(g));
            𝑁=[N′(M,I₁,I₂,i,u) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d];
            Q=[sum(𝒂[I₁,I₂,i]*𝑁[I₁,I₂,j] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂) for i ∈ 1:d, j ∈ 1:d];
            sum(
                sum(
                    C(p,q,m,n,g⁻)*𝑁[I₁,I₂,p]*Q[i,q]
                    for p ∈ 1:d, q ∈ 1:d
                )*(sum(
                    Q[o,m]*Q[o,n]
                for o ∈ 1:d)-g[m,n])
            for m ∈ 1:d, n ∈ 1:d)
        )*𝝊,k₁[ι₁]..k₁[ι₁+1], k₂[ι₂]..k₂[ι₂+1],nip=nip
    ) for ι₁ ∈ Σ₁, ι₂ ∈ Σ₂)
end

function lineup(n,I₁,I₂,i)
    n₁,n₂=n
    return (i-1)*n₁*n₂+(I₂-1)*n₁+(I₁-1)+1
end

function NewtonIteration(M::BSplineManifold,fixed;nip=NIP)
    𝒂=M.controlpoints
    P₁,P₂=P=M.bsplinespaces
    p₁,p₂=p=P₁.degree,P₂.degree
    k₁,k₂=k=P₁.knots,P₂.knots
    D₁,D₂=D=k₁[1+p₁]..k₁[end-p₁],k₂[1+p₂]..k₂[end-p₂]
    n₁,n₂=n=dim.(P)

    t₀=time()
    Ff=Array{Any}(undef,n₁,n₂,d)
    for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d
        Ff[I₁,I₂,i]=@spawn elm_F(g₍₀₎,M,I₁,I₂,i,nip=nip)
    end
    F=fetch.(Ff)
    Hf=Array{Any}(undef,n₁,n₂,d,n₁,n₂,d)
    for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d, R₁ ∈ 1:n₁, R₂ ∈ 1:n₂, r ∈ 1:d
        Hf[I₁,I₂,i,R₁,R₂,r]=@spawn elm_H(g₍₀₎,M,I₁,I₂,i,R₁,R₂,r,nip=nip)
    end
    H=fetch.(Hf)
    t₁=time()

    # H=[elm_H(g₍₀₎,M,I₁,I₂,i,R₁,R₂,r,nip=nip) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d, R₁ ∈ 1:n₁, R₂ ∈ 1:n₂, r ∈ 1:d]
    # F=[elm_F(g₍₀₎,M,I₁,I₂,i,nip=nip) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d]

    𝕟=2n₁*n₂
    Fixed=sort(collect((i->lineup(n,i...)).(fixed(n₁,n₂))))
    Unfixed=deleteat!(collect(1:𝕟),Fixed)

    F=reshape(F,𝕟)
    H=reshape(H,𝕟,𝕟)
    a=aₒ=reshape(𝒂,𝕟)
    Ȟ=H[Unfixed,Unfixed]
    ǎ=a[Unfixed]
    F̌=F[Unfixed]
    Ǧ=Ȟ\F̌
    ǎ=ǎ-Ǧ
    for i ∈ Fixed
        insert!(ǎ,i,aₒ[i])
    end
    𝒂=reshape(ǎ,n₁,n₂,d)
    M=BSplineManifold(P,𝒂)
    return (M,F,Ǧ,t₁-t₀)
end

mutable struct TreeNode
    parent::Int
    children::Vector{Int}
    comment::String
end
mutable struct Tree
    nodes::Vector{TreeNode}
end
Tree() = Tree([TreeNode(0, Vector{Int}(),"Initial Configuration")])
function addchild(tree::Tree, id::Int, comment::String)
    1 <= id <= length(tree.nodes) || throw(BoundsError(tree, id))
    push!(tree.nodes, TreeNode(id, Vector{}(),comment))
    child = length(tree.nodes)
    push!(tree.nodes[id].children, child)
    child
end
children(tree, id) = tree.nodes[id].children
parent(tree,id) = tree.nodes[id].parent
comment(tree,id)=tree.nodes[id].comment
function shownode(tree,id,depth)
    txt="  "^depth*string(id)*": "*comment(tree,id)*"\n"
    for node ∈ children(tree,id)
        txt=txt*shownode(tree,node,depth+1)
    end
    return txt
end
function showtree(tree)
    shownode(tree,1,0)
end

function Settings(name;up=5,down=-5,right=5,left=-5,mesh=(10,1),unit=100,slack=true,maximumstrain=0.0)
    global NAME=name
    global DIR=OUT_DIR*NAME
    global UP=up
    global DOWN=down
    global RIGHT=right
    global LEFT=left
    global MESH=mesh
    global UNIT=(unit,"pt")
    global SLACK=slack
    global MAXIMUMSTRAIN=maximumstrain
    return nothing
end

function Export(M::BSplineManifold,BsTree,BsJLD;comment="",maximumstrain=MAXIMUMSTRAIN)
    index=length(BsTree.nodes)
    BsJLD[string(index)]=M
    BsJLD["BsTree"]=BsTree
    save(DIR*"/"*NAME*".jld",BsJLD)
    println(showtree(BsTree))

    Name=NAME
    Dir=DIR
    Up=UP
    Down=DOWN
    Right=RIGHT
    Left=LEFT
    Mesh=MESH
    Unit=UNIT
    Slack=SLACK
    if (maximumstrain==0.0)
        MS=ComputeMaximumStrain(index=index)
        MaximumStrain=max(-MS[1],MS[2])
    else
        MaximumStrain=maximumstrain
    end

    @spawnat 1 begin
        BSplineSvg(M,filename=Dir*"/nurbs/"*Name*"-"*string(index)*"_Bspline.svg",up=Up,down=Down,right=Right,left=Left,mesh=Mesh,unitlength=Unit)
        𝒂 = M.controlpoints
        P₁,P₂=P=M.bsplinespaces
        p₁,p₂=p=P₁.degree,P₂.degree
        k₁,k₂=k=P₁.knots,P₂.knots
        D₁,D₂=D=k₁[1+p₁]..k₁[end-p₁],k₂[1+p₂]..k₂[end-p₂]
        n₁,n₂=n=dim.(P)

        𝒑₍ₜ₎(u)=Mapping(M,u)
        𝒑₁₍ₜ₎(u) = sum(BSplineBasis′(I₁,P₁,u[1])*BSplineBasis(I₂,P₂,u[2])*𝒂[I₁,I₂,:] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂)
        g₍₀₎₁₁(u)=dot(𝒑₁₍₀₎(u),𝒑₁₍₀₎(u))
        g₍ₜ₎₁₁(u)=dot(𝒑₁₍ₜ₎(u),𝒑₁₍ₜ₎(u))
        E₁₁(u)=(g₍ₜ₎₁₁(u)-g₍₀₎₁₁(u))/2
        E⁽⁰⁾₁₁(u)=E₁₁(u)/g₍₀₎₁₁(u)

        rgb(u)=E⁽⁰⁾₁₁(u)*[1,-1,-1]/(2*MaximumStrain) .+0.5
        # draw strain distribution (6000x6000)
        ParametricColor(𝒑₍ₜ₎,D,rgb=rgb,filename=Dir*"/strain/"*Name*"-"*string(index)*"_strain.png",up=Up,down=Down,right=Right,left=Left,mesh=tuple(10*[Mesh...]...),unit=5*Unit[1])
        ColorBar(max=MaximumStrain,filename=Dir*"/colorbar/"*Name*"-"*string(index)*"_colorbar.png",width=(Right-Left)*Unit[1])

        # 1200x1200
        # svg to png (1600x1600)
        run(pipeline(`convert $(Dir*"/nurbs/"*Name*"-"*string(index)*"_Bspline.svg") $(Dir*"/nurbs/"*Name*"-"*string(index)*"_Bspline.png")`, stdout=devnull, stderr=devnull))
        # add colorbar to strain distribution figure (6000x6000)
        run(pipeline(`convert $(Dir*"/strain/"*Name*"-"*string(index)*"_strain.png") $(Dir*"/colorbar/"*Name*"-"*string(index)*"_colorbar.png") -gravity southeast -compose over -composite $(Dir*"/strain/"*Name*"-"*string(index)*"_swc.png")`, stdout=devnull, stderr=devnull))
        # resize png
        # (1200x1200)
        run(pipeline(`convert -resize 75% -unsharp 2x1.4+0.5+0 -quality 100 -verbose $(Dir*"/nurbs/"*Name*"-"*string(index)*"_Bspline.png") $(Dir*"/slack/"*Name*"-"*string(index)*"_Bspline.png")`, stdout=devnull, stderr=devnull))
        # (1200x1200)
        run(pipeline(`convert -resize 20% -unsharp 2x1.4+0.5+0 -quality 100 -verbose $(Dir*"/strain/"*Name*"-"*string(index)*"_swc.png") $(Dir*"/slack/"*Name*"-"*string(index)*"_strain.png")`, stdout=devnull, stderr=devnull))
        # line up png
        run(pipeline(`convert +append $(Dir*"/slack/"*Name*"-"*string(index)*"_Bspline.png") $(Dir*"/slack/"*Name*"-"*string(index)*"_strain.png") $(Dir*"/slack/"*Name*"-"*string(index)*"_append.png")`, stdout=devnull, stderr=devnull))

        if (Slack)
            SlackString(showtree(BsTree))
            SlackFile(Dir*"/slack/"*Name*"-"*string(index)*"_append.png")
        end
    end

    return nothing
end

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

function BSpline.Refinement(;p₊::Union{Nothing,Array{Int,1}}=nothing, k₊::Union{Nothing,Array{Knots,1}}=nothing, parent=0)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    M=BsJLD[string(parent)]

    comment="refinement with "*string(k₊)
    addchild(BsTree,parent,comment)

    M=BSpline.Refinement(M,p₊=p₊,k₊=k₊)
    Export(M,BsTree,BsJLD,comment=comment)
end

function NewtonMethodIteration(;fixed=((n₁,n₂)->([(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1])),parent=0,nip=NIP)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    M=BsJLD[string(parent)]
    n₁,n₂=n=dim.(M.bsplinespaces)
    if (!isodd(n₁*n₂)) error("n₁ and n₂ should be odd numbers") end
    M=Positioning(M)
    M,F,Ǧ,Δt=NewtonIteration(M,fixed,nip=nip)
    # comment="Newton Iteration - residual norm: "*(@sprintf("%.5e",norm(F)))*", Δa norm: "*(@sprintf("%.5e",norm(Ǧ)))*", computation time: "*(@sprintf("%.5e",Δt))*" sec"
    comment="Newton Iteration - residual norm: "*(@sprintf("%.5e",norm(F)))*", Δa norm: "*(@sprintf("%.5e",norm(Ǧ)))*", computation time: "*string(Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(1000Δt÷1))))
    addchild(BsTree,parent,comment)

    Export(M,BsTree,BsJLD,comment=comment)
end

function Restoration()
    if (!isfile(DIR*"/"*NAME*".jld"))
        error("jld file doesn't exists")
    end
    BsJLD=load(DIR*"/"*NAME*".jld")
    println(showtree(BsJLD["BsTree"]))
    global EXPR=BsJLD["Expr"]
    eval(:(@everywhere $EXPR))
    return nothing
end

function FinalOutput(;index=0,unitlength=(10,"mm"),cutout=(0.1,5),mesh=60)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    # println(showtree(BsTree))

    if (index==0) index=length(BsTree.nodes) end
    M=BsJLD[string(index)]
    BSplineSvg(M,filename=DIR*"/"*NAME*"-"*string(index)*"-final.svg",up=UP,down=DOWN,right=RIGHT,left=LEFT,mesh=MESH,unitlength=unitlength,points=false)

    P₁,P₂=P=M.bsplinespaces
    p₁,p₂=p=P₁.degree,P₂.degree
    k₁,k₂=k=P₁.knots,P₂.knots
    D₁,D₂=D=k₁[1+p₁]..k₁[end-p₁],k₂[1+p₂]..k₂[end-p₂]
    n₁,n₂=n=dim.(P)

    𝒑₍ₜ₎(u)=Mapping(M,u)
    𝒑₁₍ₜ₎(u) = sum(BSplineBasis′(I₁,P₁,u[1])*BSplineBasis(I₂,P₂,u[2])*𝒂[I₁,I₂,:] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂)
    𝒆⁽⁰⁾₁(u)=normalize(𝒑₁₍ₜ₎(u))
    𝒆⁽⁰⁾₂(u)=[0.0 -1.0;1.0 0.0]*𝒆⁽⁰⁾₁(u)
    𝒑a(i,t)=𝒑₍ₜ₎([t,leftendpoint(D₂)])+𝒆⁽⁰⁾₂([t,leftendpoint(D₂)])*i*cutout[1]/unitlength[1]
    𝒑b(i,t)=𝒑₍ₜ₎([t,rightendpoint(D₂)])-𝒆⁽⁰⁾₂([t,rightendpoint(D₂)])*i*cutout[1]/unitlength[1]
    SvgCurve([[t->𝒑a(i,t) for i ∈ 0:cutout[2]]...,[t->𝒑b(i,t) for i ∈ 0:cutout[2]]...],D₁,filename=DIR*"/"*NAME*"-"*string(index)*"-cutout.svg",up=UP,down=DOWN,right=RIGHT,left=LEFT,thickness=0.1,mesh=mesh,unitlength=unitlength)

    if (SLACK)
        SlackFile(DIR*"/"*NAME*"-"*string(index)*"-final.svg")
        SlackFile(DIR*"/"*NAME*"-"*string(index)*"-cutout.svg")
    end

    return nothing
end

function ShowKnots(;index=0)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (index==0) index=length(BsTree.nodes) end
    M=BsJLD[string(index)]
    P₁,P₂=M.bsplinespaces
    p₁,p₂=P₁.degree,P₂.degree
    k₁,k₂=P₁.knots,P₂.knots
    println("k₁: ",k₁)
    println("k₂: ",k₂)
    println("Suggestion:")
    k₁′=unique(k₁)
    k₂′=unique(k₂)
    println("k₁₊: ",[(k₁′[i]+k₁′[i+1])/2 for i ∈ 1:(length(k₁′)-1)])
    println("k₂₊: ",[(k₂′[i]+k₂′[i+1])/2 for i ∈ 1:(length(k₂′)-1)])
    return nothing
end

function ComputeMaximumStrain(;index=0,mesh=tuple(20*[MESH...]...))
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (index==0) index=length(BsTree.nodes) end
    M=BsJLD[string(index)]
    𝒂=M.controlpoints
    P₁,P₂=P=M.bsplinespaces
    p₁,p₂=p=P₁.degree,P₂.degree
    k₁,k₂=k=P₁.knots,P₂.knots
    D₁,D₂=D=k₁[1+p₁]..k₁[end-p₁],k₂[1+p₂]..k₂[end-p₂]
    n₁,n₂=n=dim.(P)

    𝒑₍ₜ₎(u)=Mapping(M,u)
    𝒑₁₍ₜ₎(u) = sum(BSplineBasis′(I₁,P₁,u[1])*BSplineBasis(I₂,P₂,u[2])*𝒂[I₁,I₂,:] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂)
    g₍₀₎₁₁(u)=dot(𝒑₁₍₀₎(u),𝒑₁₍₀₎(u))
    g₍ₜ₎₁₁(u)=dot(𝒑₁₍ₜ₎(u),𝒑₁₍ₜ₎(u))
    E₁₁(u)=(g₍ₜ₎₁₁(u)-g₍₀₎₁₁(u))/2
    E⁽⁰⁾₁₁(u)=E₁₁(u)/g₍₀₎₁₁(u)

    κ₁=range(leftendpoint(D₁),stop=rightendpoint(D₁),length=mesh[1]+1)
    κ₂=range(leftendpoint(D₂),stop=rightendpoint(D₂),length=mesh[2]+1)

    E=[E⁽⁰⁾₁₁([u₁,u₂]) for u₁ ∈ κ₁, u₂ ∈ κ₂]

    return (minimum(E),maximum(E))
end

function ShowMaximumStrain(;index=0,mesh=5)
    minE,maxE=ComputeMaximumStrain(index=index,mesh=mesh)
    println("min",minE,", max",maxE)

    return nothing
end



function ReDraw()
    return nothing
end


end
