module I4SM

using LinearAlgebra
using Distributed
using IntervalSets
using ForwardDiff
using DifferentialEquations
using JLD

using Bspline
using ElementaryCalculus
using Slack
# using SvgDraw
# using POV_Ray

export Init, pRef, hRef, Newt, DrawConfig

const 𝝂=0.4 #Poisson比ν
const d=2
const 𝝀=𝝂/((1+𝝂)*(1-(d-1)*𝝂)) #Lamé定数λ
const 𝝁=1/(2+2𝝂) #Lamé定数μ

function aff(a::Array{Float64,3},A::Array{Float64,2},b::Array{Float64,1})
    #x'=Ax+b
    n₁,n₂=size(a)[1:2]
    return [sum(A[i,j]*a[I₁,I₂,j] for j ∈ 1:2)+b[i] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:2]
end

function Positioning(a::Array{Float64,3}) # 制御点の位置調整
    n₁,n₂=size(a)[1:2]
    v=a[(n₁+1)÷2,(n₂-1)÷2,:]-a[(n₁+1)÷2,(n₂+1)÷2,:]
    R=-[v[2] -v[1];v[1] v[2]]/norm(v)
    return a=aff(a,R,-R*a[(n₁+1)÷2,(n₂+1)÷2,:])
end

function Positioning(B2::Bs2mfd) # 制御点の位置調整
    p,k,a=B2.p,B2.k,B2.a
    n₁,n₂=size(a)[1:2]
    aa=Positioning(a)
    return Bs2mfd(p,k,aa)
end

function InitBs(𝒑₍₀₎,D,n₁;nip=25)
    D₁,D₂=D
    𝒑′₍₀₎(u)=ForwardDiff.jacobian(𝒑₍₀₎,u) # 接ベクトル
    𝒑₁₍₀₎(u)=ForwardDiff.derivative(u₁->𝒑₍₀₎([u₁,u[2]]),u[1])
    𝒑₂₍₀₎(u)=ForwardDiff.derivative(u₂->𝒑₍₀₎([u[1],u₂]),u[2])
    g₍₀₎(u)=𝒑′₍₀₎(u)'𝒑′₍₀₎(u) # 第一基本量
    𝝊₍₀₎(u)=norm(cross(𝒑₁₍₀₎(u),𝒑₂₍₀₎(u))) # 体積要素υ
    g⁻₍₀₎(u)=inv(g₍₀₎(u)) # 第一基本量の逆

    g′₍₀₎(u)=reshape(ForwardDiff.jacobian(g₍₀₎,u),2,2,2) # 第一基本量の微分
    c(t)=[t,sum(extrema(D[2]))/2] #中心線に沿った座標
    ṡ₍₀₎(t)=sqrt(g₍₀₎(c(t))[1,1])
    s̈₍₀₎(t)=ForwardDiff.derivative(ṡ₍₀₎,t)
    𝛤₍₀₎²₁₁(u)=(g⁻₍₀₎(u)[2,1]*g′₍₀₎(u)[1,1,1]+g⁻₍₀₎(u)[2,2]*(2g′₍₀₎(u)[2,1,1]-g′₍₀₎(u)[1,1,2]))/2 #Christoffel記号
    𝜅₍₀₎(t)=𝛤₍₀₎²₁₁(c(t))*𝝊₍₀₎(c(t))/ṡ₍₀₎(t)^3 # 測地的曲率

    function ode(𝒄̇𝒄̈,𝒄𝒄̇,par,t)
        𝒄̇𝒄̈[1]=𝒄𝒄̇[3]
        𝒄̇𝒄̈[2]=𝒄𝒄̇[4]
        𝒄̇𝒄̈[3]=dot([s̈₍₀₎(t)/ṡ₍₀₎(t),-𝜅₍₀₎(t)*ṡ₍₀₎(t)],𝒄𝒄̇[3:4])
        𝒄̇𝒄̈[4]=dot([𝜅₍₀₎(t)*ṡ₍₀₎(t),s̈₍₀₎(t)/ṡ₍₀₎(t)],𝒄𝒄̇[3:4])
    end
    𝒄𝒄̇₀=vcat([0.0,0.0],[1.,0.]*ṡ₍₀₎(minimum(D₁)))
    sol=solve(ODEProblem(ode,𝒄𝒄̇₀,extrema(D₁)))
    𝒄(t)=sol(t)[1:2] # 解となる中心曲線
    𝒄₁(t)=sol(t)[3:4] # その導関数
    𝒄₂(t)=[g₍₀₎(c(t))[1,2] -𝝊₍₀₎(c(t));𝝊₍₀₎(c(t)) g₍₀₎(c(t))[1,2]]*𝒄₁(t)/g₍₀₎(c(t))[1,1] # 中心曲線上の幅方向のベクトル場

    p₁=3
    k₁=sort(vcat(repeat(collect(extrema(D₁)),inner=p₁),collect(range(leftendpoint(D₁),stop=rightendpoint(D₁),length=n₁-2))))
    m=BsCoef(𝒄,p₁,k₁,nip=nip)
    m₂=BsCoef(𝒄₂,p₁,k₁,nip=nip)
    a1=m-width(D₂)*m₂/2
    a2=m+width(D₂)*m₂/2
    p₂=1
    k₂=repeat(collect(extrema(D₂)),inner=2)
    n₂=length(k₂)-p₂-1
    p=[p₁,p₂]
    k=[k₁,k₂]
    a=[[a1[I₁,i],a2[I₁,i]][I₂] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:2]

    return Positioning(pref(Bs2mfd(p,k,a),[0,1]))
end

function N′(B2::Bs2mfd,I₁,I₂,i,u)
    p,k,a=B2.p,B2.k,B2.a
    p₁,p₂=p
    k₁,k₂=k
    if(i==1)
        return Ḃs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])
    else
        return Bs(I₁,p₁,k₁,u[1])*Ḃs(I₂,p₂,k₂,u[2])
    end
end

function C(i,j,k,l,g⁻)
    return 𝝀*g⁻[i,j]*g⁻[k,l]+𝝁*(g⁻[i,k]*g⁻[j,l]+g⁻[i,l]*g⁻[j,k])
end

function elm_H(g₍₀₎,B2::Bs2mfd,I₁,I₂,i,R₁,R₂,r;nip=25)
    p,k,a=B2.p,B2.k,B2.a
    p₁,p₂=p
    k₁,k₂=k
    n₁,n₂=length.(k)-p.-1
    D̂₁=Bsupp(I₁,p₁,k₁)∩Bsupp(R₁,p₁,k₁)
    D̂₂=Bsupp(I₂,p₂,k₂)∩Bsupp(R₂,p₂,k₂)
    𝜹=[1.0 0.0;0.0 1.0]
    if(isnullset(D̂₁)||isnullset(D̂₂))
        return 0.0
    else
        return INT2(
            u->(
                g=g₍₀₎(u);
                g⁻=inv(g);
                𝝊=sqrt(det(g));
                𝑁=[N′(B2,I₁,I₂,i,u) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d];
                an=[sum(a[I₁,I₂,i]*𝑁[I₁,I₂,j] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂) for i ∈ 1:2, j ∈ 1:2];
                sum(
                    C(p,q,m,n,g⁻)*𝑁[I₁,I₂,p]*(𝜹[i,r]*𝑁[R₁,R₂,q]*(sum(an[o,m]*an[o,n] for o ∈ 1:2)-g[m,n])+2*𝑁[R₁,R₂,n]*an[i,q]*an[r,m])
                for p ∈ 1:2, q ∈ 1:2, m ∈ 1:2, n ∈ 1:2)
            )*𝝊,(D̂₁,D̂₂),nip=nip
        )
    end
end

function elm_F(g₍₀₎,B2::Bs2mfd,I₁,I₂,i;nip=25)
    p,k,a=B2.p,B2.k,B2.a
    p₁,p₂=p
    k₁,k₂=k
    n₁,n₂=length.(k)-p.-1
    D̂₁=Bsupp(I₁,p₁,k₁)
    D̂₂=Bsupp(I₂,p₂,k₂)
    return INT2(
        u->(
            g=g₍₀₎(u);
            g⁻=inv(g);
            𝝊=sqrt(det(g));
            𝑁=[N′(B2,I₁,I₂,i,u) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:2];
            an=[sum(a[I₁,I₂,i]*𝑁[I₁,I₂,j] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂) for i ∈ 1:2, j ∈ 1:2];
            sum(
                sum(
                    C(p,q,m,n,g⁻)*𝑁[I₁,I₂,p]*an[i,q]
                    for p ∈ 1:2, q ∈ 1:2
                )*(sum(
                    an[o,m]*an[o,n]
                for o ∈ 1:d)-g[m,n])
            for m ∈ 1:2, n ∈ 1:2)
        )*𝝊,(D̂₁,D̂₂),nip=nip
    )
end

function FFd(𝒑₍₀₎,B2::Bs2mfd;nip=25)
    𝒑′₍₀₎(u)=ForwardDiff.jacobian(𝒑₍₀₎,u) # 接ベクトル
    g₍₀₎(u)=𝒑′₍₀₎(u)'𝒑′₍₀₎(u)

    n₁,n₂=n=length.(B2.k)-B2.p.-1

    Ff=Array{Any}(undef,n₁,n₂,2)
    for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d
        Ff[I₁,I₂,i]=@spawn elm_F(g₍₀₎,B2,I₁,I₂,i,nip=nip)
    end
    @time F=fetch.(Ff)

    return F

    # F=[elm_F(g₍₀₎,B2,I₁,I₂,i,nip=nip) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:2]
    # return F
end

function lineup(n,I₁,I₂,i)
    n₁,n₂=n
    return (i-1)*n₁*n₂+(I₂-1)*n₁+(I₁-1)+1
end

function NewtonIteration(𝒑₍₀₎,B2::Bs2mfd,fixed;nip=25)
    𝒑′₍₀₎(u)=ForwardDiff.jacobian(𝒑₍₀₎,u) # 接ベクトル
    g₍₀₎(u)=𝒑′₍₀₎(u)'𝒑′₍₀₎(u)
    n₁,n₂=n=length.(B2.k)-B2.p.-1

    Ff=Array{Any}(undef,n₁,n₂,2)
    for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d
        Ff[I₁,I₂,i]=@spawn elm_F(g₍₀₎,B2,I₁,I₂,i,nip=nip)
    end
    F=fetch.(Ff)
    Hf=Array{Any}(undef,n₁,n₂,2,n₁,n₂,2)
    for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d, R₁ ∈ 1:n₁, R₂ ∈ 1:n₂, r ∈ 1:d
        Hf[I₁,I₂,i,R₁,R₂,r]=@spawn elm_H(g₍₀₎,B2,I₁,I₂,i,R₁,R₂,r,nip=nip)
    end
    H=fetch.(Hf)
    # H=[elm_H(g₍₀₎,B2,I₁,I₂,i,R₁,R₂,r,nip=nip) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:2, R₁ ∈ 1:n₁, R₂ ∈ 1:n₂, r ∈ 1:2]
    # F=[elm_F(g₍₀₎,B2,I₁,I₂,i,nip=nip) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:2]

    𝕟=2n₁*n₂
    Fixed=sort(collect((i->lineup(n,i...)).(fixed(n₁,n₂))))
    Unfixed=deleteat!(collect(1:𝕟),Fixed)

    H=reshape(H,𝕟,𝕟)
    F=reshape(F,𝕟)
    a=aₒ=reshape(B2.a,𝕟)
    Ȟ=H[Unfixed,Unfixed]
    ǎ=a[Unfixed]
    F̌=F[Unfixed]
    Ǧ=Ȟ\F̌
    ǎ=ǎ-Ǧ
    for i ∈ Fixed
        insert!(ǎ,i,aₒ[i])
    end
    a=reshape(ǎ,n₁,n₂,2)
    return (Bs2mfd(B2.p,B2.k,a),F,Ǧ)
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
    println("  "^depth*string(id)*": "*comment(tree,id))
    for node ∈ children(tree,id)
        shownode(tree,node,depth+1)
    end
end
function showtree(tree)
    shownode(tree,1,0)
end


function DrawConfig(config)
    global UP=config["UP"]
    global DOWN=config["DOWN"]
    global RIGHT=config["RIGHT"]
    global LEFT=config["LEFT"]
    global MESH=config["MESH"]
    global UNIT=config["UNIT"]
end

function ExportFig(ShapeName,𝒑₍₀₎,B2::Bs2mfd,index)
    BsDraw(B2,filename=homedir()*"/I4SM-Result/"*ShapeName*"/svg/"*ShapeName*"-"*string(index)*".svg",up=UP,down=DOWN,right=RIGHT,left=LEFT,mesh=MESH,unitlength=UNIT)
    run(`convert $(homedir()*"/I4SM-Result/"*ShapeName*"/svg/"*ShapeName*"-"*string(index)*".svg") $(homedir()*"/I4SM-Result/"*ShapeName*"/slack/"*ShapeName*"-"*string(index)*".png")`)
    SlackFile(homedir()*"/I4SM-Result/"*ShapeName*"/slack/"*ShapeName*"-"*string(index)*".png")
end

function Init(ShapeName,𝒑₍₀₎,D;n₁=15,nip=25)
    mkpath(homedir()*"/I4SM-Result/"*ShapeName)
    mkpath(homedir()*"/I4SM-Result/"*ShapeName*"/svg")
    mkpath(homedir()*"/I4SM-Result/"*ShapeName*"/slack")
    mkpath(homedir()*"/I4SM-Result/"*ShapeName*"/strain")
    B=InitBs(𝒑₍₀₎,D,n₁,nip=25)
    BsTree=Tree()
    BsJLD=Dict{String,Any}()

    index=1
    BsJLD[string(index)]=B
    BsJLD["BsTree"]=BsTree
    save(homedir()*"/I4SM-Result/"*ShapeName*"/"*ShapeName*".jld",BsJLD)
    showtree(BsTree)
    ExportFig(ShapeName,𝒑₍₀₎,B,index)
    return nothing
end
function pRef(ShapeName,𝒑₍₀₎,p₊::Array{Int64,1};parent=0,nip=25)
    BsJLD=load(homedir()*"/I4SM-Result/"*ShapeName*"/"*ShapeName*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    B=BsJLD[string(parent)]

    B=pref(B,p₊,nip=nip)
    addchild(BsTree,parent,"p-refinement with "*string(p₊))

    index=length(BsTree.nodes)
    BsJLD[string(index)]=B
    BsJLD["BsTree"]=BsTree
    save(homedir()*"/I4SM-Result/"*ShapeName*"/"*ShapeName*".jld",BsJLD)
    showtree(BsTree)
    ExportFig(ShapeName,𝒑₍₀₎,B,index)
    return nothing
end
function hRef(ShapeName,𝒑₍₀₎,h₊::Array{Array{Float64,1},1};parent=0,nip=25)
    BsJLD=load(homedir()*"/I4SM-Result/"*ShapeName*"/"*ShapeName*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    B=BsJLD[string(parent)]

    B=href(B,h₊,nip=nip)
    addchild(BsTree,parent,"h-refinement with "*string(h₊))

    index=length(BsTree.nodes)
    BsJLD[string(index)]=B
    BsJLD["BsTree"]=BsTree
    save(homedir()*"/I4SM-Result/"*ShapeName*"/"*ShapeName*".jld",BsJLD)
    showtree(BsTree)
    ExportFig(ShapeName,𝒑₍₀₎,B,index)
    return nothing
end
function Newt(ShapeName,𝒑₍₀₎;fixed=((n₁,n₂)->([(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1])),parent=0,nip=25)
    BsJLD=load(homedir()*"/I4SM-Result/"*ShapeName*"/"*ShapeName*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    B=BsJLD[string(parent)]

    B,F,Ǧ=NewtonIteration(𝒑₍₀₎,B,fixed,nip=nip)
    addchild(BsTree,parent,"Newton Iteration - Residual norm: "*string(norm(F))*", Δa norm: "*string(norm(Ǧ)))

    index=length(BsTree.nodes)
    BsJLD[string(index)]=B
    BsJLD["BsTree"]=BsTree
    save(homedir()*"/I4SM-Result/"*ShapeName*"/"*ShapeName*".jld",BsJLD)
    showtree(BsTree)
    ExportFig(ShapeName,𝒑₍₀₎,B,index)
    return nothing
end

end
