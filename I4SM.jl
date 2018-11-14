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
using SvgDraw
# using POV_Ray

export InitialConfiguration, p_Refinement, h_Refinement, NewtonMethodIteration, FinalOutput, ShowKnots, Settings

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
    ind0=[(n₁+1)÷2,(n₂+1)÷2]
    ind1=ind0-[0,1]
    v=a[ind1...,:]-a[ind0...,:]
    R=-[v[2] -v[1];v[1] v[2]]/norm(v)
    return a=aff(a,R,-R*a[ind0...,:])
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
    m=BsCoef2(𝒄,p₁,k₁,nip=nip)
    m₂=BsCoef2(𝒄₂,p₁,k₁,nip=nip)
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
    txt="  "^depth*string(id)*": "*comment(tree,id)*"\n"
    for node ∈ children(tree,id)
        txt=txt*shownode(tree,node,depth+1)
    end
    return txt
end
function showtree(tree)
    shownode(tree,1,0)
end

function Settings(name;up=5,down=-5,right=5,left=-5,mesh=(10,1),unit=100,slack=true,maximumstrain=0.01)
    global NAME=name
    global DIR=homedir()*"/I4SM-Result/"*NAME
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

function Export(𝒑₍₀₎,B2::Bs2mfd,BsTree,BsJLD;comment="",maximumstrain=MAXIMUMSTRAIN)
    index=length(BsTree.nodes)
    BsJLD[string(index)]=B2
    BsJLD["BsTree"]=BsTree
    save(DIR*"/"*NAME*".jld",BsJLD)
    println(showtree(BsTree))
    BsDraw(B2,filename=DIR*"/svg/"*NAME*"-"*string(index)*"-Bspline.svg",up=UP,down=DOWN,right=RIGHT,left=LEFT,mesh=MESH,unitlength=UNIT)
    k₁,k₂=B2.k
    D=(k₁[1]..k₁[end],k₂[1]..k₂[end])

    𝒑₁₍₀₎(u)=ForwardDiff.derivative(u₁->𝒑₍₀₎([u₁,u[2]]),u[1])
    𝒑₍ₜ₎(u)=BsMapping(B2,u)
    function 𝒑₁₍ₜ₎(u)
        p,k,a=B2.p,B2.k,B2.a
        p₁,p₂=p
        k₁,k₂=k
        n=n₁,n₂=length.(k)-p.-1
        return sum(Ḃs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])*a[I₁,I₂,:] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂)
    end
    g₍₀₎₁₁(u)=dot(𝒑₁₍₀₎(u),𝒑₁₍₀₎(u))
    g₍ₜ₎₁₁(u)=dot(𝒑₁₍ₜ₎(u),𝒑₁₍ₜ₎(u))
    E₁₁(u)=(g₍ₜ₎₁₁(u)-g₍₀₎₁₁(u))/2
    E⁽⁰⁾₁₁(u)=E₁₁(u)/g₍₀₎₁₁(u)

    rgb(u)=E⁽⁰⁾₁₁(u)*[1,-1,-1]/(2*maximumstrain) .+0.5
    ParametricColor(𝒑₍ₜ₎,D,rgb=rgb,filename=DIR*"/strain/"*NAME*"-"*string(index)*"-strain.png",up=UP,down=DOWN,right=RIGHT,left=LEFT,mesh=tuple(10*[MESH...]...),unit=5*UNIT[1])
    ColorBar(max=maximumstrain,filename=DIR*"/colorbar/"*NAME*"-"*string(index)*"-colorbar.png",unit=UNIT[1])

    run(`convert $(DIR*"/svg/"*NAME*"-"*string(index)*"-Bspline.svg") $(DIR*"/slack/"*NAME*"-"*string(index)*"-Bspline.png")`)
    run(`convert -resize 80% -unsharp 2x1.4+0.5+0 -quality 100 -verbose $(DIR*"/slack/"*NAME*"-"*string(index)*"-Bspline.png") $(DIR*"/slack/"*NAME*"-"*string(index)*"-Bspline.png")`)
    run(`convert $(DIR*"/strain/"*NAME*"-"*string(index)*"-strain.png") $(DIR*"/colorbar/"*NAME*"-"*string(index)*"-colorbar.png") -gravity southeast -compose over -composite $(DIR*"/slack/"*NAME*"-"*string(index)*"-strain.png")`)
    run(`convert -resize 20% -unsharp 2x1.4+0.5+0 -quality 100 -verbose $(DIR*"/slack/"*NAME*"-"*string(index)*"-strain.png") $(DIR*"/slack/"*NAME*"-"*string(index)*"-strain.png")`)
    run(`convert +append $(DIR*"/slack/"*NAME*"-"*string(index)*"-Bspline.png") $(DIR*"/slack/"*NAME*"-"*string(index)*"-strain.png") $(DIR*"/slack/"*NAME*"-"*string(index)*"-append.png")`)

    if (SLACK)
        SlackString(showtree(BsTree))
        SlackFile(DIR*"/slack/"*NAME*"-"*string(index)*"-append.png")
    end
    return nothing
end

function InitialConfiguration(𝒑₍₀₎,D;n₁=15,nip=25)
    if (isfile(DIR*"/"*NAME*".jld")) error("File already exists") end
    mkpath(DIR)
    mkpath(DIR*"/svg")
    mkpath(DIR*"/strain")
    mkpath(DIR*"/colorbar")
    mkpath(DIR*"/slack")
    BsJLD=Dict{String,Any}()

    B2=InitBs(𝒑₍₀₎,D,n₁,nip=25)
    comment="Initial Configuration"
    BsTree=Tree()

    Export(𝒑₍₀₎,B2,BsTree,BsJLD,comment=comment)
end

function p_Refinement(𝒑₍₀₎,p₊::Array{Int64,1};parent=0,nip=25)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    B2=BsJLD[string(parent)]

    B2=pref(B2,p₊,nip=nip)
    comment="p-refinement with "*string(p₊)
    addchild(BsTree,parent,comment)

    Export(𝒑₍₀₎,B2,BsTree,BsJLD,comment=comment)
end

function h_Refinement(𝒑₍₀₎,h₊::Array{Array{Float64,1},1};parent=0,nip=25)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    B2=BsJLD[string(parent)]

    B2=href(B2,h₊,nip=nip)
    comment="h-refinement with "*string(h₊)
    addchild(BsTree,parent,comment)

    Export(𝒑₍₀₎,B2,BsTree,BsJLD,comment=comment)
end

function NewtonMethodIteration(𝒑₍₀₎;fixed=((n₁,n₂)->([(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1])),parent=0,nip=25)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    B2=BsJLD[string(parent)]

    n₁,n₂=length.(B2.k)-B2.p.-1
    if (!isodd(n₁*n₂)) error("n₁ and n₂ should be odd numbers") end
    B2=Positioning(B2)
    B2,F,Ǧ=NewtonIteration(𝒑₍₀₎,B2,fixed,nip=nip)
    comment="Newton Iteration - Residual norm: "*string(norm(F))*", Δa norm: "*string(norm(Ǧ))
    addchild(BsTree,parent,comment)

    Export(𝒑₍₀₎,B2,BsTree,BsJLD,comment=comment)
end

function FinalOutput(;index=0,unitlength=(10,"mm"),cutout=(0.1,5),mesh=60)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    # println(showtree(BsTree))

    if (index==0) index=length(BsTree.nodes) end
    B2=BsJLD[string(index)]
    BsDraw(B2,filename=DIR*"/"*NAME*"-"*string(index)*"-final.svg",up=UP,down=DOWN,right=RIGHT,left=LEFT,mesh=MESH,unitlength=unitlength,points=false)

    k₁,k₂=B2.k
    D₁=ClosedInterval(k₁[1],k₁[end])
    D₂=ClosedInterval(k₂[1],k₂[end])

    𝒑₍ₜ₎(u)=BsMapping(B2,u)
    function 𝒑₁₍ₜ₎(u)
        p,k,a=B2.p,B2.k,B2.a
        p₁,p₂=p
        k₁,k₂=k
        n=n₁,n₂=length.(k)-p.-1
        return sum(Ḃs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])*a[I₁,I₂,:] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂)
    end
    𝒆⁽⁰⁾₁(u)=𝒑₁₍ₜ₎(u)
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
    B2=BsJLD[string(index)]
    k₁,k₂=B2.k
    println("k₁: ",k₁)
    println("k₂: ",k₂)
    println("Suggestion:")
    k₁′=DelDpl(k₁)
    k₂′=DelDpl(k₂)
    println("k₁₊: ",[(k₁′[i]+k₁′[i+1])/2 for i ∈ 1:(length(k₁′)-1)])
    println("k₂₊: ",[(k₂′[i]+k₂′[i+1])/2 for i ∈ 1:(length(k₂′)-1)])
    return nothing
end


end
