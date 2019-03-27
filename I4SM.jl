module I4SM

using LinearAlgebra
using Printf
using Distributed
using IntervalSets
using ForwardDiff
using Dates
using DifferentialEquations
using JLD

using Bspline
using ElementaryCalculus
using Slack
using ParametricDraw
# using POV_Ray

export @DefineShape, InitialConfiguration, p_Refinement, h_Refinement, NewtonMethodIteration, FinalOutput, ShowKnots, ShowMaximumStrain, Settings, Restoration

const d=2 #Dimension
const 𝝂=0.25 #Poisson比ν
const Y=1.0 #Young率Y
const 𝝀=𝝂*Y/((1+𝝂)*(1-(d-1)*𝝂)) #Lamé定数λ
const 𝝁=1/2(1+𝝂) #Lamé定数μ
const NIP=25 # Number of Integration Points

macro DefineShape(ex)
    global EXPR=ex
    return :(@everywhere $ex)
end

𝒑′₍₀₎(u)=ForwardDiff.jacobian(Main.𝒑₍₀₎,u) # 接ベクトル
𝒑₁₍₀₎(u)=ForwardDiff.derivative(u₁->Main.𝒑₍₀₎([u₁,u[2]]),u[1])
𝒑₂₍₀₎(u)=ForwardDiff.derivative(u₂->Main.𝒑₍₀₎([u[1],u₂]),u[2])
g₍₀₎(u)=𝒑′₍₀₎(u)'𝒑′₍₀₎(u) # 第一基本量
𝝊₍₀₎(u)=norm(cross(𝒑₁₍₀₎(u),𝒑₂₍₀₎(u))) # 体積要素υ
g⁻₍₀₎(u)=inv(g₍₀₎(u)) # 第一基本量の逆

function aff(a::Array{Float64,3},A::Array{Float64,2},b::Array{Float64,1})
    #x'=Ax+b
    n₁,n₂=size(a)[1:d]
    return [sum(A[i,j]*a[I₁,I₂,j] for j ∈ 1:d)+b[i] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d]
end

function Positioning(a::Array{Float64,3}) # 制御点の位置調整
    n₁,n₂=size(a)[1:d]
    ind0=[(n₁+1)÷2,(n₂+1)÷2]
    ind1=ind0-[0,1]
    v=a[ind1...,:]-a[ind0...,:]
    R=-[v[2] -v[1];v[1] v[2]]/norm(v)
    return a=aff(a,R,-R*a[ind0...,:])
end

function Positioning(B2::Bs2mfd) # 制御点の位置調整
    p,k,a=B2.p,B2.k,B2.a
    n₁,n₂=size(a)[1:d]
    aa=Positioning(a)
    return Bs2mfd(p,k,aa)
end

function InitBs(D,n₁;nip=NIP)
    D₁,D₂=D

    g′₍₀₎(u)=reshape(ForwardDiff.jacobian(g₍₀₎,u),d,d,d) # 第一基本量の微分
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
    𝒄(t)=sol(t)[1:d] # 解となる中心曲線
    𝒄₁(t)=sol(t)[(d+1):(2d)] # その導関数
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
    a=[[a1[I₁,i],a2[I₁,i]][I₂] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d]

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

function elm_H(g₍₀₎,B2::Bs2mfd,I₁,I₂,i,R₁,R₂,r;nip=NIP)
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
                Q=[sum(a[I₁,I₂,i]*𝑁[I₁,I₂,j] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂) for i ∈ 1:d, j ∈ 1:d];
                sum(
                    C(p,q,m,n,g⁻)*𝑁[I₁,I₂,p]*(𝜹[i,r]*𝑁[R₁,R₂,q]*(sum(Q[o,m]*Q[o,n] for o ∈ 1:d)-g[m,n])+2*𝑁[R₁,R₂,n]*Q[i,q]*Q[r,m])
                for p ∈ 1:d, q ∈ 1:d, m ∈ 1:d, n ∈ 1:d)
            )*𝝊,(D̂₁,D̂₂),nip=nip
        )
    end
end

function elm_F(g₍₀₎,B2::Bs2mfd,I₁,I₂,i;nip=NIP)
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
            𝑁=[N′(B2,I₁,I₂,i,u) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d];
            Q=[sum(a[I₁,I₂,i]*𝑁[I₁,I₂,j] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂) for i ∈ 1:d, j ∈ 1:d];
            sum(
                sum(
                    C(p,q,m,n,g⁻)*𝑁[I₁,I₂,p]*Q[i,q]
                    for p ∈ 1:d, q ∈ 1:d
                )*(sum(
                    Q[o,m]*Q[o,n]
                for o ∈ 1:d)-g[m,n])
            for m ∈ 1:d, n ∈ 1:d)
        )*𝝊,(D̂₁,D̂₂),nip=nip
    )
end

function lineup(n,I₁,I₂,i)
    n₁,n₂=n
    return (i-1)*n₁*n₂+(I₂-1)*n₁+(I₁-1)+1
end

function NewtonIteration(B2::Bs2mfd,fixed;nip=NIP)
    n₁,n₂=n=length.(B2.k)-B2.p.-1

    t₀=time()
    Ff=Array{Any}(undef,n₁,n₂,d)
    for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d
        Ff[I₁,I₂,i]=@spawn elm_F(g₍₀₎,B2,I₁,I₂,i,nip=nip)
    end
    F=fetch.(Ff)
    Hf=Array{Any}(undef,n₁,n₂,d,n₁,n₂,d)
    for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d, R₁ ∈ 1:n₁, R₂ ∈ 1:n₂, r ∈ 1:d
        Hf[I₁,I₂,i,R₁,R₂,r]=@spawn elm_H(g₍₀₎,B2,I₁,I₂,i,R₁,R₂,r,nip=nip)
    end
    H=fetch.(Hf)
    t₁=time()

    # H=[elm_H(g₍₀₎,B2,I₁,I₂,i,R₁,R₂,r,nip=nip) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d, R₁ ∈ 1:n₁, R₂ ∈ 1:n₂, r ∈ 1:d]
    # F=[elm_F(g₍₀₎,B2,I₁,I₂,i,nip=nip) for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d]

    𝕟=2n₁*n₂
    Fixed=sort(collect((i->lineup(n,i...)).(fixed(n₁,n₂))))
    Unfixed=deleteat!(collect(1:𝕟),Fixed)

    F=reshape(F,𝕟)
    H=reshape(H,𝕟,𝕟)
    a=aₒ=reshape(B2.a,𝕟)
    Ȟ=H[Unfixed,Unfixed]
    ǎ=a[Unfixed]
    F̌=F[Unfixed]
    Ǧ=Ȟ\F̌
    ǎ=ǎ-Ǧ
    for i ∈ Fixed
        insert!(ǎ,i,aₒ[i])
    end
    a=reshape(ǎ,n₁,n₂,d)
    return (Bs2mfd(B2.p,B2.k,a),F,Ǧ,t₁-t₀)
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

function Export(B2::Bs2mfd,BsTree,BsJLD;comment="",maximumstrain=MAXIMUMSTRAIN)
    index=length(BsTree.nodes)
    BsJLD[string(index)]=B2
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
        BsDraw(B2,filename=Dir*"/nurbs/"*Name*"-"*string(index)*"_Bspline.svg",up=Up,down=Down,right=Right,left=Left,mesh=Mesh,unitlength=Unit)
        k₁,k₂=B2.k
        D=(k₁[1]..k₁[end],k₂[1]..k₂[end])

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
    if (isfile(DIR*"/"*NAME*".jld")) error("File already exists") end
    mkpath(DIR)
    mkpath(DIR*"/nurbs")
    mkpath(DIR*"/strain")
    mkpath(DIR*"/colorbar")
    mkpath(DIR*"/slack")
    BsJLD=Dict{String,Any}("Expr"=>EXPR)

    B2=InitBs(D,n₁,nip=nip)
    comment="Initial Configuration"
    BsTree=Tree()

    Export(B2,BsTree,BsJLD,comment=comment)
end

function p_Refinement(p₊::Array{Int64,1};parent=0)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    B2=BsJLD[string(parent)]

    B2=pref(B2,p₊)
    comment="p-refinement with "*string(p₊)
    addchild(BsTree,parent,comment)

    Export(B2,BsTree,BsJLD,comment=comment)
end

function h_Refinement(k₊::Array{Array{Float64,1},1};parent=0)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    B2=BsJLD[string(parent)]

    k₁₊,k₂₊=k₊
    k̄₁,k̄₂=[min(B2.k[i]...)..max(B2.k[i]...) for i in 1:2]

    if (*(Bool[kᵢ ∈ k̄₁ for kᵢ in k₁₊]...,Bool[kᵢ ∈ k̄₂ for kᵢ in k₂₊]...))
        B2=href(B2,k₊)
        comment="h-refinement with "*string(k₊)
        addchild(BsTree,parent,comment)
    else
        comment="error: can't compute h-refinement with "*string(k₊)
        addchild(BsTree,parent,comment)
    end
    Export(B2,BsTree,BsJLD,comment=comment)
end

function NewtonMethodIteration(;fixed=((n₁,n₂)->([(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1])),parent=0,nip=NIP)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    B2=BsJLD[string(parent)]

    n₁,n₂=length.(B2.k)-B2.p.-1
    if (!isodd(n₁*n₂)) error("n₁ and n₂ should be odd numbers") end
    B2=Positioning(B2)
    B2,F,Ǧ,Δt=NewtonIteration(B2,fixed,nip=nip)
    # comment="Newton Iteration - residual norm: "*(@sprintf("%.5e",norm(F)))*", Δa norm: "*(@sprintf("%.5e",norm(Ǧ)))*", computation time: "*(@sprintf("%.5e",Δt))*" sec"
    comment="Newton Iteration - residual norm: "*(@sprintf("%.5e",norm(F)))*", Δa norm: "*(@sprintf("%.5e",norm(Ǧ)))*", computation time: "*string(Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(1000Δt÷1))))
    addchild(BsTree,parent,comment)

    Export(B2,BsTree,BsJLD,comment=comment)
end

function Restoration()
    if (!isfile(DIR*"/"*NAME*".jld")) error("File doen't exists") end
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

function ComputeMaximumStrain(;index=0,mesh=tuple(20*[MESH...]...))
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (index==0) index=length(BsTree.nodes) end
    B2=BsJLD[string(index)]
    k₁,k₂=B2.k
    D=(k₁[1]..k₁[end],k₂[1]..k₂[end])

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

    κ₁=range(leftendpoint(D[1]),stop=rightendpoint(D[1]),length=mesh[1]+1)
    κ₂=range(leftendpoint(D[2]),stop=rightendpoint(D[2]),length=mesh[2]+1)

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
