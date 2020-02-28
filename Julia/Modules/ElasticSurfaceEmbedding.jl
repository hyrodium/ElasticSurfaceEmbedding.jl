module ElasticSurfaceEmbedding

using LinearAlgebra
using Printf
using Distributed
using IntervalSets
using ForwardDiff
using Dates
using DifferentialEquations
using JLD
using JSON

using BSpline
using ParametricDraw

include("Constants.jl")
include("NumericalIntegral.jl")
include("Slack.jl")
include("TreeStructure.jl")
include("GeometryAndElasticity.jl")

export @DefineShape
macro DefineShape(ex)
    global EXPR=ex
    return :(@everywhere $ex)
end

export Settings
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

export Restoration
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

include("InitialConfiguration.jl")
include("NewtonRaphsonMethod.jl")

function affine(𝒂::Array{Float64,3},A::Array{Float64,2},b::Array{Float64,1})::Array{Float64,3}
    #x'=Ax+b
    n₁,n₂,_=size(𝒂)
    return [(A*𝒂[I₁,I₂,:]+b)[i] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d]
end

function Positioning(𝒂::Array{Float64,3})::Array{Float64,3} # 制御点の位置調整
    n₁,n₂,_=size(𝒂)
    ind0=[(n₁+1)÷2,(n₂+1)÷2]
    ind1=ind0-[0,1]
    v=𝒂[ind1...,:]-𝒂[ind0...,:]
    R=-[v[2] -v[1];v[1] v[2]]/norm(v)
    return 𝒂=affine(𝒂,R,-R*𝒂[ind0...,:])
end

function Positioning(M::BSplineManifold)::BSplineManifold # 制御点の位置調整
    𝒫s = M.bsplinespaces
    𝒂 = M.controlpoints
    if (length(𝒫s) ≠ d)
        error("dimension does not match")
    end

    p¹,p²=p=[M.bsplinespaces[i].degree for i ∈ 1:2]
    k¹,k²=k=[M.bsplinespaces[i].knots for i ∈ 1:2]

    n₁,n₂,_=size(𝒂)
    𝒂′=Positioning(𝒂)
    return BSplineManifold(𝒫s,𝒂′)
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

        rgb(u)=E⁽⁰⁾₁₁(M,u)*[1,-1,-1]/(2*MaximumStrain) .+0.5
        # draw strain distribution (6000x6000)
        ParametricColor(u->𝒑₍ₜ₎(M,u),D,rgb=rgb,filename=Dir*"/strain/"*Name*"-"*string(index)*"_strain.png",up=Up,down=Down,right=Right,left=Left,mesh=tuple(10*[Mesh...]...),unit=5*Unit[1])
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

export Refinement
function BSpline.Refinement(;p₊::Union{Nothing,Array{Int,1}}=nothing, k₊::Union{Nothing,Array{Knots,1}}=nothing, parent=0)
    BsJLD=load(DIR*"/"*NAME*".jld")
    BsTree=BsJLD["BsTree"]
    if (parent==0) parent=length(BsTree.nodes) end
    M=BsJLD[string(parent)]

    comment="refinement with "*string(p₊)*", "*string(k₊)
    addchild(BsTree,parent,comment)

    M=BSpline.Refinement(M,p₊=p₊,k₊=k₊)
    Export(M,BsTree,BsJLD,comment=comment)
end

export FinalOutput
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

    𝒆⁽⁰⁾₁(u)=normalize(𝒑₁₍ₜ₎(M,u))
    𝒆⁽⁰⁾₂(u)=[0.0 -1.0;1.0 0.0]*𝒆⁽⁰⁾₁(M,u)
    𝒑a(i,t)=𝒑₍ₜ₎(M,[t,leftendpoint(D₂)])+𝒆⁽⁰⁾₂([t,leftendpoint(D₂)])*i*cutout[1]/unitlength[1]
    𝒑b(i,t)=𝒑₍ₜ₎(M,[t,rightendpoint(D₂)])-𝒆⁽⁰⁾₂([t,rightendpoint(D₂)])*i*cutout[1]/unitlength[1]
    SvgCurve([[t->𝒑a(i,t) for i ∈ 0:cutout[2]]...,[t->𝒑b(i,t) for i ∈ 0:cutout[2]]...],D₁,filename=DIR*"/"*NAME*"-"*string(index)*"-cutout.svg",up=UP,down=DOWN,right=RIGHT,left=LEFT,thickness=0.1,mesh=mesh,unitlength=unitlength)

    if (SLACK)
        SlackFile(DIR*"/"*NAME*"-"*string(index)*"-final.svg")
        SlackFile(DIR*"/"*NAME*"-"*string(index)*"-cutout.svg")
    end

    return nothing
end

export ShowKnots
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

    κ₁=range(leftendpoint(D₁),stop=rightendpoint(D₁),length=mesh[1]+1)
    κ₂=range(leftendpoint(D₂),stop=rightendpoint(D₂),length=mesh[2]+1)

    E=[E⁽⁰⁾₁₁(M,[u₁,u₂]) for u₁ ∈ κ₁, u₂ ∈ κ₂]

    return (minimum(E),maximum(E))
end

export ShowMaximumStrain
function ShowMaximumStrain(;index=0,mesh=5)
    minE,maxE=ComputeMaximumStrain(index=index,mesh=mesh)
    println("min",minE,", max",maxE)

    return nothing
end

export ReDraw
function ReDraw()
    return nothing
end

end # module
