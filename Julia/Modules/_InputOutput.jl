using ParametricDraw

export @ParametricMapping
macro ParametricMapping(ex)
    global EXPR=ex
    if startswith(repr(EXPR),":(function 𝒑₍₀₎(u)\n") || startswith(repr(EXPR),":(𝒑₍₀₎(u) =")
        return :(@everywhere $EXPR)
    else
        error("Symbol of parametric mapping must be 𝒑₍₀₎(u)")
    end
end

export Settings
function Settings(name;up=5,down=-5,right=5,left=-5,mesh=(10,1),unit=100,slack=true,maximumstrain=0.0,overwrite=false)
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
    global OVERWRITE=overwrite
    return nothing
end

function isTheShapeComputed()
    return isfile(DIR*"/"*NAME*".jld")
end

function loadEMT(;index=0)
    if !isTheShapeComputed()
        error("jld file doesn't exists")
    end
    BsJLD=load(DIR*"/"*NAME*".jld")
    Expr=BsJLD["Expr"]
    BsTree=BsJLD["BsTree"]
    if index==0
        index=length(BsTree.nodes)
    end
    M=BsJLD[string(index)]
    return Expr, M, BsTree
end

export Restoration
function Restoration()
    Expr, _, BsTree =loadEMT()

    println(showtree(BsTree))
    global EXPR=BsJLD["Expr"]
    eval(:(@everywhere $EXPR))
    return nothing
end

function Export(M::BSplineManifold;comment="",maximumstrain=MAXIMUMSTRAIN,parent=0)
    if isTheShapeComputed()
        BsJLD=load(DIR*"/"*NAME*".jld")
        BsTree=BsJLD["BsTree"]
        addchild(BsTree,parent,comment)
    else
        BsJLD=Dict{String,Any}("Expr"=>EXPR)
        BsTree=Tree()
    end

    index=length(BsTree.nodes)
    BsJLD[string(index)]=M
    BsJLD["BsTree"]=BsTree
    save(DIR*"/"*NAME*".jld",BsJLD)
    println(showtree(BsTree))

    if maximumstrain==0.0
        MS=ComputeMaximumStrain(index=index)
        MaximumStrain=max(-MS[1],MS[2])
    else
        MaximumStrain=maximumstrain
    end

    if distributed
        @spawnat 1 ExportFiles(M,MaximumStrain,BsTree,index)
    else
        ExportFiles(M,MaximumStrain,BsTree,index)
    end

    return nothing
end

function ExportFiles(M::BSplineManifold,MaximumStrain,BsTree,index;Name=NAME,Dir=DIR,Up=UP,Down=DOWN,Right=RIGHT,Left=LEFT,Mesh=MESH,Unit=UNIT,Slack=SLACK)
    BSplineSvg(M,filename=Dir*"/nurbs/"*Name*"-"*string(index)*"_Bspline.svg",up=Up,down=Down,right=Right,left=Left,mesh=Mesh,unitlength=Unit)
    𝒂 = M.controlpoints
    P₁,P₂=P=M.bsplinespaces
    p₁,p₂=p=P₁.degree,P₂.degree
    k₁,k₂=k=P₁.knots,P₂.knots
    D₁,D₂=D=k₁[1+p₁]..k₁[end-p₁],k₂[1+p₂]..k₂[end-p₂]

    Width=(Right-Left)*Unit[1]
    Height=(Up-Down)*Unit[1]

    rgb(u)=E⁽⁰⁾₁₁(M,u)*[1,-1,-1]/(2*MaximumStrain) .+0.5
    # draw strain distribution (6000x6000)
    ParametricColor(u->𝒑₍ₜ₎(M,u),D,rgb=rgb,filename=Dir*"/strain/"*Name*"-"*string(index)*"_strain.png",up=Up,down=Down,right=Right,left=Left,mesh=tuple(10*[Mesh...]...),unit=5*Unit[1])
    ColorBar(max=MaximumStrain,filename=Dir*"/colorbar/"*Name*"-"*string(index)*"_colorbar.png",width=Width)

    # 1200x1200
    # svg to png (1600x1600)
    run(pipeline(`convert $(Dir*"/nurbs/"*Name*"-"*string(index)*"_Bspline.svg") $(Dir*"/nurbs/"*Name*"-"*string(index)*"_Bspline.png")`, stdout=devnull, stderr=devnull))
    # add colorbar to strain distribution figure (6000x6000)
    run(pipeline(`convert $(Dir*"/strain/"*Name*"-"*string(index)*"_strain.png") $(Dir*"/colorbar/"*Name*"-"*string(index)*"_colorbar.png") -gravity southeast -compose over -composite $(Dir*"/strain/"*Name*"-"*string(index)*"_swc.png")`, stdout=devnull, stderr=devnull))
    # resize png
    # (1200x1200)
    run(pipeline(`convert -resize $(Width)x$(Height) -unsharp 2x1.4+0.5+0 -quality 100 -verbose $(Dir*"/nurbs/"*Name*"-"*string(index)*"_Bspline.png") $(Dir*"/slack/"*Name*"-"*string(index)*"_Bspline.png")`, stdout=devnull, stderr=devnull))
    # (1200x1200)
    run(pipeline(`convert -resize $(Width)x$(Height) -unsharp 2x1.4+0.5+0 -quality 100 -verbose $(Dir*"/strain/"*Name*"-"*string(index)*"_swc.png") $(Dir*"/slack/"*Name*"-"*string(index)*"_strain.png")`, stdout=devnull, stderr=devnull))
    # ()
    run(pipeline(`convert $(Dir*"/slack/"*Name*"-"*string(index)*"_strain.png")  \( +clone -alpha opaque -fill white -colorize 100% \) +swap -geometry +0+0 -compose Over -composite -alpha off $(Dir*"/slack/"*Name*"-"*string(index)*"_strain.png")`, stdout=devnull, stderr=devnull))
    # append png imgs
    run(pipeline(`convert +append $(Dir*"/slack/"*Name*"-"*string(index)*"_Bspline.png") $(Dir*"/slack/"*Name*"-"*string(index)*"_strain.png") $(Dir*"/slack/"*Name*"-"*string(index)*"_append.png")`, stdout=devnull, stderr=devnull))

    if Slack
        SlackString(showtree(BsTree))
        SlackFile(Dir*"/slack/"*Name*"-"*string(index)*"_append.png")
    end
end

export FinalOutput
function FinalOutput(;index=0,unitlength=(10,"mm"),cutout=(0.1,5),mesh=60)
    _, M, BsTree =loadEMT(index=index)
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

    if SLACK
        SlackFile(DIR*"/"*NAME*"-"*string(index)*"-final.svg")
        SlackFile(DIR*"/"*NAME*"-"*string(index)*"-cutout.svg")
    end

    return nothing
end

export ReDraw
function ReDraw()
    return nothing
end
