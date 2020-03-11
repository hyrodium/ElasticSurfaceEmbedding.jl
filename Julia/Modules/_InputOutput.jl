using ParametricDraw

export @ParametricMapping
macro ParametricMapping(ex)
    expr=toJSON(ex)
    if startswith(expr,"function 𝒑₍₀₎(u)\n") || startswith(expr,"𝒑₍₀₎(u) =")
        global EXPR=expr
    else
        error("Symbol of parametric mapping must be 𝒑₍₀₎(u)")
    end
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
    if isTheShapeComputed()
        dict=LoadResultDict()
        println(TreeString(dict["Result"]))
        global EXPR=dict["Expr"]
    elseif !(@isdefined EXPR)
        error("use @ParametricMapping or input name of computed shape")
    end
    eval(:(@everywhere $(Meta.parse(EXPR))))
    return nothing
end

function toJSON(ex::Expr) ::String
    return replace(string(ex), r"#=.*=#\n" => s"")
end
function toJSON(k::Knots)
    return k.vector
end
function toJSON(P::BSplineSpace)
    return Dict("degree"=>P.degree, "knots"=>toJSON(P.knots))
end
function toJSON(P::Array{BSplineSpace,1})
    return convert(Array{Any,1}, toJSON.(P))
end
function toJSON(𝒂::Array{Float64})
    return convert(Array{Any,1}, 𝒂[:])
end
function toJSON(M::BSplineManifold)
    return Dict("bsplinespaces"=>toJSON(M.bsplinespaces), "controlpoints"=>toJSON(M.controlpoints))
end

function JSONtoBSplineSpace(jP::Dict)
    return BSplineSpace(jP["degree"],Knots(jP["knots"]))
end
function JSONtoBSplineSpaces(jPs::Array)
    return [JSONtoBSplineSpace(jP) for jP ∈ jPs]
end
function JSONtoControlPoints(j𝒂, dims)
    n = prod(dims)
    d = length(j𝒂) ÷ n
    return reshape(convert(Array{Float64},j𝒂),dims...,d)
end
function JSONtoBSplineManifold(dict::Dict)
    P = JSONtoBSplineSpaces(dict["bsplinespaces"])
    dims = dim.(P)
    𝒂 = JSONtoControlPoints(dict["controlpoints"], dims)
    return BSplineManifold(P,𝒂)
end

function NodeSeries(tree::Dict,node)
    Nodes=[node]
    while Nodes[end] ≠ "0"
        push!(Nodes, tree[Nodes[end]]["parent"])
    end
    return Nodes
end

function TreeString(tree::Dict)
    serieses=Array{Int,1}[]
    for key in keys(tree)
        push!(serieses,(s->parse(Int,s)).(reverse(NodeSeries(tree,key))))
    end
    sort!(serieses)
    lowstrings=String[]
    n = length(serieses)
    for i in 1:n
        l=length(serieses[i])
        key=string(serieses[i][end])
        comment=tree[key]["comment"]
        if l == 2
            lowstring=key*": "*comment
            push!(lowstrings,lowstring)
        elseif l ≥ 3
            lowstring="  "^(l-3)*"└─"*key*": "*comment
            push!(lowstrings,lowstring)
            for j in 1:(i-1)
                chars=collect(lowstrings[end-j])
                if chars[2(l-3)+1]==' '
                    lowstrings[end-j]=join(chars[1:2(l-3)])*"│"*join(chars[2(l-3)+2:end])
                elseif chars[2(l-3)+1]=='└'
                    lowstrings[end-j]=join(chars[1:2(l-3)])*"├"*join(chars[2(l-3)+2:end])
                    break
                else
                    break
                end
            end
        end
    end
    outsting=""
    for s in lowstrings
        outsting=outsting*s*"\n"
    end
    return outsting
end


function isTheShapeComputed()
    return isfile(DIR*"/"*NAME*".json")
end

function Parent(index::Int)
    dict=LoadResultDict()
    if index==0
        return length(keys(dict["Result"]))
    else
        return index
    end
end

function loadM(;index=0)
    if !isTheShapeComputed()
        error("Result file doesn't exists")
    end
    dict=LoadResultDict()
    if EXPR ≠ dict["Expr"]
        println(EXPR)
        println(dict["Expr"])
        # error("The definition of 𝒑₍₀₎(u) has been changed")
    end
    index=Parent(index)
    M=JSONtoBSplineManifold(dict["Result"][string(index)]["bsplinemanifold"])
    return M
end

function LoadResultDict()
    f=open(DIR*"/"*NAME*".json","r")
    dict=JSON.parse(f)
    close(f)
    return dict
end

function Export(M::BSplineManifold,parent::Int;comment="",maximumstrain=MAXIMUMSTRAIN)
    if isTheShapeComputed()
        dict=LoadResultDict()
        index=Parent(0)+1
    else
        dict=Dict{String,Any}("Expr" => EXPR)
        dict["Result"]=Dict{String,Any}()
        index=1
    end

    dict["Result"][string(index)]=Dict{String,Any}("parent" => string(parent))
    dict["Result"][string(index)]["bsplinemanifold"]=toJSON(M)
    dict["Result"][string(index)]["comment"]=comment
    open(DIR*"/"*NAME*".json","w") do f
        JSON.print(f, dict,4)
    end
    println(TreeString(dict["Result"]))

    if maximumstrain==0.0
        MS=ComputeMaximumStrain(index=index)
        MaximumStrain=max(-MS[1],MS[2])
    else
        MaximumStrain=maximumstrain
    end

    if distributed
        @spawnat 1 ExportFiles(M,MaximumStrain,index)
    else
        ExportFiles(M,MaximumStrain,index)
    end

    return nothing
end

function ExportFiles(M::BSplineManifold,MaximumStrain,index;Name=NAME,Dir=DIR,Up=UP,Down=DOWN,Right=RIGHT,Left=LEFT,Mesh=MESH,Unit=UNIT,Slack=SLACK)
    dict=LoadResultDict()
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
        SlackString("```\n"*TreeString(dict["Result"])*"```")
        SlackFile(Dir*"/slack/"*Name*"-"*string(index)*"_append.png")
    end
end

export FinalOutput
function FinalOutput(;index=0,unitlength=(10,"mm"),cutout=(0.1,5),mesh=60)
    M=loadM(index=index)
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
