export PinState
function PinState(;parent::Int=0, tag::String="")
    if tag==""
        tag=Dates.format(now(),"yyyy-mm-dd_H-M-S")
    end
    CheckAsFineName(tag)
    parent=Parent(parent)
    M=loadM(index=parent)
    comment="📌 - tag: "*tag
    Export(M,parent,comment=comment)
end

function FindPinnedStates()
    dict=LoadResultDict()
    PinnedStates = String[]
    for i in keys(dict["Result"])
        comment = dict["Result"][i]["comment"]
        if startswith(comment, "📌 ")
            push!(PinnedStates,i)
        end
    end
    return PinnedStates
end

export ExportPinnedStates
function ExportPinnedStates(;unitlength=(10,"mm"),cutout=(0.1,5),mesh=60)
    mkpath(DIR*"/pinned")
    PinnedStates = FindPinnedStates()

    for i in PinnedStates
        index = parse(Int, i)

        M=loadM(index=index)
        BSplineSvg(M,filename=DIR*"/pinned/"*NAME*"-"*string(index)*"-final.svg",up=UP,down=DOWN,right=RIGHT,left=LEFT,mesh=MESH,unitlength=unitlength,points=false)

        P₁,P₂=P=M.bsplinespaces
        p₁,p₂=p=P₁.degree,P₂.degree
        k₁,k₂=k=P₁.knots,P₂.knots
        D₁,D₂=D=k₁[1+p₁]..k₁[end-p₁],k₂[1+p₂]..k₂[end-p₂]
        n₁,n₂=n=dim.(P)

        𝒆⁽⁰⁾₁(u)=normalize(𝒑₁₍ₜ₎(M,u))
        𝒆⁽⁰⁾₂(u)=[0.0 -1.0;1.0 0.0]*𝒆⁽⁰⁾₁(u)
        𝒑a(i,t)=𝒑₍ₜ₎(M,[t,leftendpoint(D₂)])+𝒆⁽⁰⁾₂([t,leftendpoint(D₂)])*i*cutout[1]/unitlength[1]
        𝒑b(i,t)=𝒑₍ₜ₎(M,[t,rightendpoint(D₂)])-𝒆⁽⁰⁾₂([t,rightendpoint(D₂)])*i*cutout[1]/unitlength[1]
        SvgCurve([[t->𝒑a(i,t) for i ∈ 0:cutout[2]]...,[t->𝒑b(i,t) for i ∈ 0:cutout[2]]...],D₁,filename=DIR*"/pinned/"*NAME*"-"*string(index)*"-cutout.svg",up=UP,down=DOWN,right=RIGHT,left=LEFT,thickness=0.1,mesh=mesh,unitlength=unitlength)
    end

    return nothing
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
