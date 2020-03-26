export PinState
function PinState(; parent::Int=0, tag::String="")
    if tag==""
        tag=Dates.format(now(),"yyyy-mm-dd_H-M-S")
    end
    CheckAsFineName(tag)
    if TagExists(tag)
        error("The tag $(tag) is already exists.")
    end

    parent = Parent(parent)
    M = loadM(index=parent)
    comment = "📌 - tag: "*tag
    Export(M,parent,comment=comment)
    return nothing
end

function TagExists(tag, dict::Union{Dict,Nothing}=nothing)::Bool
    if dict isa Nothing
        dict = LoadResultDict()
    end
    PinnedStates = FindPinnedStates()
    for i_key in PinnedStates
        index = parse(Int, i_key)
        if GetTag(index) == tag
            return true
        end
    end
    return false
end

function GetTag(index; dict::Union{Dict,Nothing}=nothing)::String
    if dict isa Nothing
        dict = LoadResultDict()
    end
    comment = dict["Result"][repr(index)]["comment"]
    if startswith(comment, "📌 ")
        return replace(comment, "📌 - tag: " => "")
    else
        error("The index $(index) is not pinned.")
    end
end

export RemovePin
function RemovePin(index)::Nothing
    dict = LoadResultDict()
    tag=GetTag(index)
    comment = dict["Result"][repr(index)]["comment"]
    comment = replace(comment, "📌" => "💨")
    dict["Result"][repr(index)]["comment"] = comment
    SaveResultDict(dict)
    return nothing
end

function FindPinnedStates(; dict::Union{Dict,Nothing}=nothing)::Array{String}
    if dict isa Nothing
        dict = LoadResultDict()
    end
    PinnedStates = String[]
    for i_key in keys(dict["Result"])
        comment = dict["Result"][i_key]["comment"]
        if startswith(comment, "📌 ")
            push!(PinnedStates,i_key)
        end
    end
    return PinnedStates
end

export ExportPinnedStates
function ExportPinnedStates(; unitlength=(10,"mm"),cutout=(0.1,5),mesh::Int=60)
    mkpath(DIR*"/pinned")
    dict = LoadResultDict()
    PinnedStates = FindPinnedStates(dict=dict)

    for i_key in PinnedStates
        index = parse(Int, i_key)

        M = loadM(index=index)
        BSplineSvg(M,filename=DIR*"/pinned/"*GetTag(index,dict=dict)*".svg",up=UP,down=DOWN,right=RIGHT,left=LEFT,mesh=MESH,unitlength=unitlength,points=false)

        P₁,P₂ = P = M.bsplinespaces
        p₁,p₂ = p = P₁.degree,P₂.degree
        k₁,k₂ = k = P₁.knots,P₂.knots
        D₁,D₂ = D = k₁[1+p₁]..k₁[end-p₁],k₂[1+p₂]..k₂[end-p₂]
        n₁,n₂ = n = dim.(P)

        𝒆⁽⁰⁾₁(u) = normalize(𝒑₁₍ₜ₎(M,u))
        𝒆⁽⁰⁾₂(u) = [0.0 -1.0;1.0 0.0]*𝒆⁽⁰⁾₁(u)
        𝒑a(i,t) = 𝒑₍ₜ₎(M,[t,leftendpoint(D₂)])+𝒆⁽⁰⁾₂([t,leftendpoint(D₂)])*i*cutout[1]/unitlength[1]
        𝒑b(i,t) = 𝒑₍ₜ₎(M,[t,rightendpoint(D₂)])-𝒆⁽⁰⁾₂([t,rightendpoint(D₂)])*i*cutout[1]/unitlength[1]
        SvgCurve([[t->𝒑a(i,t) for i ∈ 0:cutout[2]]...,[t->𝒑b(i,t) for i ∈ 0:cutout[2]]...],D₁,filename=DIR*"/pinned/"*GetTag(index,dict=dict)*"-cutout.svg",up=UP,down=DOWN,right=RIGHT,left=LEFT,thickness=0.1,mesh=mesh,unitlength=unitlength)
    end

    return nothing
end
