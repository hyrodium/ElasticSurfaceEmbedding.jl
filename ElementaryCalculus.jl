module ElementaryCalculus

using IntervalSets
using FastGaussQuadrature

export INT, INT2, isnullset, DelDpl

function INT(f,D::ClosedInterval;nip=25)
    nodes, weights = gausslegendre(nip)
    return sum(
    weights.*
    [f(x) for
            x in (width(D)*nodes.+sum(extrema(D)))/2
    ])*width(D)/2
end

function INT2(f,D;nip=25)
    nodes, weights = gausslegendre(nip)
    return sum(
    (weights*weights').*
    [f([x,y]) for
            x in (width(D[1])*nodes.+sum(extrema(D[1])))/2,
            y in (width(D[2])*nodes.+sum(extrema(D[2])))/2
    ])*prod(width.(D))/4
end

function isnullset(I::ClosedInterval)::Bool
    leftendpoint(I)-rightendpoint(I)≥0
end

function DelDpl(v::Array{T,1})::Array{Float64,1} where T<:Real
    w=Array{eltype(v),1}()
    if (length(v)>0)
        for e in v
            if e ∉ w
                push!(w,e)
            end
        end
    end
    return w
end

function CompareArray(A,B)
    dim=size(A)
    if(dim≠size(B))
        DimensionMismatch("dimensions must match")
    end
    dict=Dict{String,Any}()
    dict["Difference"]=A-B
    #dict["EuclidDistance"]=norm(reshape(A-B,prod(dim)))
    #dict["ChebyshevDistance"]=maximum(abs.(A-B))
    dict["NormalizedEuclidDistance"]=norm(reshape(A-B,prod(dim)))/norm(reshape(A,prod(dim)))
    dict["NormalizedChebyshevDistance"]=maximum(abs.(A-B))/norm(reshape(A,prod(dim)))
    ind=findall(elm->elm==maximum(abs.(A-B)),A-B)∪findall(elm->elm==maximum(abs.(B-A)),B-A)
    #println(ind)
    dict["ChebyshevDistance"]=(maximum(abs.(A-B)),(A[ind],B[ind]))
    return dict
end

end
