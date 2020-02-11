module ElementaryCalculus

using IntervalSets
using FastGaussQuadrature

export INT, INT₊, INT2, INT2₊, isnullset, DelDpl

const NIP=25 # Number of Integration Points

function INT(f,D::ClosedInterval;nip=NIP)
    nodes, weights = gausslegendre(nip)
    return sum(
    weights.*
    [f(x) for
            x ∈ (width(D)*nodes.+sum(extrema(D)))/2
    ])*width(D)/2
end

function INT₊(f,D;nip=NIP)
    if !isnullset(D)
        return INT(f,D,nip=nip)
    else
        return 0.0
    end
end

function INT2(f,D;nip=NIP)
    nodes, weights = gausslegendre(nip)
    return sum(
    (weights*weights').*
    [f([x,y]) for
            x ∈ (width(D[1])*nodes.+sum(extrema(D[1])))/2,
            y ∈ (width(D[2])*nodes.+sum(extrema(D[2])))/2
    ])*prod(width.(D))/4
end

function INT2₊(f,D;nip=NIP)
    if (*(((!).([isnullset.(D)...]))...))
        return INT2(f,D,nip=nip)
    else
        return 0.0
    end
end

function isnullset(I::ClosedInterval)::Bool
    leftendpoint(I)-rightendpoint(I)≥0
end

function DelDpl(v::Array{T,1})::Array{Float64,1} where T<:Real
    w=Array{eltype(v),1}()
    if (length(v)>0)
        for e ∈ v
            if e ∉ w
                push!(w,e)
            end
        end
    end
    return w
end

# function CompareArray(A,B)
#     dim=size(A)
#     if(dim≠size(B))
#         DimensionMismatch("dimensions must match")
#     end
#     dict=Dict{String,Any}()
#     dict["Difference"]=A-B
#     #dict["EuclidDistance"]=norm(reshape(A-B,prod(dim)))
#     #dict["ChebyshevDistance"]=maximum(abs.(A-B))
#     dict["NormalizedEuclidDistance"]=norm(reshape(A-B,prod(dim)))/norm(reshape(A,prod(dim)))
#     dict["NormalizedChebyshevDistance"]=maximum(abs.(A-B))/norm(reshape(A,prod(dim)))
#     ind=findall(elm->elm==maximum(abs.(A-B)),A-B)∪findall(elm->elm==maximum(abs.(B-A)),B-A)
#     #println(ind)
#     dict["ChebyshevDistance"]=(maximum(abs.(A-B)),(A[ind],B[ind]))
#     return dict
# end

end
