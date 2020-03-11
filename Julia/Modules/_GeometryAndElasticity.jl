using ForwardDiff

# BSpline
function FittingBSpline(f, P::BSplineSpace; nip=NIP) # 1-dimensional
    p=P.degree
    k=P.knots
    D=k[1+p]..k[end-p]
    function a(i,j)
        D′=(max(k[i],k[j])..min(k[i+p+1],k[j+p+1])) ∩ D
        if width(D′)==0
            return 0
        else
            return GaussianQuadrature(t->BSplineBasis(i,P,t)*BSplineBasis(j,P,t), D′)
        end
    end
    n=dim(P)
    A=[a(i,j) for i ∈ 1:n, j ∈ 1:n]
    b=[GaussianQuadrature(t->BSplineBasis(i,P,t)*f(t), ((k[i]..k[i+p+1]) ∩ D)) for i ∈ 1:n]
    return inv(A)*b
end

function N′(P₁::BSplineSpace,P₂::BSplineSpace,I₁,I₂,i,u)
    if i==1
        return BSplineBasis′(I₁,P₁,u[1])*BSplineBasis(I₂,P₂,u[2])
    else
        return BSplineBasis(I₁,P₁,u[1])*BSplineBasis′(I₂,P₂,u[2])
    end
end

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
    if length(𝒫s) ≠ d
        error("dimension does not match")
    end

    p¹,p²=p=[M.bsplinespaces[i].degree for i ∈ 1:2]
    k¹,k²=k=[M.bsplinespaces[i].knots for i ∈ 1:2]

    n₁,n₂,_=size(𝒂)
    𝒂′=Positioning(𝒂)
    return BSplineManifold(𝒫s,𝒂′)
end

export Refinement
function BSpline.Refinement(;p₊::Union{Nothing,Array{Int,1}}=nothing, k₊::Union{Nothing,Array{Knots,1}}=nothing, parent=0)
    parent=Parent(parent)

    M=loadM(index=parent)

    comment="Refinement with p₊:"*string(p₊)*", k₊:"*string(k₊)

    M=BSpline.Refinement(M,p₊=p₊,k₊=k₊)
    Export(M,parent,comment=comment)
end

export ShowKnots
function ShowKnots(;index=0)
    M=loadM(index=index)

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


# Reference State
𝒑′₍₀₎(u) = ForwardDiff.jacobian(Main.𝒑₍₀₎,u) # Tangent vector
𝒑₁₍₀₎(u) = ForwardDiff.derivative(u₁->Main.𝒑₍₀₎([u₁,u[2]]),u[1])
𝒑₂₍₀₎(u) = ForwardDiff.derivative(u₂->Main.𝒑₍₀₎([u[1],u₂]),u[2])
𝒑₁₁₍₀₎(u) = ForwardDiff.derivative(u₁->𝒑₁₍₀₎([u₁,u[2]]),u[1])
𝒑₁₂₍₀₎(u) = ForwardDiff.derivative(u₂->𝒑₁₍₀₎([u[1],u₂]),u[2])
𝒑₂₁₍₀₎(u) = ForwardDiff.derivative(u₁->𝒑₂₍₀₎([u₁,u[2]]),u[1])
𝒑₂₂₍₀₎(u) = ForwardDiff.derivative(u₂->𝒑₂₍₀₎([u[1],u₂]),u[2])
𝒆₍₀₎(u) = normalize(cross(𝒑₁₍₀₎(u),𝒑₂₍₀₎(u))) # Normal vector
g₍₀₎(u) = 𝒑′₍₀₎(u)'𝒑′₍₀₎(u) # 第1基本量
g₍₀₎₁₁(u) = 𝒑₁₍₀₎(u)'𝒑₁₍₀₎(u)
g₍₀₎₁₂(u) = 𝒑₁₍₀₎(u)'𝒑₂₍₀₎(u)
g₍₀₎₂₁(u) = 𝒑₂₍₀₎(u)'𝒑₁₍₀₎(u)
g₍₀₎₂₂(u) = 𝒑₂₍₀₎(u)'𝒑₂₍₀₎(u)
h₍₀₎(u) = [(𝒆₍₀₎(u)'*𝒑₁₁₍₀₎(u)) (𝒆₍₀₎(u)'*𝒑₁₂₍₀₎(u)) ; (𝒆₍₀₎(u)'*𝒑₂₁₍₀₎(u)) (𝒆₍₀₎(u)'*𝒑₂₂₍₀₎(u))] # 第2基本量
K₍₀₎(u) = det(h₍₀₎(u))/det(g₍₀₎(u)) # Gaussian curvature
𝝊₍₀₎(u) = norm(cross(𝒑₁₍₀₎(u),𝒑₂₍₀₎(u))) # volume form
g⁻₍₀₎(u) = inv(g₍₀₎(u)) # 第1基本量の逆
g′₍₀₎(u) = reshape(ForwardDiff.jacobian(g₍₀₎,u),d,d,d) # 第1基本量の微分
𝛤₍₀₎²₁₁(u) = (g⁻₍₀₎(u)[2,1]*g′₍₀₎(u)[1,1,1]+g⁻₍₀₎(u)[2,2]*(2g′₍₀₎(u)[2,1,1]-g′₍₀₎(u)[1,1,2]))/2 # Christoffel symbol
e⁽⁰⁾₁(u)=normalize(𝒑₁₍₀₎(u))
e⁽⁰⁾₂(u)=normalize(𝒑₂₍₀₎(u) - (g₍₀₎₁₂(u)/g₍₀₎₁₁(u))*𝒑₁₍₀₎(u))

c(D₂,t)=[t,sum(extrema(D₂))/2] # 中心線に沿った座標
ṡ₍₀₎(D₂,t)=sqrt(g₍₀₎₁₁(c(D₂,t)))
s̈₍₀₎(D₂,t)=(1/2)*(g′₍₀₎(c(D₂,t)))[1,1,1]/sqrt(g₍₀₎₁₁(c(D₂,t)))
𝜅₍₀₎(D₂,t)=𝛤₍₀₎²₁₁(c(D₂,t))*𝝊₍₀₎(c(D₂,t))/ṡ₍₀₎(D₂,t)^3 # Geodesic curvature
K₍₀₎(D₂,t)=K₍₀₎(c(D₂,t)) # Gaussian curvature
B̃(D₂,t)=dot(e⁽⁰⁾₂(c(D₂,t)),𝒑₂₍₀₎(c(D₂,t)))*width(D₂)/2 # Breadth of the piece of surface


# Current State
𝒑₍ₜ₎(M,u)=Mapping(M,u)
function 𝒑′₍ₜ₎(M::BSplineManifold,u)
    P₁,P₂=M.bsplinespaces
    𝒂=M.controlpoints
    n₁,n₂,_=size(𝒂)
    return [sum(N′(P₁,P₂,I₁,I₂,j,u)*𝒂[I₁,I₂,i] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂) for i ∈ 1:d, j ∈ 1:d]
end
function 𝒑₁₍ₜ₎(M::BSplineManifold,u)
    P₁,P₂=M.bsplinespaces
    𝒂=M.controlpoints
    n₁,n₂,_=size(𝒂)
    return sum(N′(P₁,P₂,I₁,I₂,1,u)*𝒂[I₁,I₂,:] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂)
end
function 𝒑₂₍ₜ₎(M::BSplineManifold,u)
    P₁,P₂=M.bsplinespaces
    𝒂=M.controlpoints
    n₁,n₂,_=size(𝒂)
    return sum(N′(P₁,P₂,I₁,I₂,2,u)*𝒂[I₁,I₂,:] for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂)
end
g₍ₜ₎(M,u)=𝒑′₍ₜ₎(M,u)'𝒑′₍ₜ₎(M,u) # 第1基本量
g₍ₜ₎₁₁(M,u)=𝒑₁₍ₜ₎(M,u)'𝒑₁₍ₜ₎(M,u) # 第1基本量
g₍ₜ₎₁₂(M,u)=𝒑₁₍ₜ₎(M,u)'𝒑₂₍ₜ₎(M,u) # 第1基本量
g₍ₜ₎₂₁(M,u)=𝒑₂₍ₜ₎(M,u)'𝒑₁₍ₜ₎(M,u) # 第1基本量
g₍ₜ₎₂₂(M,u)=𝒑₂₍ₜ₎(M,u)'𝒑₂₍ₜ₎(M,u) # 第1基本量


# Strain
E(M,u)=(g₍ₜ₎(M,u)-g₍₀₎(u))/2
E₁₁(M::BSplineManifold,u)=(g₍ₜ₎₁₁(M,u)-g₍₀₎₁₁(u))/2
E⁽⁰⁾₁₁(M::BSplineManifold,u)=E₁₁(M,u)/g₍₀₎₁₁(u)

function Ẽ⁽⁰⁾₁₁(D₂::ClosedInterval,u)
    b=width(D₂)/2
    c=sum(extrema(D₂))/2
    r=(u[2]-c)/b
    return (1/2)*K₍₀₎(D₂,u[1])*B̃(D₂,u[1])^2*(r^2-1/3)
end

function Ẽ⁽⁰⁾₁₁(M::BSplineManifold,u)
    P₁,P₂=M.bsplinespaces
    p₁,p₂=P₁.degree,P₂.degree
    k₁,k₂=P₁.knots,P₂.knots
    D₂=k₂[1+p₂]..k₂[end-p₂]
    return Ẽ⁽⁰⁾₁₁(D₂,u)
end

function ComputeMaximumStrain(;index=0,mesh=tuple(20*[MESH...]...))
    M=loadM(index=index)
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

function PredictMaximumStrain(D;mesh=tuple(20*[MESH...]...))
    D₁,D₂=D

    κ₁=range(leftendpoint(D₁),stop=rightendpoint(D₁),length=mesh[1]+1)
    κ₂=range(leftendpoint(D₂),stop=rightendpoint(D₂),length=mesh[2]+1)

    E=[Ẽ⁽⁰⁾₁₁(D₂,[u₁,u₂]) for u₁ ∈ κ₁, u₂ ∈ κ₂]

    return (minimum(E),maximum(E))
end

export ShowMaximumStrain
function ShowMaximumStrain(D;index=0)
    minE,maxE=PredictMaximumStrain(D)

    println("Predicted: (min: ",minE,", max: ",maxE,")")

    if isTheShapeComputed()
        minE,maxE=ComputeMaximumStrain(index=index)
        println("Computed: (min: ",minE,", max: ",maxE,")")
    end

    return nothing
end


# Elastic Modulus
function C(i,j,k,l,g⁻)
    return 𝝀*g⁻[i,j]*g⁻[k,l]+𝝁*(g⁻[i,k]*g⁻[j,l]+g⁻[i,l]*g⁻[j,k])
end
