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
    if (i==1)
        return BSplineBasis′(I₁,P₁,u[1])*BSplineBasis(I₂,P₂,u[2])
    else
        return BSplineBasis(I₁,P₁,u[1])*BSplineBasis′(I₂,P₂,u[2])
    end
end

# Reference State
𝒑′₍₀₎(u)=ForwardDiff.jacobian(Main.𝒑₍₀₎,u) # Tangent vector
𝒑₁₍₀₎(u)=ForwardDiff.derivative(u₁->Main.𝒑₍₀₎([u₁,u[2]]),u[1])
𝒑₂₍₀₎(u)=ForwardDiff.derivative(u₂->Main.𝒑₍₀₎([u[1],u₂]),u[2])
𝒑₁₁₍₀₎(u)=ForwardDiff.derivative(u₁->Main.𝒑₁₍₀₎([u₁,u[2]]),u[1])
𝒑₁₂₍₀₎(u)=ForwardDiff.derivative(u₂->Main.𝒑₁₍₀₎([u[1],u₂]),u[2])
𝒑₂₁₍₀₎(u)=ForwardDiff.derivative(u₁->Main.𝒑₂₍₀₎([u₁,u[2]]),u[1])
𝒑₂₂₍₀₎(u)=ForwardDiff.derivative(u₂->Main.𝒑₂₍₀₎([u[1],u₂]),u[2])
𝒆₍₀₎(u)=normalize(cross(𝒑₁₍₀₎(u),𝒑₂₍₀₎(u))) # Normal vector
g₍₀₎(u)=𝒑′₍₀₎(u)'𝒑′₍₀₎(u) # 第1基本量
g₍₀₎₁₁(u)=𝒑₁₍₀₎(u)'𝒑₁₍₀₎(u)
g₍₀₎₁₂(u)=𝒑₁₍₀₎(u)'𝒑₂₍₀₎(u)
g₍₀₎₂₁(u)=𝒑₂₍₀₎(u)'𝒑₁₍₀₎(u)
g₍₀₎₂₂(u)=𝒑₂₍₀₎(u)'𝒑₂₍₀₎(u)
h₍₀₎(u)=[(𝒆₍₀₎(u)'*𝒑₁₁₍₀₎(u)) (𝒆₍₀₎(u)'*𝒑₁₂₍₀₎(u)) ; (𝒆₍₀₎(u)'*𝒑₂₁₍₀₎(u)) (𝒆₍₀₎(u)'*𝒑₂₂₍₀₎(u))] # 第2基本量
K₍₀₎(u::Array{Float64,1})=det(h₍₀₎(u))/det(g₍₀₎(u)) # Gaussian curvature
𝝊₍₀₎(u)=norm(cross(𝒑₁₍₀₎(u),𝒑₂₍₀₎(u))) # volume form
g⁻₍₀₎(u)=inv(g₍₀₎(u)) # 第1基本量の逆
g′₍₀₎(u)=reshape(ForwardDiff.jacobian(g₍₀₎,u),d,d,d) # 第1基本量の微分
𝛤₍₀₎²₁₁(u)=(g⁻₍₀₎(u)[2,1]*g′₍₀₎(u)[1,1,1]+g⁻₍₀₎(u)[2,2]*(2g′₍₀₎(u)[2,1,1]-g′₍₀₎(u)[1,1,2]))/2 # Christoffel symbol

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
E₁₁(M,u)=(g₍ₜ₎₁₁(M,u)-g₍₀₎₁₁(u))/2
E⁽⁰⁾₁₁(M,u)=E₁₁(M,u)/g₍₀₎₁₁(u)

# Elastic Modulus
function C(i,j,k,l,g⁻)
    return 𝝀*g⁻[i,j]*g⁻[k,l]+𝝁*(g⁻[i,k]*g⁻[j,l]+g⁻[i,l]*g⁻[j,k])
end
