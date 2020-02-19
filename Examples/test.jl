using Distributed
addprocs(1);
@everywhere push!(LOAD_PATH, "Modules")
@everywhere push!(LOAD_PATH, ".")
using IntervalSets
using Printf
using Bspline
using I4SM

@DefineShape 𝒑₍₀₎(u)=[u...,u'*(u.*[-1.0,1.0])]
D=((-1.0)..1.0,0.0..0.1)
Settings("Paraboloid-1d",up=3,down=-3,right=3,left=-3,mesh=(20,1),unit=200,slack=true)

InitialConfiguration(D)

fixed(n₁,n₂)=[[1,(n₂+1)÷2,1],[1,(n₂+1)÷2,2],[n₁,(n₂+1)÷2,1],[n₁,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1]]
NewtonMethodIteration(nip=25,fixed=fixed)
NewtonMethodIteration(nip=45)
p_Refinement([0,1])
h_Refinement([Float64[],[0.05]])
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)
NewtonMethodIteration(nip=45)


# I4SM.N

using JLD
pnt=0

BsJLD=load(I4SM.DIR*"/"*I4SM.NAME*".jld")
BsTree=BsJLD["BsTree"]
B2=BsJLD["1"]



n₁,n₂=length.(B2.k)-B2.p.-1
if (!isodd(n₁*n₂)) error("n₁ and n₂ should be odd numbers") end
B2=I4SM.Positioning(B2)
fixed=((n₁,n₂)->([(n₁+1)÷2,(n₂+1)÷2,1],[(n₁+1)÷2,(n₂+1)÷2,2],[(n₁+1)÷2,(n₂+1)÷2-1,1]))
B2,F,Ǧ,Δt=I4SM.NewtonIteration(B2,fixed,nip=25)


Hf=Array{Any}(undef,n₁,n₂,d,n₁,n₂,d)
for I₁ ∈ 1:n₁, I₂ ∈ 1:n₂, i ∈ 1:d, R₁ ∈ 1:n₁, R₂ ∈ 1:n₂, r ∈ 1:d
    Hf[I₁,I₂,i,R₁,R₂,r]=@spawn I4SM.elm_H(I4SM.g₍₀₎,B2,I₁,I₂,i,R₁,R₂,r,nip=25)
end
H=fetch.(Hf)

using ElementaryCalculus
using LinearAlgebra
using BenchmarkTools
I4SM.elm_H(I4SM.g₍₀₎,B2,1,1,1,1,1,1,nip=25)
@benchmark I4SM.elm_H(I4SM.g₍₀₎,B2,1,1,1,1,1,1,nip=25)
@benchmark elm_H3(I4SM.g₍₀₎,B2,1,1,1,1,1,1,nip=25)
@benchmark elm_H3(I4SM.g₍₀₎,B2,1,1,1,1,1,1,nip=25)

function N′(B2::Bs2mfd,I₁,I₂,i,u)
    p,k,a=B2.p,B2.k,B2.a
    p₁,p₂=p
    k₁,k₂=k
    if (i==1)
        return Ḃs(I₁,p₁,k₁,u[1])*Bs(I₂,p₂,k₂,u[2])
    else
        return Bs(I₁,p₁,k₁,u[1])*Ḃs(I₂,p₂,k₂,u[2])
    end
end

function C(i,j,k,l,g⁻)
    return 𝝀*g⁻[i,j]*g⁻[k,l]+𝝁*(g⁻[i,k]*g⁻[j,l]+g⁻[i,l]*g⁻[j,k])
end

function elm_H3(g₍₀₎,B2::Bs2mfd,I₁,I₂,i,R₁,R₂,r;nip=NIP)
    p,k,a=B2.p,B2.k,B2.a
    p₁,p₂=p
    k₁,k₂=k
    n₁,n₂=length.(k)-p.-1
    D̂₁=Bsupp(I₁,p₁,k₁)∩Bsupp(R₁,p₁,k₁)
    D̂₂=Bsupp(I₂,p₂,k₂)∩Bsupp(R₂,p₂,k₂)
    𝜹=[1.0 0.0;0.0 1.0]
    if (isnullset(D̂₁)||isnullset(D̂₂))
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


1

# comment="Newton Iteration - residual norm: "*(@sprintf("%.5e",norm(F)))*", Δa norm: "*(@sprintf("%.5e",norm(Ǧ)))*", computation time: "*(@sprintf("%.5e",Δt))*" sec"
comment="Newton Iteration - residual norm: "*(@sprintf("%.5e",norm(F)))*", Δa norm: "*(@sprintf("%.5e",norm(Ǧ)))*", computation time: "*string(Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(1000Δt÷1))))
addchild(BsTree,parent,comment)




Export(B2,BsTree,BsJLD,comment=comment)
