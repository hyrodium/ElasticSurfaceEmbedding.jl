using IntervalSets
using FastGaussQuadrature


const NIP=45

function INT(f,D::ClosedInterval;nip=NIP)
    nodes, weights = gausslegendre(nip)
    return sum(
    weights.*
    [f(x) for
            x ∈ (width(D)*nodes.+sum(extrema(D)))/2
    ])*width(D)/2
end

function integrate0(f,D::ClosedInterval;n=1)
    a,b=endpoints(D)
    k=range(endpoints(D)...,length=n+1)
    Δ=width(D)/n
    S=0.0
    for i ∈ 1:n
        S+=f((k[i]+k[i+1])/2)
    end
    return Δ*S
end

function integrate1(f,D::ClosedInterval;n=1)
    a,b=endpoints(D)
    k=range(endpoints(D)...,length=n+1)
    S=0.0
    for i ∈ 1:n
        m=1
        Δ=k[i+1]-k[i]
        S+=Δ*(f(k[i])+f(k[i+1]))
    end
    return S/2
end

function integrate2(f,D::ClosedInterval;n=1)
    a,b=endpoints(D)
    k=range(endpoints(D)...,length=n+1)
    S=0.0
    for i ∈ 1:n
        m=2
        Δ=k[i+1]-k[i]
        nodes=range(k[i],k[i+1],length=m+1)
        weights=[1,4,1]
        S+=Δ*sum(weights.*f.(nodes))
    end
    return S/6
end


f(x)=x*x


integrate0(f,0..2,n=15)
integrate1(f,0..2,n=15)
integrate2(f,0..2,n=15)
INT(f,0..π)

integrate0(f,-π..π)
INT(f,-π..π)

range(3,5,length=8)
