using LinearAlgebra
using ForwardDiff
using Plots
using FastGaussQuadrature
using Statistics

# 基本的に等長的に等間隔
# ただし、曲率が急なところは密に
# 幅が狭いところも密に
# いい感じに無次元になるようにする

# 分割数は決まっていない
# initial_stateと同じ間隔になるように調整

# v(t) = ṗ(t)
# w(t) = v(t)*B(t)

# Function definitions
p(t) = 2*t + sinpi(t/2)
ṗ(t) = ForwardDiff.derivative(p, t)
B(t) = 0.5
ω(t) = abs(ṗ(t))/B(t)

# caluculate the length
t_min = 0
t_max = 3.9
n = 1001
sum([ω(t) for t in range(t_min,t_max,length=n)])*(t_max-t_min)/n

nodes, weights = gausslegendre(10)
nodes_shifted = t_min .+ (nodes .+ 1) ./ 2 .* (t_max-t_min)
dot(ω.(nodes_shifted), weights)*(t_max-t_min)/2

# first step
nodes, weights = gausslegendre(10)
t2 = Float64(t_min)
ts = [t2]
Ls = Float64[]

for _ in 1:100
    t1 = t2
    t2 = t1+1/ω(t1)
    for _ in 1:10
        nodes_shifted = t1 .+ (nodes .+ 1) ./ 2 .* (t2-t1)
        L12 = dot(ω.(nodes_shifted), weights)*(t2-t1)/2
        t2 += (1-L12)/ω(t2)
    end
    if t2 < t_max
        push!(ts, t2)
        push!(Ls, L12)
    else
        t2 = t_max
        nodes_shifted = t1 .+ (nodes .+ 1) ./ 2 .* (t2-t1)
        L12 = dot(ω.(nodes_shifted), weights)*(t2-t1)/2
        push!(ts, t2)
        push!(Ls, L12)
        break
    end
end

for _ in 1:10
    L̄ = mean(Ls)
    for i in 2:l-1
        ΔL = sum(Ls[1:i-1]) - L̄*(i-1)
        ts[i] -= ΔL / ω(ts[i])
    end
    for i in 1:l-1
        t1 = ts[i]
        t2 = ts[i+1]
        nodes_shifted = t1 .+ (nodes .+ 1) ./ 2 .* (t2-t1)
        L12 = dot(ω.(nodes_shifted), weights)*(t2-t1)/2
        Ls[i] = L12
    end
end

Ls

sum(Ls)

ts



