using Plots
using IntervalSets
using BasicBSpline
using BasicBSplineExporter
using StaticArrays
using ElasticSurfaceEmbedding
using Rotations
using Statistics

const N = 6
const J = 1
f0(s) = max(-abs(s+1/2N-1)-(1/2N-1), 0)
f1(s) = -1/2+f0(mod(s-J/N, 2))
f2(s) = 1/2-f0(mod(s-1-J/N, 2))
# plot(f0, -2, 3)

plot(aspectratio=1)
plot!(f1, 0, 2)
plot!(f2, 0, 2)

plot(aspectratio=1)
plot!(f1, -1/N, 2-1/N)
plot!(f2, -1/N, 2-1/N)

plot(aspectratio=1)
plot!(f1, -2/N, 2-2/N)
plot!(f2, -2/N, 2-2/N)



function merge(M1, M2)
    p₁ = 3
    k₁ = M1.bsplinespaces[1].knotvector[1:end-p₁-1] + M2.bsplinespaces[1].knotvector[2:end]
    P₁ = BSplineSpace{p₁}(k₁)
    P₂ = M1.bsplinespaces[2]

    v1 = M1.controlpoints[end,:]
    v2 = M2.controlpoints[1,:]
    r = rotation_between(v2[end] - v2[1], v1[end] - v1[1])
    w2 = [r*v for v in v2]
    c1 = sum(v1)/length(v1)
    c2 = sum(w2)/length(w2)
    a1 = M1.controlpoints
    a2 = [r*a-c2+c1 for a in M2.controlpoints]
    a = vcat(a1[1:end-1, :], (a1[end:end, :]+a2[1:1, :])/2, a2[2:end, :])
    M = BSplineManifold(a, P₁, P₂)
    return M
end

function merge(Ms::Vector{<:BSplineManifold})
    M = Ms[1]
    for i in 2:length(Ms)
        M = merge(M, Ms[i])
    end
    return M
end

u(s,t) = π*s
v(s,t) = π*(f1(s)*(1-t) + t*f2(s))

# 0≤u≤2π, -π/2≤v≤π/2
# 0≤s≤2, 0≤t≤1

catenoid(u,v) = SVector(cos(u)*cosh(v),sin(u)*cosh(v),v)
ElasticSurfaceEmbedding.𝒑₍₀₎(s,t) = catenoid(u(s,t), v(s,t))

steptree = StepTree()

# strip-01
plot(aspectratio=1)
plot!(f1, 0, 2)
plot!(f2, 0, 2)
intervals = [0..1/2N, 1/2N..1/N, 1/N..1, 1..(1+1/2N), (1+1/2N)..(1+1/N), (1+1/N)..2]
initial_state!(steptree, (intervals[1], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[2], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[3], 0..1), n₁=43)
newton_onestep!(steptree, fixingmethod=:fix5points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[4], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[5], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[6], 0..1), n₁=43)
newton_onestep!(steptree, fixingmethod=:fix5points)
for _ in 1:5 newton_onestep!(steptree) end
manifolds = [step.manifold for step in steptree.steps[(-length(intervals)+1:0) .* 7 .+ end]]
M = merge(manifolds)
info = Dict(["type" => "custom"])
step = ElasticSurfaceEmbedding.Step(M, "Custom", info)
ElasticSurfaceEmbedding.addstep!(steptree, step, 0)
for _ in 1:2 newton_onestep!(steptree) end
refinement!(steptree, p₊=(0,1), k₊=(EmptyKnotVector(), KnotVector([0.5])))
for _ in 1:6 newton_onestep!(steptree) end
pin!(steptree)

# strip-02
plot(aspectratio=1)
plot!(f1, -1/N, 2-1/N)
plot!(f2, -1/N, 2-1/N)
intervals = [-1/N..0, 0..1/2N, 1/2N..1/N, 1/N..1, 1..(1+1/2N), (1+1/2N)..(1+1/N), (1+1/N)..(2-1/N)]
initial_state!(steptree, (intervals[1], 0..1), n₁=7)
newton_onestep!(steptree, fixingmethod=:fix5points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[2], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[3], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[4], 0..1), n₁=43)
newton_onestep!(steptree, fixingmethod=:fix5points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[5], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[6], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[7], 0..1), n₁=39)
newton_onestep!(steptree, fixingmethod=:fix5points)
for _ in 1:5 newton_onestep!(steptree) end
manifolds = [step.manifold for step in steptree.steps[(-length(intervals)+1:0) .* 7 .+ end]]
M = merge(manifolds)
info = Dict(["type" => "custom"])
step = ElasticSurfaceEmbedding.Step(M, "Custom", info)
ElasticSurfaceEmbedding.addstep!(steptree, step, 0)
for _ in 1:2 newton_onestep!(steptree) end
refinement!(steptree, p₊=(0,1), k₊=(EmptyKnotVector(), KnotVector([0.5])))
for _ in 1:6 newton_onestep!(steptree) end
pin!(steptree)


# strip-03
plot(aspectratio=1)
plot!(f1, -2/N, 2-2/N)
plot!(f2, -2/N, 2-2/N)
intervals = [-2/N..0, 0..1/2N, 1/2N..1/N, 1/N..1, 1..(1+1/2N), (1+1/2N)..(1+1/N), (1+1/N)..(2-2/N)]
initial_state!(steptree, (intervals[1], 0..1), n₁=11)
newton_onestep!(steptree, fixingmethod=:fix5points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[2], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[3], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[4], 0..1), n₁=43)
newton_onestep!(steptree, fixingmethod=:fix5points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[5], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[6], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[7], 0..1), n₁=35)
newton_onestep!(steptree, fixingmethod=:fix5points)
for _ in 1:5 newton_onestep!(steptree) end
manifolds = [step.manifold for step in steptree.steps[(-length(intervals)+1:0) .* 7 .+ end]]
M = merge(manifolds)
info = Dict(["type" => "custom"])
step = ElasticSurfaceEmbedding.Step(M, "Custom", info)
ElasticSurfaceEmbedding.addstep!(steptree, step, 0)
for _ in 1:2 newton_onestep!(steptree) end
refinement!(steptree, p₊=(0,1), k₊=(EmptyKnotVector(), KnotVector([0.5])))
for _ in 1:6 newton_onestep!(steptree) end
pin!(steptree)

export_pinned_steps("fff", steptree, xlims=(-5,5), ylims=(-5,5), unitlength=(100/π, "mm"))


export_all_steps("hogehoge5", steptree)


intervals = [0..1/2N, 1/2N..1/N, 1/N..1, 1..(1+1/2N), (1+1/2N)..(1+1/N), (1+1/N)..2]
initial_state!(steptree, (intervals[1], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[2], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[3], 0..1), n₁=43)
newton_onestep!(steptree, fixingmethod=:fix5points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[4], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[5], 0..1), n₁=5)
newton_onestep!(steptree, fixingmethod=:fix3points)
for _ in 1:5 newton_onestep!(steptree) end
initial_state!(steptree, (intervals[6], 0..1), n₁=43)
newton_onestep!(steptree, fixingmethod=:fix5points)
for _ in 1:5 newton_onestep!(steptree) end

manifolds = [step.manifold for step in steptree.steps[(-length(intervals)+1:0) .* 7 .+ end]]
M = merge(manifolds)
info = Dict(["type" => "custom"])
step = ElasticSurfaceEmbedding.Step(M, "Custom", info)
ElasticSurfaceEmbedding.addstep!(steptree, step, 0)

for _ in 1:5 newton_onestep!(steptree) end
pin!(steptree)


# export_one_step




save_png("a3.png", M, ylims=(-5,10))



c1
sum(w2)/length(w2)

v1[end] - v1[1]









export_all_steps("hogehoge5", steptree, unitlength=(100, "mm"), xlims=(-5,5), ylims=(-5,5))








steptree.steps[end].manifold.bsplinespaces[1].knotvector.vector

