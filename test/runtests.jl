using Test
using IntervalSets
using StaticArrays
using Images
using LinearAlgebra
using BasicBSpline
using ElasticSurfaceEmbedding
import ElasticSurfaceEmbedding.𝝂
import ElasticSurfaceEmbedding.𝒑₁₍ₜ₎
import ElasticSurfaceEmbedding.𝒑₂₍ₜ₎

function L²(f, B)
    n = 240
    𝟙 = 0.99999
    xs = range(-B*𝟙, stop=B*𝟙, length=n+1)
    Δ = 2*B*𝟙/n
    return sqrt(Δ*(2*sum(f.(xs).^2)-f(xs[begin])^2-f(xs[end])^2)/2)
end

function L²(f, g, B)
    return L²(x->f(x)-g(x),B)
end

function delta(f, B)
    n = 10
    𝟙 = 1 - 1e-8
    xs = range(-B*𝟙, stop=B*𝟙, length=n+1)
    return maximum(f.(xs))-minimum(f.(xs))
end

dir_result = joinpath(@__DIR__, "result")

rm(dir_result, recursive=true, force=true)

@testset "Rhomboid" begin
    ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(u¹,u²,u¹+u²)
    D = (-1.0..1.0, -1.0..1.0)
    show_strain(D)
    @test_logs (:info, "Strain - domain: [-1.0, 1.0]×[-1.0, 1.0]\nPredicted: (min: -0.0, max: 0.0)\n") show_strain(D)

    result = initial_state(D, n₁=5)
    M = ElasticSurfaceEmbedding.loadM(result)
    𝒂 = controlpoints(M)
    @test 𝒂[1,1] ≈ [-√(3/2), -3/√(2)]
    @test 𝒂[1,2] ≈ [-√(3/2), -1/√(2)]
    @test 𝒂[1,3] ≈ [-√(3/2), 1/√(2)]
    @test 𝒂[3,1] ≈ [0, -2/√(2)]
    @test 𝒂[3,2] ≈ [0, 0]  atol=1e-14
    @test 𝒂[3,3] ≈ [0, 2/√(2)]
    @test 𝒂[5,1] ≈ [√(3/2), -1/√(2)]
    @test 𝒂[5,2] ≈ [√(3/2), 1/√(2)]
    @test 𝒂[5,3] ≈ [√(3/2), 3/√(2)]

    newton_onestep!(result)
    M = ElasticSurfaceEmbedding.loadM(result)
    𝒂 = controlpoints(M)
    @test 𝒂[1,1] ≈ [-√(3/2), -3/√(2)]
    @test 𝒂[1,2] ≈ [-√(3/2), -1/√(2)]
    @test 𝒂[1,3] ≈ [-√(3/2), 1/√(2)]
    @test 𝒂[3,1] ≈ [0, -2/√(2)]
    @test norm(𝒂[3,2]) < 1e-14
    @test 𝒂[3,3] ≈ [0, 2/√(2)]
    @test 𝒂[5,1] ≈ [√(3/2), -1/√(2)]
    @test 𝒂[5,2] ≈ [√(3/2), 1/√(2)]
    @test 𝒂[5,3] ≈ [√(3/2), 3/√(2)]
end

@testset "Planar" begin
    ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(sin(u¹)*u², u²+cos(u¹)-u¹^2/5, 0.0)
    # See https://www.desmos.com/calculator/4usvqpr0iu
    D = (-1.0..2.0, 1.0..1.2)

    show_strain(D)
    result = initial_state(D, n₁=35)
    M = ElasticSurfaceEmbedding.loadM(result)
    @test norm([ElasticSurfaceEmbedding.E(M, u¹, u²) for u¹ in -0.9:0.1:1.9, u² in 1.05:0.05:1.15], Inf) < 1e-5

    newton_onestep!(result)
    M = ElasticSurfaceEmbedding.loadM(result)
    @test norm([ElasticSurfaceEmbedding.E(M, u¹, u²) for u¹ in -0.9:0.1:1.9, u² in 1.05:0.05:1.15], Inf) < 1e-5
end

@testset "Sphere-thin" begin
    # For deriving analytical solution, see https://hackmd.io/@hyrodium/r1sCtEsLX
    L = 20
    B = 1/8

    ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(cos(u¹)*cos(u²), sin(u¹)*cos(u²), sin(u²))
    D = (-L..L, -B..B)

    show_strain(D)
    result = initial_state(D, n₁=5)
    newton_onestep!(result)
    newton_onestep!(result)
    refinement!(result, p₊=(0,1),k₊=(KnotVector([-L+B,-L+2B,-L+3B,L-3B,L-2B,L-B]),KnotVector([-B/2, 0., B/2])))
    newton_onestep!(result)
    newton_onestep!(result)

    M = ElasticSurfaceEmbedding.loadM(result)
    𝒂 = controlpoints(M)

    # Analytical
    k = sqrt(4atanh(tan(B/2))/(sin(B)/cos(B)^2+2atanh(tan(B/2))))
    # Numerical computed
    k̃ = 𝒑₁₍ₜ₎(M,0,0)[1]
    # Approximated
    k̂ = 1-B^2/6

    # If the strip is thin, the analytical result k can be approximated with k̂.
    @test abs(log(k̃/k)) < 1e-4
    @test abs(log(k̂/k)) < 1e-4

    # Analytical
    h′(u²) = √(1-𝝂*(k^2/cos(u²)^2-1))
    # Numerical computed
    h̃′(u²) = 𝒑₂₍ₜ₎(M,0,u²)[2]
    # Approximated
    ĥ′(u²) = √(1+𝝂*(1-k̂^2))-(𝝂*k̂^2*u²^2)/(2*√(1+𝝂*(1-k̂^2)))

    # If the strip is thin, the analytical result h′ can be approximated with ĥ′.
    @test L²(h′,h̃′,B)/delta(h′,B) < 1e-2
    @test L²(h′,ĥ′,B)/delta(h′,B) < 1e-2
end


@testset "Sphere-thick" begin
    # For deriving analytical solution, see https://hackmd.io/@hyrodium/r1sCtEsLX
    L = 20
    B = 2/3

    ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(cos(u¹)*cos(u²), sin(u¹)*cos(u²), sin(u²))
    D = (-L..L, -B..B)

    show_strain(D)
    result = initial_state(D, n₁=5)
    newton_onestep!(result)
    newton_onestep!(result)
    refinement!(result, p₊=(0,1),k₊=(KnotVector([-L+B,-L+2B,-L+3B,L-3B,L-2B,L-B]),KnotVector([-B/2, 0., B/2])))
    newton_onestep!(result)
    newton_onestep!(result)

    M = ElasticSurfaceEmbedding.loadM(result)
    𝒂 = controlpoints(M)

    # Analytical
    k = sqrt(4atanh(tan(B/2))/(sin(B)/cos(B)^2+2atanh(tan(B/2))))
    # Numerical computed
    k̃ = 𝒑₁₍ₜ₎(M,0,0)[1]
    # Approximated
    k̂ = 1-B^2/6

    # If the strip is thick, the analytical result k cannot be approximated with k̂.
    @test abs(log(k̃/k)) < 1e-4
    @test abs(log(k̂/k)) > 1e-4

    # Analytical
    h′(u²) = √(1-𝝂*(k^2/cos(u²)^2-1))
    # Numerical computed
    h̃′(u²) = 𝒑₂₍ₜ₎(M,0,u²)[2]
    # Approximated
    ĥ′(u²) = √(1+𝝂*(1-k̂^2))-(𝝂*k̂^2*u²^2)/(2*√(1+𝝂*(1-k̂^2)))

    # If the strip is thick, the analytical result h′ cannot be approximated with ĥ′.
    @test L²(h′,h̃′,B)/delta(h′,B) < 1e-2
    @test L²(h′,ĥ′,B)/delta(h′,B) > 1e-2

    ## Note
    # Try the following script to check the difference between analytical solution and numerical solution.
    # using Plots
    # 𝟙 = 1 - 1e-8
    # plot(h′,-B*𝟙,B*𝟙)
    # plot!(h̃′,-B*𝟙,B*𝟙)
    # plot!(ĥ′,-B*𝟙,B*𝟙)
end

@testset "Paraboloid" begin
    ElasticSurfaceEmbedding.𝒑₍₀₎(u¹,u²) = SVector(u¹,u²,u¹^2+u²^2)
    name = "Paraboloid"

    N = 10
    result = AllSteps()
    for i in 1:N
        D = (-1.0..1.0, (i-1)/N..i/N)
        show_strain(D)
        result = initial_state!(result, D, n₁=25)
        newton_onestep!(result, fixingmethod=:fix3points)
        newton_onestep!(result)
        refinement!(result, p₊=(0,1), k₊=(EmptyKnotVector(),KnotVector([(i-1/2)/10])))
        newton_onestep!(result)
        newton_onestep!(result)
        pin!(result)
    end

    export_all_steps(joinpath(dir_result, "Paraboloid"), result)
    files_pinned =  readdir(joinpath(dir_result,"Paraboloid","pinned"))

    @test length(files_pinned) == N

    # img_b = load(joinpath(dir_result,"Paraboloid","append","Paraboloid-5_append.png"))
    # d = Euclidean()
    # @test d(RGB.(img_a), RGB.(img_b)) < 0.0001
end
