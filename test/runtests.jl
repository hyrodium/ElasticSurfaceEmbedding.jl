using Test
using IntervalSets
using StaticArrays
using Images
using LinearAlgebra
using BasicBSpline
using ElasticSurfaceEmbedding
using Aqua
import ElasticSurfaceEmbedding.𝝂
import ElasticSurfaceEmbedding.𝒑₁₍ₜ₎
import ElasticSurfaceEmbedding.𝒑₂₍ₜ₎

Aqua.test_all(ElasticSurfaceEmbedding; ambiguities=false)

function L²(f, B)
    n = 240
    𝟙 = 0.99999
    xs = range(-B * 𝟙, stop = B * 𝟙, length = n + 1)
    Δ = 2 * B * 𝟙 / n
    return sqrt(Δ * (2 * sum(f.(xs) .^ 2) - f(xs[begin])^2 - f(xs[end])^2) / 2)
end

function L²(f, g, B)
    return L²(x -> f(x) - g(x), B)
end

function delta(f, B)
    n = 10
    𝟙 = 1 - 1e-8
    xs = range(-B * 𝟙, stop = B * 𝟙, length = n + 1)
    return maximum(f.(xs)) - minimum(f.(xs))
end

DIR_RESULT = joinpath(@__DIR__, "result")

rm(DIR_RESULT, recursive = true, force = true)

@testset "README example" begin
    # Overload the shape definition
    ElasticSurfaceEmbedding.surface(x,y) = SVector(x, y, x^2+y^2)
    # (1) split the surface into strips
    dom = [(-1..1, (i-1)/10..i/10) for i in 1:10]
    # (2) Embed the strips onto a plane
    res = auto_allsteps(dom)
    export_pinned_steps(joinpath(DIR_RESULT, "paraboloid"), res)

    @test isfile(joinpath(DIR_RESULT, "paraboloid", "pinned", "pinned-7.svg"))
    @test isfile(joinpath(DIR_RESULT, "paraboloid", "pinned", "pinned-14.svg"))
    @test isfile(joinpath(DIR_RESULT, "paraboloid", "pinned", "pinned-21.svg"))
    @test isfile(joinpath(DIR_RESULT, "paraboloid", "pinned", "pinned-28.svg"))
    @test isfile(joinpath(DIR_RESULT, "paraboloid", "pinned", "pinned-35.svg"))
    @test isfile(joinpath(DIR_RESULT, "paraboloid", "pinned", "pinned-42.svg"))
    @test isfile(joinpath(DIR_RESULT, "paraboloid", "pinned", "pinned-49.svg"))
    @test isfile(joinpath(DIR_RESULT, "paraboloid", "pinned", "pinned-56.svg"))
    @test isfile(joinpath(DIR_RESULT, "paraboloid", "pinned", "pinned-63.svg"))
    @test isfile(joinpath(DIR_RESULT, "paraboloid", "pinned", "pinned-70.svg"))
end

@testset "Rhomboid" begin
    ElasticSurfaceEmbedding.𝒑₍₀₎(u¹, u²) = SVector(u¹, u², u¹ + u²)
    D = (-1.0 .. 1.0, -1.0 .. 1.0)
    show_strain(D)
    @test_logs (:info, "Strain - domain: [-1.0, 1.0]×[-1.0, 1.0]\nPredicted: (min: -0.0, max: 0.0)\n") show_strain(D)

    result = initial_state(D)
    M = ElasticSurfaceEmbedding.loadM(result)
    𝒂 = controlpoints(M)
    M, N = size(𝒂)
    m = M ÷ 2 + 1
    n = N ÷ 2 + 1
    @test 𝒂[1, 1] ≈ [-√(3 / 2), -3 / √(2)]
    @test 𝒂[1, n] ≈ [-√(3 / 2), -1 / √(2)]
    @test 𝒂[1, N] ≈ [-√(3 / 2), 1 / √(2)]
    @test 𝒂[m, 1] ≈ [0, -2 / √(2)]
    @test 𝒂[m, n] ≈ [0, 0] atol = 1e-14
    @test 𝒂[m, N] ≈ [0, 2 / √(2)]
    @test 𝒂[M, 1] ≈ [√(3 / 2), -1 / √(2)]
    @test 𝒂[M, n] ≈ [√(3 / 2), 1 / √(2)]
    @test 𝒂[M, N] ≈ [√(3 / 2), 3 / √(2)]

    newton_onestep!(result)
    M = ElasticSurfaceEmbedding.loadM(result)
    𝒂 = controlpoints(M)
    M, N = size(𝒂)
    m = M ÷ 2 + 1
    n = N ÷ 2 + 1
    @test 𝒂[1, 1] ≈ [-√(3 / 2), -3 / √(2)]
    @test 𝒂[1, n] ≈ [-√(3 / 2), -1 / √(2)]
    @test 𝒂[1, N] ≈ [-√(3 / 2), 1 / √(2)]
    @test 𝒂[m, 1] ≈ [0, -2 / √(2)]
    @test 𝒂[m, n] ≈ [0, 0] atol = 1e-14
    @test 𝒂[m, N] ≈ [0, 2 / √(2)]
    @test 𝒂[M, 1] ≈ [√(3 / 2), -1 / √(2)]
    @test 𝒂[M, n] ≈ [√(3 / 2), 1 / √(2)]
    @test 𝒂[M, N] ≈ [√(3 / 2), 3 / √(2)]
end

@testset "Planar" begin
    ElasticSurfaceEmbedding.𝒑₍₀₎(u¹, u²) = SVector(sin(u¹) * u², u² + cos(u¹) - u¹^2 / 5, 0.0)
    # See https://www.desmos.com/calculator/4usvqpr0iu
    D = (-1.0 .. 2.0, 1.0 .. 1.2)

    show_strain(D)
    result = initial_state(D)
    M = ElasticSurfaceEmbedding.loadM(result)
    @test norm([ElasticSurfaceEmbedding.E(M, u¹, u²) for u¹ in -0.9:0.1:1.9, u² in 1.05:0.05:1.15], Inf) < 1e-4

    newton_onestep!(result)
    refinement!(result, p₊=(0,1), k₊=suggest_knotvector(result))
    newton_onestep!(result)
    M = ElasticSurfaceEmbedding.loadM(result)
    @test norm([ElasticSurfaceEmbedding.E(M, u¹, u²) for u¹ in -0.9:0.1:1.9, u² in 1.05:0.05:1.15], Inf) < 1e-5

    @test result.pinned[2] == false
    pin!(result, 2)
    @test result.pinned[2] == true
    unpin!(result, 2)
    @test result.pinned[2] == false
end

@testset "Sphere-thin" begin
    # For deriving analytical solution, see https://hackmd.io/@hyrodium/r1sCtEsLX
    L = 20
    B = 1 / 8

    ElasticSurfaceEmbedding.𝒑₍₀₎(u¹, u²) = SVector(cos(u¹) * cos(u²), sin(u¹) * cos(u²), sin(u²))
    D = (-L .. L, -B .. B)

    show_strain(D)
    result = initial_state(D)
    newton_onestep!(result)
    newton_onestep!(result)
    refinement!(
        result,
        p₊ = (0, 1),
        k₊ = (KnotVector([-L + B, -L + 2B, -L + 3B, L - 3B, L - 2B, L - B]), KnotVector([-B / 2, 0.0, B / 2])),
    )
    newton_onestep!(result)
    newton_onestep!(result)

    M = ElasticSurfaceEmbedding.loadM(result)
    𝒂 = controlpoints(M)

    # Analytical
    k = sqrt(4atanh(tan(B / 2)) / (sin(B) / cos(B)^2 + 2atanh(tan(B / 2))))
    # Numerical computed
    k̃ = 𝒑₁₍ₜ₎(M, 0, 0)[1]
    # Approximated
    k̂ = 1 - B^2 / 6

    # If the strip is thin, the analytical result k can be approximated with k̂.
    @test abs(log(k̃ / k)) < 1e-4
    @test abs(log(k̂ / k)) < 1e-4

    # Analytical
    h′(u²) = √(1 - 𝝂 * (k^2 / cos(u²)^2 - 1))
    # Numerical computed
    h̃′(u²) = 𝒑₂₍ₜ₎(M, 0, u²)[2]
    # Approximated
    ĥ′(u²) = √(1 + 𝝂 * (1 - k̂^2)) - (𝝂 * k̂^2 * u²^2) / (2 * √(1 + 𝝂 * (1 - k̂^2)))

    # If the strip is thin, the analytical result h′ can be approximated with ĥ′.
    @test L²(h′, h̃′, B) / delta(h′, B) < 1e-2
    @test L²(h′, ĥ′, B) / delta(h′, B) < 1e-2
end


@testset "Sphere-thick" begin
    # For deriving analytical solution, see https://hackmd.io/@hyrodium/r1sCtEsLX
    L = 20
    B = 2 / 3

    ElasticSurfaceEmbedding.𝒑₍₀₎(u¹, u²) = SVector(cos(u¹) * cos(u²), sin(u¹) * cos(u²), sin(u²))
    D = (-L .. L, -B .. B)

    show_strain(D)
    result = initial_state(D)
    newton_onestep!(result)
    newton_onestep!(result)
    refinement!(
        result,
        p₊ = (0, 1),
        k₊ = (KnotVector([-L + B, -L + 2B, -L + 3B, L - 3B, L - 2B, L - B]), KnotVector([-B / 2, 0.0, B / 2])),
    )
    newton_onestep!(result)
    newton_onestep!(result)

    M = ElasticSurfaceEmbedding.loadM(result)
    𝒂 = controlpoints(M)

    # Analytical
    k = sqrt(4atanh(tan(B / 2)) / (sin(B) / cos(B)^2 + 2atanh(tan(B / 2))))
    # Numerical computed
    k̃ = 𝒑₁₍ₜ₎(M, 0, 0)[1]
    # Approximated
    k̂ = 1 - B^2 / 6

    # If the strip is thick, the analytical result k cannot be approximated with k̂.
    @test abs(log(k̃ / k)) < 1e-4
    @test abs(log(k̂ / k)) > 1e-4

    # Analytical
    h′(u²) = √(1 - 𝝂 * (k^2 / cos(u²)^2 - 1))
    # Numerical computed
    h̃′(u²) = 𝒑₂₍ₜ₎(M, 0, u²)[2]
    # Approximated
    ĥ′(u²) = √(1 + 𝝂 * (1 - k̂^2)) - (𝝂 * k̂^2 * u²^2) / (2 * √(1 + 𝝂 * (1 - k̂^2)))

    # If the strip is thick, the analytical result h′ cannot be approximated with ĥ′.
    @test L²(h′, h̃′, B) / delta(h′, B) < 1e-2
    @test L²(h′, ĥ′, B) / delta(h′, B) > 1e-2

    ## Note
    # Try the following script to check the difference between analytical solution and numerical solution.
    # using Plots
    # 𝟙 = 1 - 1e-8
    # plot(h′,-B*𝟙,B*𝟙)
    # plot!(h̃′,-B*𝟙,B*𝟙)
    # plot!(ĥ′,-B*𝟙,B*𝟙)
end

@testset "Paraboloid" begin
    ElasticSurfaceEmbedding.𝒑₍₀₎(u¹, u²) = SVector(u¹, u², u¹^2 + u²^2)
    name = "Paraboloid"

    N = 10
    result = StepTree()
    for i in 1:N
        D = (-1.0 .. 1.0, (i - 1) / N .. i / N)
        show_strain(D)
        result = initial_state!(result, D)
        newton_onestep!(result, fixingmethod = :fix3points)
        newton_onestep!(result)
        refinement!(result, p₊ = (0, 1), k₊ = (EmptyKnotVector(), KnotVector([(i - 1 / 2) / 10])))
        newton_onestep!(result)
        newton_onestep!(result)
        pin!(result)
    end

    export_all_steps(joinpath(DIR_RESULT, "Paraboloid"), result)
    files_pinned = readdir(joinpath(DIR_RESULT, "Paraboloid", "pinned"))

    @test length(files_pinned) == N

    # img_b = load(joinpath(DIR_RESULT,"Paraboloid","append","Paraboloid-5_append.png"))
    # d = Euclidean()
    # @test d(RGB.(img_a), RGB.(img_b)) < 0.0001
end
