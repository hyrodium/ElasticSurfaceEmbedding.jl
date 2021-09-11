using Test
using IntervalSets
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
    𝟙 = 0.99999
    xs = range(-B*𝟙, stop=B*𝟙, length=n+1)
    return maximum(f.(xs))-minimum(f.(xs))
end

dir_result_a = joinpath(@__DIR__, "result_a")
dir_result_b = joinpath(@__DIR__, "result_b")

rm(dir_result_b, recursive=true, force=true)
config_dir(dir_result_b)

@testset "Rhomboid" begin
    @parametric_mapping 𝒑₍₀₎(u) = [u...,u[1]+u[2]]
    D = (-1.0..1.0, -1.0..1.0)
    name = "Rhomboid"
    settings(name,canvas=(3,5),mesh=(20,1),unit=200,colorbarsize=0.3)

    initial_state(D, n₁=5)
    M = ElasticSurfaceEmbedding.loadM()
    𝒂 = controlpoints(M)
    @test 𝒂[1,1,:] ≈ [-√(3/2), -3/√(2)]
    @test 𝒂[1,2,:] ≈ [-√(3/2), -1/√(2)]
    @test 𝒂[1,3,:] ≈ [-√(3/2), 1/√(2)]
    @test 𝒂[3,1,:] ≈ [0, -2/√(2)]
    @test norm(𝒂[3,2,:]) < 1e-14
    @test 𝒂[3,3,:] ≈ [0, 2/√(2)]
    @test 𝒂[5,1,:] ≈ [√(3/2), -1/√(2)]
    @test 𝒂[5,2,:] ≈ [√(3/2), 1/√(2)]
    @test 𝒂[5,3,:] ≈ [√(3/2), 3/√(2)]

    newton_onestep()
    M = ElasticSurfaceEmbedding.loadM()
    𝒂 = controlpoints(M)
    @test 𝒂[1,1,:] ≈ [-√(3/2), -3/√(2)]
    @test 𝒂[1,2,:] ≈ [-√(3/2), -1/√(2)]
    @test 𝒂[1,3,:] ≈ [-√(3/2), 1/√(2)]
    @test 𝒂[3,1,:] ≈ [0, -2/√(2)]
    @test norm(𝒂[3,2,:]) < 1e-14
    @test 𝒂[3,3,:] ≈ [0, 2/√(2)]
    @test 𝒂[5,1,:] ≈ [√(3/2), -1/√(2)]
    @test 𝒂[5,2,:] ≈ [√(3/2), 1/√(2)]
    @test 𝒂[5,3,:] ≈ [√(3/2), 3/√(2)]
end

@testset "Planar" begin
    @parametric_mapping 𝒑₍₀₎(u) = [sin(u[1])*u[2], u[2]+cos(u[1])-u[1]^2/5, 0.0]
    # See https://www.desmos.com/calculator/4usvqpr0iu
    D = (-1.0..2.0, 1.0..1.2)
    name = "Planar"
    settings(name,canvas=(4,4),mesh=(20,1),unit=200,colorbarsize=0.3)

    initial_state(D, n₁=35)
    M = ElasticSurfaceEmbedding.loadM()
    @test norm([ElasticSurfaceEmbedding.E(M, [u¹, u²]) for u¹ in -0.9:0.1:1.9, u² in 1.05:0.05:1.15], Inf) < 1e-5

    newton_onestep()
    M = ElasticSurfaceEmbedding.loadM()
    @test norm([ElasticSurfaceEmbedding.E(M, [u¹, u²]) for u¹ in -0.9:0.1:1.9, u² in 1.05:0.05:1.15], Inf) < 1e-5
end

@testset "Sphere-thin" begin
    # For deriving analytical solution, see https://hackmd.io/@hyrodium/r1sCtEsLX
    L = 20
    B = 1/8

    @parametric_mapping 𝒑₍₀₎(u) = [cos(u[1])*cos(u[2]), sin(u[1])*cos(u[2]), sin(u[2])]
    D = (-L..L, -B..B)
    name = "Sphere-thin"
    settings(name,canvas=(2L,2),mesh=(L,1),unit=50,colorbarsize=0.05)

    initial_state(D, n₁=5)

    newton_onestep()
    newton_onestep()
    spline_refinement(p₊=[0,1],k₊=[Knots(-L+B,-L+2B,-L+3B,L-3B,L-2B,L-B),Knots(-B/2, 0, B/2)])
    newton_onestep()
    newton_onestep()

    M = ElasticSurfaceEmbedding.loadM()
    𝒂 = controlpoints(M)

    # Analytical
    k = sqrt(4atanh(tan(B/2))/(sin(B)/cos(B)^2+2atanh(tan(B/2))))
    # Numerical computed
    k̃ = 𝒑₁₍ₜ₎(M, [0,0])[1]
    # Approximated
    k̂ = 1-B^2/6

    # If the strip is thin, the analytical result k can be approximated with k̂.
    @test abs(log(k̃/k)) < 1e-4
    @test abs(log(k̂/k)) < 1e-4

    # Analytical
    h′(u²) = √(1-𝝂*(k^2/cos(u²)^2-1))
    # Numerical computed
    h̃′(u²) = 𝒑₂₍ₜ₎(M, [0,u²])[2]
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

    @parametric_mapping 𝒑₍₀₎(u) = [cos(u[1])*cos(u[2]), sin(u[1])*cos(u[2]), sin(u[2])]
    D = (-L..L, -B..B)
    name = "Sphere-thick"
    settings(name,canvas=(2L,2),mesh=(L,1),unit=50,colorbarsize=0.05)

    initial_state(D, n₁=5)

    newton_onestep()
    newton_onestep()
    spline_refinement(p₊=[0,1],k₊=[Knots(-L+B,-L+2B,-L+3B,L-3B,L-2B,L-B),Knots(-B/2, 0, B/2)])
    newton_onestep()
    newton_onestep()

    M = ElasticSurfaceEmbedding.loadM()
    𝒂 = controlpoints(M)

    # Analytical
    k = sqrt(4atanh(tan(B/2))/(sin(B)/cos(B)^2+2atanh(tan(B/2))))
    # Numerical computed
    k̃ = 𝒑₁₍ₜ₎(M, [0,0])[1]
    # Approximated
    k̂ = 1-B^2/6

    # If the strip is thick, the analytical result k cannot be approximated with k̂.
    @test abs(log(k̃/k)) < 1e-4
    @test abs(log(k̂/k)) > 1e-4

    # Analytical
    h′(u²) = √(1-𝝂*(k^2/cos(u²)^2-1))
    # Numerical computed
    h̃′(u²) = 𝒑₂₍ₜ₎(M, [0,u²])[2]
    # Approximated
    ĥ′(u²) = √(1+𝝂*(1-k̂^2))-(𝝂*k̂^2*u²^2)/(2*√(1+𝝂*(1-k̂^2)))

    # If the strip is thick, the analytical result h′ cannot be approximated with ĥ′.
    @test L²(h′,h̃′,B)/delta(h′,B) < 1e-2
    @test L²(h′,ĥ′,B)/delta(h′,B) > 1e-2

    ## Note
    # Try the following script to check the difference between analytical solution and numerical solution.
    # using Plots
    # 𝟙 = 0.99999
    # plot(h′,-B*𝟙,B*𝟙)
    # plot!(h̃′,-B*𝟙,B*𝟙)
    # plot!(ĥ′,-B*𝟙,B*𝟙)
end

@testset "Paraboloid" begin
    @parametric_mapping 𝒑₍₀₎(u) = [u...,u'*u]
    D(i,n) = (-1.0..1.0, (i-1)/n..i/n)
    name = "Paraboloid"
    settings(name,canvas=(4,4),mesh=(20,1),unit=200,colorbarsize=0.3)

    i=3
    initial_state(D(i,10))

    newton_onestep(fixingmethod=:fix3points)
    newton_onestep()
    spline_refinement(p₊=(0,1),k₊=(Knots(),Knots([(i-1/2)/10])))
    newton_onestep()
    newton_onestep()
    add_pin(tag="$name-$i")

    img_a = load(joinpath(dir_result_a,"Paraboloid","append","Paraboloid-5_append.png"))
    img_b = load(joinpath(dir_result_b,"Paraboloid","append","Paraboloid-5_append.png"))
    d = Euclidean()
    @test d(img_a, img_b) < 0.0001
end
