using Test
using IntervalSets
using BasicBSpline
using ElasticSurfaceEmbedding
using Images
using LinearAlgebra

dir_result_a = joinpath(@__DIR__, "result_a")
dir_result_b = joinpath(@__DIR__, "result_b")

rm(dir_result_b, recursive=true, force=true)
config_dir(dir_result_b)

@testset "Rhomboid" begin
    @parametric_mapping 𝒑₍₀₎(u) = [u...,u[1]+u[2]]
    D = (-1.0..1.0, -1.0..1.0)
    name = "Rhomboid"
    settings(name,canvas=(4,4),mesh=(20,1),unit=200,colorbarsize=0.3)

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
    D = (-1.0..2.0, 1.0..1.2)
    name = "Planar"
    settings(name,canvas=(4,4),mesh=(20,1),unit=200,colorbarsize=0.3)

    initial_state(D, n₁=35)
    M = ElasticSurfaceEmbedding.loadM()
    @test norm([tr(ElasticSurfaceEmbedding.E(M, [u¹, u²])) for u¹ in -0.9:0.1:1.9, u² in 1.05:0.05:1.15], Inf) < 1e-5

    newton_onestep()
    M = ElasticSurfaceEmbedding.loadM()
    @test norm([tr(ElasticSurfaceEmbedding.E(M, [u¹, u²])) for u¹ in -0.9:0.1:1.9, u² in 1.05:0.05:1.15], Inf) < 1e-5
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
    spline_refinement(p₊=[0,1],k₊=[Knots(),Knots([(i-1/2)/10])])
    newton_onestep()
    newton_onestep()
    add_pin(tag="$name-$i")

    img_a = load(joinpath(dir_result_a,"Paraboloid","append","Paraboloid-5_append.png"))
    img_b = load(joinpath(dir_result_b,"Paraboloid","append","Paraboloid-5_append.png"))
    d = Euclidean()
    @test d(img_a, img_b) < 0.0001
end
