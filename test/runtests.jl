using Test
using IntervalSets
using BasicBSpline
using ElasticSurfaceEmbedding
using Images

dir_result_a = joinpath(@__DIR__, "result_a")
dir_result_b = joinpath(@__DIR__, "result_b")

rm(dir_result_b, recursive=true, force=true)
config_dir(dir_result_b)

@testset "ElasticSurfaceEmbedding.jl" begin
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
end
