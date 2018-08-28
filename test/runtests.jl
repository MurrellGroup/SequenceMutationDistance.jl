using SequenceMutationDistance
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
	srand(42)
else
    using Test
	using LinearAlgebra
	using Random
	Random.seed!(42)
end

@testset "SMD" begin
	@testset "simple" begin
		include("test_simple.jl")
	end

	@testset "precomp" begin
		include("test_precomp.jl")
	end
end
