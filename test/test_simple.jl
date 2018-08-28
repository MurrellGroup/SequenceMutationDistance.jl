@testset "no-freq" begin
	s = 100
	
	matrix = ones(s,s)
	@test isapprox(smd_mf(matrix)[1], 1)
	
	matrix = Array{Float64,2}(I,s,s)
	@test isapprox(smd_mf(matrix)[1], 0)
	
	matrix = ones(s,s) .- Array{Float64,2}(I,s,s)
	@test isapprox(smd_mf(matrix)[1], 0)
	@test all(isapprox(smd_mf(matrix)[2],
		(1/s) .* Array{Float64,2}(I,s,s)	
		))
end
