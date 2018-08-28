mat = rand(100,100)
f1 = rand(100)
f2 = rand(100)

@testset "no_freq" begin
	@test isapprox(smd_mf(mat)[1], 0.013642842651339401)
end

@testset "freq" begin
	@test isapprox(smd_mf(mat, freq1=f1)[1], 0.017396290381031448)
	@test isapprox(smd_mf(mat, freq2=f2)[1], 0.01785180369745115)
	@test isapprox(smd_mf(mat, freq1=f1, freq2=f2)[1], 0.020111466838585175)
end

@testset "freq-unbounded" begin
	@test isapprox(smd_mf(mat, freq1=f1, freq2=f2, unbounded_first = true)[1], 0.008888605782330361)
	@test isapprox(smd_mf(mat, freq1=f1, freq2=f2, unbounded_second = true)[1], 0.008261507631168017)
	@test isapprox(smd_mf(mat, freq1=f1, freq2=f2, unbounded_first = true, unbounded_second = true)[1], 0.0)
end
