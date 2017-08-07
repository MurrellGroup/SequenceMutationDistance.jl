"""distances[i, j] is distance from distr1[i] to distr2[j]"""
function default_dist_matrix(distr1, distr2; k=6)
    # create distance matrix using distmat function
    distances = zeros(length(distr1), length(distr2))
    for i in 1:length(distr1)
        for j in 1:length(distr2)
            distances[i, j] = kmer_seeded_edit_dist(distr1[i], distr2[j])
        end
    end
    return distances
end

"""Mutation distance from one population to another"""
function smd_mf(distances_matrix::Array{Float64,2};
                freq1::Vector{Float64}=Float64[],
                freq2::Vector{Float64}=Float64[],
                unbounded_first::Bool=false,
                unbounded_second::Bool=false)
    if (size(distances_matrix)[1] != length(freq1))
        freq1 = ones(size(distances_matrix)[1])
    end
    if (size(distances_matrix)[2] != length(freq2))
        freq2 = ones(size(distances_matrix)[2])
    end

    # init solver
    m = Model(solver = ClpSolver(SolveType=5))

    # flow is a matrix with how much of each distr1 elem maps to how
    # much of each distr2 elem
    @variable(m, flow[1:length(freq1), 1:length(freq2)] >= 0 )

    # constraint: The occurrence of a element "i" in flow is equal to
    # its frequency (normalized)
    if !unbounded_first
        for i = 1:length(freq1)
            @constraint(m, sum(flow[i, :]) == freq1[i] / sum(freq1))
        end
    end

    # constraint: The occurence of a element "j" in flow is equal to
    # its frequency (normalized)
    if !unbounded_second
        for j = 1:length(freq2)
            @constraint(m, sum(flow[:, j]) == freq2[j] / sum(freq2))
        end
    end

    # minimize sum of distance flows
    @objective(m, Min, sum(flow .* distances_matrix))

    status = solve(m)
    return (getobjectivevalue(m), (getvalue(flow)))
end


