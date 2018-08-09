
<a id='SMD-Functions-1'></a>

# SMD Functions

- [`SMD.dist_matrix`](smdfuncs.md#SMD.dist_matrix)
- [`SMD.smd_mf`](smdfuncs.md#SMD.smd_mf)

<a id='SMD.dist_matrix' href='#SMD.dist_matrix'>#</a>
**`SMD.dist_matrix`** &mdash; *Function*.



```
dist_matrix(distr1, distr2; dist_met = kmer_seeded_edit_dist)
```

distances[i, j] is distance from distr1[i] to distr2[j]


<a target='_blank' href='https://github.com/MurrellGroup/SMD.jl/blob/ed7e08f87219aa3039e88980799fcaafb61ef954/src/smdfuncs.jl#L1-L5' class='documenter-source'>source</a><br>

<a id='SMD.smd_mf' href='#SMD.smd_mf'>#</a>
**`SMD.smd_mf`** &mdash; *Function*.



smd_mf(distances_matrix::Array{Float64,2};                 freq1::Vector{Float64}=Float64[],                 freq2::Vector{Float64}=Float64[],                 unbounded_first::Bool=false,                 unbounded_second::Bool=false)

Mutation distance from one population to another


<a target='_blank' href='https://github.com/MurrellGroup/SMD.jl/blob/ed7e08f87219aa3039e88980799fcaafb61ef954/src/smdfuncs.jl#L70-L78' class='documenter-source'>source</a><br>

