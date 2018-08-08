
<a id='SMD-Functions-1'></a>

# SMD Functions

- [`SMD.dist_matrix`](smdfuncs.md#SMD.dist_matrix)
- [`SMD.dist_matrix_mt`](smdfuncs.md#SMD.dist_matrix_mt)
- [`SMD.fst_wrapper`](smdfuncs.md#SMD.fst_wrapper)
- [`SMD.permutation_test`](smdfuncs.md#SMD.permutation_test)
- [`SMD.smd_distance_wrapper`](smdfuncs.md#SMD.smd_distance_wrapper)
- [`SMD.smd_mf`](smdfuncs.md#SMD.smd_mf)
- [`SMD.smd_sum`](smdfuncs.md#SMD.smd_sum)
- [`SMD.symmetric_dist_matrix`](smdfuncs.md#SMD.symmetric_dist_matrix)
- [`SMD.symmetric_dist_matrix_mt`](smdfuncs.md#SMD.symmetric_dist_matrix_mt)

<a id='SMD.dist_matrix' href='#SMD.dist_matrix'>#</a>
**`SMD.dist_matrix`** &mdash; *Function*.



```
dist_matrix(distr1, distr2; dist_met = kmer_seeded_edit_dist)
```

distances[i, j] is distance from distr1[i] to distr2[j]


<a target='_blank' href='https://github.com/MurrellGroup/SMD.jl/blob/ac382b9d2bc86379be1cc84952bffba507ea6c35/src/smdfuncs.jl#L1-L5' class='documenter-source'>source</a><br>

<a id='SMD.dist_matrix_mt' href='#SMD.dist_matrix_mt'>#</a>
**`SMD.dist_matrix_mt`** &mdash; *Function*.



```
dist_matrix_mt(distr1, distr2; dist_met = kmer_seeded_edit_dist)
```

dist_matrix but multithreaded instead


<a target='_blank' href='https://github.com/MurrellGroup/SMD.jl/blob/ac382b9d2bc86379be1cc84952bffba507ea6c35/src/smdfuncs.jl#L20-L24' class='documenter-source'>source</a><br>

<a id='SMD.symmetric_dist_matrix' href='#SMD.symmetric_dist_matrix'>#</a>
**`SMD.symmetric_dist_matrix`** &mdash; *Function*.



```
symmetric_dist_matrix(distr; dist_met = kmer_seeded_edit_dist)
```

similar to dist_matrix but where distr1 == distr2. Automatically called by dist_matrix


<a target='_blank' href='https://github.com/MurrellGroup/SMD.jl/blob/ac382b9d2bc86379be1cc84952bffba507ea6c35/src/smdfuncs.jl#L40-L44' class='documenter-source'>source</a><br>

<a id='SMD.symmetric_dist_matrix_mt' href='#SMD.symmetric_dist_matrix_mt'>#</a>
**`SMD.symmetric_dist_matrix_mt`** &mdash; *Function*.



symmetric_dist_matrix_mt(distr; dist_met = kmer_seeded_edit_dist)

similar to dist_matrix_mt but where distr1 == distr2. Automatically called by dist_matrix_mt


<a target='_blank' href='https://github.com/MurrellGroup/SMD.jl/blob/ac382b9d2bc86379be1cc84952bffba507ea6c35/src/smdfuncs.jl#L55-L59' class='documenter-source'>source</a><br>

<a id='SMD.smd_mf' href='#SMD.smd_mf'>#</a>
**`SMD.smd_mf`** &mdash; *Function*.



smd_mf(distances_matrix::Array{Float64,2};                 freq1::Vector{Float64}=Float64[],                 freq2::Vector{Float64}=Float64[],                 unbounded_first::Bool=false,                 unbounded_second::Bool=false)

Mutation distance from one population to another


<a target='_blank' href='https://github.com/MurrellGroup/SMD.jl/blob/ac382b9d2bc86379be1cc84952bffba507ea6c35/src/smdfuncs.jl#L70-L78' class='documenter-source'>source</a><br>

<a id='SMD.smd_sum' href='#SMD.smd_sum'>#</a>
**`SMD.smd_sum`** &mdash; *Function*.



smd_sum(distances_matrix::Array{Float64,2};                 freq1::Vector{Float64}=Float64[],                 freq2::Vector{Float64}=Float64[])


<a target='_blank' href='https://github.com/MurrellGroup/SMD.jl/blob/ac382b9d2bc86379be1cc84952bffba507ea6c35/src/smdfuncs.jl#L121-L125' class='documenter-source'>source</a><br>

<a id='SMD.smd_distance_wrapper' href='#SMD.smd_distance_wrapper'>#</a>
**`SMD.smd_distance_wrapper`** &mdash; *Function*.



smd_distance_wrapper(distmat::Array{Float64,2}, inds1, inds2)

Wrapper for the smd_mf call


<a target='_blank' href='https://github.com/MurrellGroup/SMD.jl/blob/ac382b9d2bc86379be1cc84952bffba507ea6c35/src/smdfuncs.jl#L133-L137' class='documenter-source'>source</a><br>

<a id='SMD.fst_wrapper' href='#SMD.fst_wrapper'>#</a>
**`SMD.fst_wrapper`** &mdash; *Function*.



fst_wrapper(distmat, inds1, inds2)

Wrapper function


<a target='_blank' href='https://github.com/MurrellGroup/SMD.jl/blob/ac382b9d2bc86379be1cc84952bffba507ea6c35/src/smdfuncs.jl#L142-L146' class='documenter-source'>source</a><br>

<a id='SMD.permutation_test' href='#SMD.permutation_test'>#</a>
**`SMD.permutation_test`** &mdash; *Function*.



permutation_test(distmat::Array{Float64,2}; l1 = nothing, l2 = nothing, tests=10000, dist_func = smd_distance_wrapper, randvariation=true)


<a target='_blank' href='https://github.com/MurrellGroup/SMD.jl/blob/ac382b9d2bc86379be1cc84952bffba507ea6c35/src/smdfuncs.jl#L156-L158' class='documenter-source'>source</a><br>

