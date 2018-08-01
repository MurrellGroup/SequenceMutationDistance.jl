## Installation
```julia
Pkg.clone("SMD")

```

## Run tests
```julia
Pkg.test("SMD")

```

## Set paths
```julia
using SMD
```

<a id='smdfuncs-1'></a>
# smdfuncs
*`dist_matrix`* &mdash; *Function*
```julia
dist_matrix(distr1, distr2; dist_met = kmer_seeded_edit_dist)
```
distances[i, j] is distance from distr1[i] to distr2[j]

*`dist_matrix_mt`* &mdash; *Function*
```julia
dist_matrix_mt(distr1, distr2; dist_met = kmer_seeded_edit_dist)
```
dist_matrix but multithreaded instead


*`symmetric_dist_matrix`* &mdash; *Function*
```julia
symmetric_dist_matrix(distr; dist_met = kmer_seeded_edit_dist)
```
similar to dist_matrix but where distr1 == distr2. Automatically called by dist_matrix

*`symmetric_dist_matrix_mt`* &mdash; *Function*
```julia
symmetric_dist_matrix_mt(distr; dist_met = kmer_seeded_edit_dist)
```
similar to dist_matrix_mt but where distr1 == distr2. Automatically called by dist_matrix_mt


*`smd_mf`* &mdash; *Function*
```julia
smd_mf(distances_matrix::Array{Float64,2};
                freq1::Vector{Float64}=Float64[],
                freq2::Vector{Float64}=Float64[],
                unbounded_first::Bool=false,
                unbounded_second::Bool=false)
```
Mutation distance from one population to another


*`smd_distance_wrapper`* &mdash; *Function*
```julia
smd_distance_wrapper(distmat::Array{Float64,2}, inds1, inds2)
```
Wrapper fir the smd_mf call

*`fst_wrapper`* &mdash; *Function*
```julia
fst_wrapper(distmat, inds1, inds2)
```
Wrapper function

*`permutation_test`* &mdash; *Function*
```julia
permutation_test(distmat::Array{Float64,2}; l1 = nothing, l2 = nothing, tests=10000, dist_func = smd_distance_wrapper, randvariation=true)
```
