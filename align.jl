"""Needleman-Wunch, with end gaps penalized slightly less."""
function nw_align(s1::String, s2::String; edge_reduction = 0.99)
    # edge_reduction is a multiplicative score that gets multiplied to
    # the penalties along the edges, to prefer terminal gaps.
    ins_cost = -1.0  # -> : horizontal
    del_cost = -1.0  # V : vertical
    mismatch_cost = -1.0
    match_cost = 1.0

    s1arr = collect(s1)  # vertical
    s2arr = collect(s2)  # horizontal
    arr = zeros(length(s1arr)+1, length(s2arr)+1)
    traceArr = zeros(Int, length(s1arr), length(s2arr))

    # this will need to be generalized when we want to allow overhang
    arr[:, 1] = edge_reduction*del_cost*(0:length(s1arr))
    arr[1,:] = edge_reduction*ins_cost*(0:length(s2arr))

    for i in 2:length(s1arr)+1
        for j in 2:length(s2arr)+1
            if s1arr[i-1] == s2arr[j-1]
                diag = arr[i-1, j-1] + match_cost
            else
                diag = arr[i-1, j-1] + mismatch_cost
            end

            # to handle the lower edge penalties.
            delMult = 1
            if i == length(s1arr)+1
                delMult = edge_reduction
            end
            insMult=1
            if j == length(s2arr)+1
                insMult = edge_reduction
            end

            ins = arr[i-1, j]+(ins_cost*insMult)
            del = arr[i, j-1]+(del_cost*delMult)
            scores = [diag, del, ins]
            best = indmax(scores)
            if best == 1
                arr[i, j] = diag
            elseif best ==2
                arr[i, j] = del
            elseif best ==3
                arr[i, j] = ins
            end
            traceArr[i-1, j-1] = best
        end
    end
    alignedScore = arr[end, end]

    # return arr

    # the trace endpoing will need to be generalized if we want to
    # allow overhang.
    trI, trJ = length(s1arr), length(s2arr)

    # First compute the trace running backwards. Initialized to the
    # maximum size. With unit penalties, is it possible to tell how
    # big this should be in advance?
    backtrace = Array{Int}(length(s1arr) + length(s2arr))
    # If you hit any boundary, you have to run all the way to the side.
    btInd = 1
    while (trI > 0) && (trJ > 0)
        backtrace[btInd] = traceArr[trI, trJ]
        if backtrace[btInd] == 1
            trI += -1
            trJ += -1
        elseif backtrace[btInd] == 2
            trJ += -1
        elseif backtrace[btInd] == 3
            trI += -1
        end
        btInd += 1
    end

    # If you hit the boundaries not at the top left corner.
    while trI > 0
        backtrace[btInd] = 3
        btInd += 1
        trI += -1
    end

    while trJ > 0
        backtrace[btInd] = 2
        btInd += 1
        trJ += -1
    end

    backtrace = backtrace[1:btInd-1]

    # This will need to be generalized to work on non-strings. Not
    # important.
    ali1arr = Array{Char}(length(backtrace))
    ali2arr = Array{Char}(length(backtrace))
    ind1 = 1
    ind2 = 1
    for i in 1:length(backtrace)
        tr = backtrace[length(backtrace)-(i-1)]
        if tr == 1
            ali1arr[i] = s1arr[ind1]
            ali2arr[i] = s2arr[ind2]
            ind1 += 1
            ind2 += 1
        elseif tr == 2
            ali2arr[i] = s2arr[ind2]
            ali1arr[i] = '-'
            ind2 += 1
        elseif tr == 3
            ali1arr[i] = s1arr[ind1]
            ali2arr[i] = '-'
            ind1 += 1
        end
    end
    return join(ali1arr), join(ali2arr)
end


#--------Banded alignment---------

"""Wrapper for nw_align and banded_nw_align."""
function nw_align(s1::String, s2::String, banded::Float64)
    if banded <= 0
        return nw_align(s1, s2)
    else
        return banded_nw_align(s1, s2, band_coeff = banded)
    end
end

"""Set value in band where `i` and `j` are in square matrix coords, `dim_diff` = ncols - nrows"""
function add_to_band!(band, val, i::Int, j::Int, bandwidth::Int, dim_diff::Int)
    # smallest coord is slice, largest is +/- offset from middle of slice + offset for diff in string lengths
    ii = min(i, j)
    jj = j - i + bandwidth + 1 + max(0, dim_diff)
    band[ii, jj] = val
end

"""Get value from band where `i` and `j` are in square matrix coords, `dim_diff` = ncols - nrows"""
function get_band_val(band, i::Int, j::Int, bandwidth::Int, dim_diff::Int)
    # smallest coord is slice, largest is +/- offset from middle of slice + offset for diff in string lengths
    ii = min(i, j)
    jj = j - i + bandwidth + 1 + max(0, dim_diff)
    return band[ii, jj]
end

"""Check if within width of band (ignore length)"""
function in_band(i, j, bandwidth, dim_diff)
    ii = min(i, j)
    jj = j - i + bandwidth + 1 + max(0, dim_diff)
    return !( ii < 1 || jj < 1 || jj > 2*bandwidth + 1 + abs(dim_diff) )
end


"""
Like nw_align, but sub quadratic by only computing values within a band around the center diagonal.
One "band" of radius 3 = (4,1), (3,1), (2,1), (1,1), (1,2), (1,3), (1,4), aka upside-down L shape.
band_coeff = 1 is sufficient to get same alignments as nw_align for 10% diverged sequences ~97% of the time;
increase for more conservative alignment with longer computation time.
Radius of band = `bandwidth` = `band_coeff` * sqrt(avg seq length)
"""
function banded_nw_align(s1::String, s2::String; edge_reduction = 0.99, band_coeff = 1)
    # calculate band width
    avg_len = (length(s1) + length(s2)) / 2
    bandwidth = Int(ceil(band_coeff * sqrt(avg_len)))
    # edge_reduction is a multiplicative score that gets multiplied to
    # the penalties along the edges, to prefer terminal gaps.
    ins_cost = -1.0  # -> : horizontal
    del_cost = -1.0  # V : vertical
    mismatch_cost = -1.0
    match_cost = 1.0
    
    s1arr = collect(s1)  # vertical
    s2arr = collect(s2)  # horizontal
    
    # dimension difference added to band length on one side
    # to give padding for large sequence length differences
    dim_diff = length(s1arr) - length(s2arr)
    maxlen = max(length(s1), length(s2))
    arr = zeros(maxlen + 1, 2*bandwidth + 1 + abs(dim_diff))
    traceArr = zeros(Int, size(arr))
    for i in 1:bandwidth+1
        add_to_band!(arr, edge_reduction*del_cost*(i-1), i, 1, bandwidth, dim_diff)
        add_to_band!(arr, edge_reduction*ins_cost*(i-1), 1, i, bandwidth, dim_diff)
    end

    for i in 2:length(s1arr)+1
        # keep j within 2*bandwidth + 1 band of i, but keep bandwidth sized buffer on edges
        lo = 0
        # to keep a big enough padding around the edges when there's a big length mismatch
        if i <= bandwidth + 1 + max(0, dim_diff)
            lo = 2
        else
            lo = i - bandwidth - max(0, dim_diff)
        end
        for j in lo:length(s2arr)+1
            # keep j within band
            if !(in_band(i, j, bandwidth, dim_diff))
                break
            end
            
            if in_band(i-1, j-1, bandwidth, dim_diff) && s1arr[i-1] == s2arr[j-1]
                diag = get_band_val(arr, i-1, j-1, bandwidth, dim_diff) + match_cost
            else
                diag = get_band_val(arr, i-1, j-1, bandwidth, dim_diff) + mismatch_cost
            end
            
            # to handle the lower edge penalties.
            delMult = 1
            if i == length(s1arr)+1
                delMult = edge_reduction
            end
            insMult=1
            if j == length(s2arr)+1
                insMult = edge_reduction
            end

            ins = in_band(i-1, j, bandwidth, dim_diff) ? 
                    (get_band_val(arr, i-1, j, bandwidth, dim_diff)+(ins_cost*insMult)) : NaN
            del = in_band(i, j-1, bandwidth, dim_diff) ? 
                    (get_band_val(arr, i, j-1, bandwidth, dim_diff)+(del_cost*delMult)) : NaN

            scores = [diag, del, ins]
            best = indmax(scores)
            add_to_band!(arr, scores[best], i, j, bandwidth, dim_diff)
            add_to_band!(traceArr, best, i-1, j-1, bandwidth, dim_diff)                                         
        end
    end
    alignedScore = get_band_val(arr, length(s1arr)+1, length(s2arr)+1, bandwidth, dim_diff)

    trI, trJ = length(s1arr), length(s2arr)
    backtrace = Array{Int}(length(s1arr) + length(s2arr))
    btInd = 1
    while (trI > 0) && (trJ > 0)
        backtrace[btInd] = get_band_val(traceArr, trI, trJ, bandwidth, dim_diff)
        if backtrace[btInd] == 1
            trI += -1
            trJ += -1
        elseif backtrace[btInd] == 2
            trJ += -1
        elseif backtrace[btInd] == 3
            trI += -1
        else 
            error("Bad trace value: $(backtrace[btInd]): ($trI, $trJ)")
        end
        btInd += 1
    end
    
    # If you hit the boundaries not at the top left corner.
    while trI > 0
        backtrace[btInd] = 3                                                    
        btInd += 1
        trI += -1                                                               
    end                                                                         
    while trJ > 0
        backtrace[btInd] = 2                                                    
        btInd += 1                                                              
        trJ += -1                                                               
    end                                                                         
                                                                                
    backtrace = backtrace[1:btInd-1]           
    ali1arr = Array{Char}(length(backtrace))                                    
    ali2arr = Array{Char}(length(backtrace))                                    
    ind1 = 1                                                                    
    ind2 = 1                                                                    
    for i in 1:length(backtrace)                                                
        tr = backtrace[length(backtrace)-(i-1)]                                 
        if tr == 1                                                              
            ali1arr[i] = s1arr[ind1]                                            
            ali2arr[i] = s2arr[ind2]                                            
            ind1 += 1                                                           
            ind2 += 1                                                           
        elseif tr == 2                                                          
            ali2arr[i] = s2arr[ind2]                                            
            ali1arr[i] = '-'                                                    
            ind2 += 1                                                           
        elseif tr == 3                                                          
            ali1arr[i] = s1arr[ind1]                                            
            ali2arr[i] = '-'                                                    
            ind1 += 1                                                           
        end                                                                     
    end                                                                         
    return join(ali1arr), join(ali2arr)                                         
end

#--------Kmer internals---------

function unique_key(dicto::Dict{String, Int}, keyo::String, indo::Int)
    if haskey(dicto, keyo)
        dicto[keyo] = -1
    else
        dicto[keyo] = indo
    end
end

function sorted_matches(s1, s2, wordlength, skip, aligncodons)
    word_dict1 = Dict{String, Int}()
    word_dict2 = Dict{String, Int}()
    # we can space one of these out, but not both
    for i in 1:skip:(length(s1)-(wordlength-1))  # bounds checked
        if aligncodons
            unique_key(word_dict1, translate_to_aa(s1[i:i+(wordlength-1)]), i)
        else
            unique_key(word_dict1, s1[i:i+(wordlength-1)], i)  # bounds checked
        end
    end
    for i in 1:(length(s2)-(wordlength-1))  # bounds checked
        if aligncodons
            unique_key(word_dict2, translate_to_aa(s2[i:i+(wordlength-1)]), i)
        else
            unique_key(word_dict2, s2[i:i+(wordlength-1)], i)  # bounds checked
        end
    end
    
    intersection = Dict{String, UInt8}()
    for (word, ind) in word_dict1
        if ind > 0 && haskey(word_dict2, word) && word_dict2[word] > 0
            intersection[word] = 0
        end
    end
    common = collect(keys(intersection))
    matches = zeros(Int, length(common), 2)       
    for i in 1:length(common)
        matches[i, 1] = word_dict1[common[i]]
        matches[i, 2] = word_dict2[common[i]]
    end
    return sortrows(matches, by=x->(x[1]))
end

"""Create sorted list of matches of amino acid kmers in any reference frame."""
function sorted_aa_matches(str1, str2, wordlength)
    word_dict1 = Dict{String, Int}()
    word_dict2 = Dict{String, Int}()
    # get amino acid sequence in each reference frame
    s1s = generate_aa_seqs(str1)
    s2s = generate_aa_seqs(str2)
    for (offset, s1) in enumerate(s1s)
        for i in 1:(length(s1) - (div(wordlength, 3) - 1))
            # offset acounts for the shift in index due to being in a different reference frame.
            unique_key(word_dict1, s1[i:i+(div(wordlength, 3) - 1)], i*3 + offset - 3)
        end
    end
    for (offset, s2) in enumerate(s2s)
        for i in 1:(length(s2) - (div(wordlength, 3) - 1))
            unique_key(word_dict2, s2[i:i+(div(wordlength, 3) - 1)], i*3 + offset - 3)
        end
    end    
    intersection = Dict{String, UInt8}()
    for (word, ind) in word_dict1
        if ind > 0 && haskey(word_dict2, word) && word_dict2[word] > 0
            intersection[word] = 0
        end
    end
    common = collect(keys(intersection))
    matches = zeros(Int, length(common), 2)       
    for i in 1:length(common)
        matches[i, 1] = word_dict1[common[i]]
        matches[i, 2] = word_dict2[common[i]]
    end
    return sortrows(matches, by=x->(x[1]))
end


function matches_are_inconsistent(matches)
    # `matches` must already be sorted
    if size(matches)[1] <= 1
        return false
    end
    return minimum(matches[2:length(matches[:, 2]), 2] -
                   matches[1:(length(matches[:, 2]) - 1), 2]) < 1
end


"""Some edge cases commonly arise where a kmer match starts before
the previous kmer match ends, but the two sequences still
mismatch.

To handle these, we need to make sure that all runs of kmer
matches are spaced out by exactly `skip`, or by an amount that
is greater than wordlength. So I just check for this and scrub
these from `sortedmatches`

Kmer matches in second sequence may not overlap at all.

"""
function clean_matches(matches, wordlength, skip)
    # `matches` must already by sorted
    cleanedmatches = Vector{Int}[]
    push!(cleanedmatches, matches[1,:])
    currentvec = matches[1,:]
    for i in 2:length(matches[:, 2])
        # The condition below involves a kludge. It should really be
        # considering the case that skip == diff, rather than skip <
        # diff. Fix later.
        if !(skip < matches[i, 1] - currentvec[1] < wordlength ||
             matches[i, 2] - currentvec[2] < wordlength)
            push!(cleanedmatches, matches[i,:])
            currentvec = matches[i,:]
        end
    end
    return hcat(cleanedmatches...)'
end

function merge_overlapping(matches, wordlength, skip)
    range_inds = Vector{Int}[]
    push!(range_inds,[1, 1])
    current_range_ind_ind = 1
    for i in 2:length(matches[:, 2])
        # TODO: THESE BOUNDS (especially +wordlength) NEED TO BE CHECKED.
        # if ((matches[i, 1] < matches[i-1, 1] + (wordlength - 1) &&
        #      (matches[i, 2] < matches[i-1, 2] + (wordlength - 1))))
        if (matches[i, 1] - matches[i-1, 1] == skip) &&
            (matches[i, 2] - matches[i-1, 2] == skip)
            range_inds[current_range_ind_ind][2] = i
        else
            push!(range_inds,[i, i])
            current_range_ind_ind += 1
        end
    end
    return range_inds
end

# get matches
function get_matches(s1, s2, clean, range_inds, wordlength)
    s1matches = [clean[range_inds[i], 1] + [0, wordlength - 1]
                 for i in 1:length(range_inds)]
    s2matches = [clean[range_inds[i], 2] + [0, wordlength - 1]
                 for i in 1:length(range_inds)]
    match1 = [s1[i[1]:i[2]] for i in s1matches]
    match2 = [s2[i[1]:i[2]] for i in s2matches]
    return match1, match2
end

function get_mismatches(s1, s2, matches, range_inds, wordlength)
    s1mismatches = [matches[(range_inds[i-1][2]):range_inds[i][1], 1] .+ [wordlength, -1]
                    for i in 2:length(range_inds)]
    s2mismatches = [matches[(range_inds[i-1][2]):range_inds[i][1], 2] .+ [wordlength,-1]
                    for i in 2:length(range_inds)]
    s1mismatches = vcat([[1, matches[1, 1]-1]],
                        s1mismatches,
                        [[matches[range_inds[length(range_inds)][2], 1]+wordlength, length(s1)]])
    s2mismatches = vcat([[1, matches[1, 2]-1]],
                        s2mismatches,
                        [[matches[range_inds[length(range_inds)][2], 2]+wordlength, length(s2)]])

    mismatch1 = [s1[i[1]:i[2]] for i in s1mismatches]
    mismatch2 = [s2[i[1]:i[2]] for i in s2mismatches]
    return mismatch1, mismatch2
end

"""Find longest increasing subsequence of second column of given array.
Used to resolve bad orders of word matches while preserving as many matches as possible"""
function longest_incr_subseq(arr::Array{Int, 2})
    len = size(arr)[1]
    prevs = ones(Int, len)
    curr = ones(Int, len)
    maxlen = 0
    for i in 1:len
        # binary search for longest seq so far that's less than current value.
        lo = 1
        hi = maxlen
        while lo <= hi
            mid = Int(ceil((lo+hi)/2))
            if arr[curr[mid], 2] < arr[i, 2]
                lo = mid + 1
            else
                hi = mid - 1
            end
        end
        # update sequences seen and corresponding previous indices
        if lo > 1
            prevs[i] = curr[lo-1]
        else
            prevs[i] = 1
        end
        curr[lo] = i
        if lo > maxlen
            maxlen = lo
        end
    end
    # build subsequence from indices
    subseq = zeros(Int, maxlen, 2)
    k = curr[maxlen]
    for i in maxlen:-1:1
        subseq[i, :] = arr[k, :]
        k = prevs[k]
    end
    return subseq
end


"""returns aligned strings"""
function kmer_seeded_align(s1::String, s2::String;
                           wordlength = 30,
                           skip = 10,
                           aligncodons = false,
                           banded = 1.0,
                           debug::Bool = false)
    # ToDo:
    # 1: Make this recurse. So instead fo calling nw_align on the set
    # of mismatches strings, it calls itself, but with a lower word
    # length. This will help for really noisy sequences.

    if (length(s1) < wordlength || length(s2) < wordlength || s1 == "" || s2 == "")
        return nw_align(s1, s2)
    end

    sorted = sorted_matches(s1, s2, wordlength, skip, aligncodons)
    if size(sorted)[1] == 0
        return nw_align(s1, s2, banded)
    end
    if matches_are_inconsistent(sorted)
        # try to extract an in-order subsequence to resolve inconsistency
        sorted = longest_incr_subseq(sorted)
        if matches_are_inconsistent(sorted)
            println("Notice: Word matching produced inconsistent ordering.",
                    " Consider a larger word size. Returning full DP alignment.",
                    "\nNote: this message shouldn't print.")
            return nw_align(s1, s2, banded)
        end
    end

    clean = clean_matches(sorted, wordlength, skip)
    range_inds = merge_overlapping(clean, wordlength, skip)
    # get matches
    match1, match2 = get_matches(s1, s2, clean, range_inds, wordlength)
    # get mismatches
    mismatch1, mismatch2 = get_mismatches(s1, s2, clean, range_inds, wordlength)

    alignedStrings = collect(nw_align(mismatch1[1], mismatch2[1]))
    for i in 1:(length(match1))
        alignedStrings[1] = alignedStrings[1] * match1[i]
        alignedStrings[2] = alignedStrings[2] * match2[i]
        # al1, al2 = align(mismatch1[i+1], mismatch2[i+1],-1.0,-1.0, 20)
        al1, al2 = nw_align(mismatch1[i+1], mismatch2[i+1])
        alignedStrings[1] = alignedStrings[1] * al1
        alignedStrings[2] = alignedStrings[2] * al2
    end

    if debug
        if (degap(alignedStrings[1]) != degap(s1) || degap(alignedStrings[2]) != degap(s2))
            error("Aligned strings do not match original strings")
        end
    end
    return alignedStrings
end

"""If aa_matches = true, will attempt to find amino acid matches in any reference frame, 
and add the nucleotide Hamming distance of these matches to Levenshtein distances of mismatches."""
function kmer_seeded_edit_dist(s1::String , s2::String;
                               wordlength = 30,
                               skip = 5,
                               aa_matches = false)
    if (length(s1) < wordlength || length(s2) < wordlength || s1 == "" || s2 == "")
        return levenshtein(s1, s2)
    end
    if aa_matches
        sorted = sorted_aa_matches(s1, s2, wordlength)
        skip = 1
    else
        sorted = sorted_matches(s1, s2, wordlength, skip, false)
    end
    if size(sorted)[1] == 0
        return levenshtein(s1, s2)
    end
    if matches_are_inconsistent(sorted)
        # try to extract an in-order subsequence to resolve inconsistency
        sorted = longest_incr_subseq(sorted)
        if matches_are_inconsistent(sorted)
            println("Notice: Word matching produced inconsistent ordering.",
                    " Consider a larger word size. Returning full edit distance.")
            return levenshtein(s1, s2)
        end
    end
    trim = 0
    if aa_matches
        # trim word matches to avoid double counting some errors on ends of matches
        trim = max(div(wordlength, 5), 1)
        for i in 1:size(sorted)[1]
            sorted[i, 1] += trim
            sorted[i, 2] += trim
        end
    end
    clean = clean_matches(sorted, wordlength, skip)
    range_inds = merge_overlapping(clean, wordlength, skip)
    # 2*trim to accomodate trimmed word matches on each side of the word match
    mismatch1, mismatch2 = get_mismatches(s1, s2, clean, range_inds, wordlength - 2*trim)
    matched_diffs = 0
    if aa_matches
        match1, match2 = get_matches(s1, s2, clean, range_inds, wordlength - 2*trim)
        for i in 1:length(match1)
            for j in 1:length(match1[i])
                if match1[i][j] != match2[i][j]
                    matched_diffs += 1
                end
            end
        end
    end
    return matched_diffs + sum([levenshtein(mismatch1[i], mismatch2[i]) for i in 1:length(mismatch1)])
end

"""Wrapper for kmer_seeded_edit_dist for aligning amino acids"""
function aa_kmer_seeded_edit_dist(s1::String, s2::String; wordlength = 30, skip = 5)
    return kmer_seeded_edit_dist(s1, s2, wordlength=wordlength, skip=skip, aa_matches=true)
end
