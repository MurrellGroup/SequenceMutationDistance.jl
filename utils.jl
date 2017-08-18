function degap(s::String)
    return replace(s, "-", "")
end

"""Return amino acid string given sequence using BioJulia"""
function translate_to_aa(s::String)
    rna = convert(RNASequence, DNASequence(s))
    return string(translate(rna))
end

"""Return sequence translated to amino acids in each reference frame
(returns three amino acid sequences)."""
function generate_aa_seqs(str::String)
    aa1 = translate_to_aa(str[1:length(str) - length(str)%3])
    aa2 = translate_to_aa(str[2:(length(str) - (length(str) - 1) % 3)])
    aa3 = translate_to_aa(str[3:(length(str) - (length(str) - 2) % 3)])
    return aa1, aa2, aa3
end


function local_edit_dist(s1::String, s2::String)
    str1, str2 = s1, s2
    if length(s1) > length(s2)
        str1, str2 = s2, s1
    end
    a, b = loc_kmer_seeded_align(str1, str2)
    dst = 0
    for i in 1:length(a)
        if a[i] != b[i] && a[i] != 'N' && b[i] != 'N'
            dst += 1
        end
    end
    return dst
end

#-------PHREDS AND PROBS--------
const Phred = Int8
const Prob = Float64
const LogProb = Float64

const MIN_PHRED = Phred(1)
const MAX_PHRED = Phred(Int('~') - 33)

@generated function phred_to_log_p(x)
    return quote
        return x / (-10.0)
    end
end

function phred_to_p(q::Phred)
    return exp10(phred_to_log_p(q))
end

function phred_to_p(x::Vector{Phred})
    return exp10.(phred_to_log_p(x))
end

function p_to_phred(p::Prob)
    return Phred(min(round(-10.0 * log10(p)), MAX_PHRED))
end

function p_to_phred(x::Vector{LogProb})
    return Phred[p_to_phred(p) for p in x]
end

#-------FILTER FUNCTION----------

function quality_filter(infile,outfile=join(split(infile, ".")[1:end-1], ".") * ".filt.fasta" ; errorRate=0.01,minLength=0,labelPrefix="seq",errorOut = true)
    seqs, scores, names = read_fastq(infile)
    
    inds = quality_filter_inds(seqs, scores, errorRate=errorRate, minLength=minLength)
        
    #names = names[inds]
    if errorOut == true
        names = ["$labelPrefix$(i)|ee=$(mean(p_vals[i]))" for i in inds]
    else
        names = ["$labelPrefix$(i)" for i in 1:length(inds)]
    end
    
    if last(outfile) == 'a'
        write_fasta(outfile, seqs[inds], names=names)
    else
        write_fastq(outfile, seqs[inds], scores[inds], names=names)
    end
end

function quality_filter(seqs::Array{String, 1}, scores::Array, names::Array{String, 1}; errorRate=0.01,minLength=0)
    inds = quality_filter_inds(seqs, scores, errorRate=errorRate, minLength=minLength)
    return seqs[inds], scores[inds], names[inds]
end

function quality_filter_inds(seqs::Array{String, 1}, scores::Array; errorRate=0.01,minLength=100)
    p_vals = [phred_to_p(score) for score in scores]
    inds = filter!(x -> (mean(p_vals[x]) < errorRate && (length(seqs[x]) >= minLength)) , collect(1:length(seqs)))
    return inds
end

#------KMERS-------

const NUCLEOTIDE_BITS = Dict('A' => unsigned(0),
                       'C' => unsigned(1),
                       'G' => unsigned(2),
                       'T' => unsigned(3))

const KmerType = Array{UInt32, 1}

"""Count kmers in string"""
function kmer_count(str::String, k::Int)
    # TODO: could directly encode `str` as 2-bit BioSequence
    bins = zeros(eltype(KmerType), 4^k)

    mask = unsigned(4^k - 1)  # all ones
    kmer = unsigned(0)
    for c in str[1:k-1]
        kmer = (kmer << 2) + NUCLEOTIDE_BITS[c]
    end
    for c in str[k:end]
        kmer = ((kmer << 2) & mask) + NUCLEOTIDE_BITS[c]
        bins[kmer + 1] += 1
    end
    return bins
end

"""Compute distance function that is correct for small differences."""
function corrected_kmer_dist(kmers1::KmerType, kmers2::KmerType; k = nothing)
    if k == nothing
        k = Int(log(4, length(kmers1)))
    end
    return sqeuclidean(kmers1, kmers2)/ (k*(sum(kmers1) + sum(kmers2)))
end

function sorted_freqs(vec)
    propDict = proportionmap(vec)
    seqkeys = keys(propDict)
    return reverse(sort([(propDict[k],k) for k in seqkeys]))
end
