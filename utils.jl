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


function do_some_stuff(s1::String, s2::String)
    a, b = loc_kmer_seeded_align(s1, s2)
    dst = kmer_seeded_edit_dist(degap(a), degap(b))
    for c in s1
        if c == 'N'
            dst -= 1
        end
    end
    for c in s2
        if c == 'N'
            dst -= 1
        end
    end
    return dst
end

