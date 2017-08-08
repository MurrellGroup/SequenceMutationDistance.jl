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

