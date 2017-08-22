"""Denoising algorithm designed for Illumina data"""
function menoise2(fl::String; error_rate = 0.005, out = fl[1:end-6]*"-menoiseed.fasta", fastmode = false)
    println("Reading from ", fl)
    if fl[end-5:end] == ".fasta"
        seqs = read_fasta(fl, seqtype=String)
    elseif fl[end-5:end] == ".fastq"
        seqs, _, _ = read_fastq(fl, seqtype=String)
    else
        error("bad file type: $fl -- must be type .fasta or .fastq")
    end
    println("Read ", length(seqs), " seqs")
    
    numseqs = length(seqs)
    cutoff_freq = 5
    cutoff_freq /= numseqs
    k = 5
    centroids = Array{Any, 1}[]
    not_centroids = Array{Any, 1}[]
    freqs = proportionmap(seqs)
    for (s, f) in freqs
        if f >= cutoff_freq
            push!(centroids, [s, f])
        elseif fastmode == false
            push!(not_centroids, [s, f])
        end
    end
    centrkmers = [kmer_count(q[1], k) for q in centroids]
    for (s, f) in not_centroids
        skc = kmer_count(s, k)
        dists = map(q -> corrected_kmer_dist(skc, q, k=k), centrkmers)
        #dists = map(q -> kmer_seeded_edit_dist(s, q[1]), centroids)
        closest = indmin(dists)
        centroids[closest][2] += f
    end

    # determine actual cutoff
    seqfreqs = reverse(sort(centroids, by=x->x[2]))
    toplen = length(seqfreqs[1][1])
    numoccurs = seqfreqs[1][2] *numseqs
    println("numoccurs: ", numoccurs)
    numchange = (1 - ((1 - error_rate) ^ toplen)) * numoccurs
    println("numchange: ", numchange)
    lambda = numchange / (3 * toplen)
    cutoff_freq = 1
    while (lambda^cutoff_freq) * exp(-cutoff_freq) / factorial(cutoff_freq) > 0.0001
        cutoff_freq += 1
    end
    println("Using cutoff of $cutoff_freq occurrences")

    cutoff_freq /= numseqs
    k = 5
    centroids = Array{Any, 1}[]
    not_centroids = Array{Any, 1}[]
    freqs = proportionmap(seqs)
    for (s, f) in freqs
        if f >= cutoff_freq
            push!(centroids, [s, f])
        elseif fastmode == false
            push!(not_centroids, [s, f])
        end
    end
    if fastmode == false
        centrkmers = [kmer_count(q[1], k) for q in centroids]
        for (s, f) in not_centroids
            skc = kmer_count(s, k)
            dists = map(q -> corrected_kmer_dist(skc, q, k=k), centrkmers)
            #dists = map(q -> kmer_seeded_edit_dist(s, q[1]), centroids)
            closest = indmin(dists)
            centroids[closest][2] += f
        end
    end
    
    println("Found ", length(centroids), " centroids")
    println("Writing to ", out)
    write_fasta(out, [s[1] for s in centroids],
        names=["seq"*string(i)*"_"*string(Int(round(f[2] * numseqs))) 
            for (i, f) in enumerate(centroids)])
end

"""Denoising algorithm designed for Illumina data"""
function menoise(fl::String; error_rate = 0.005, out = fl[1:end-6]*"-menoiseed.fasta", fastmode = false)
    println("Reading from ", fl)
    if fl[end-5:end] == ".fasta"
        seqs = read_fasta(fl, seqtype=String)
    elseif fl[end-5:end] == ".fastq"
        seqs, _, _ = read_fastq(fl, seqtype=String)
    else
        error("bad file type: $fl -- must be type .fasta or .fastq")
    end
    println("Read ", length(seqs), " seqs")
    
    # determine actual cutoff
    numseqs = length(seqs)
    seqfreqs = sorted_freqs(seqs)
    toplen = length(seqfreqs[1][2])
    numstay = seqfreqs[1][1] *numseqs
    numchange = numstay / ((1 - error_rate) ^ toplen) - numstay
    println("numchange: ", numchange)
    lambda = numchange / (3 * toplen)
    cutoff_freq = 1
    while (lambda^cutoff_freq) * exp(-cutoff_freq) / factorial(cutoff_freq) > 0.0001
        cutoff_freq += 1
    end
    println("Using cutoff of $cutoff_freq occurrences")
    
    # account for weird observed distributions
    cutoff_freq = min(numstay, cutoff_freq)
    cutoff_freq /= numseqs
    k = 5
    centroids = Array{Any, 1}[]
    not_centroids = Array{Any, 1}[]
    freqs = proportionmap(seqs)
    for (s, f) in freqs
        if f >= cutoff_freq
            push!(centroids, [s, f])
        elseif fastmode == false
            push!(not_centroids, [s, f])
        end
    end
    if fastmode == false
        centrkmers = [kmer_count(q[1], k) for q in centroids]
        for (s, f) in not_centroids
            skc = kmer_count(s, k)
            dists = map(q -> corrected_kmer_dist(skc, q, k=k), centrkmers)
            #dists = map(q -> kmer_seeded_edit_dist(s, q[1]), centroids)
            closest = indmin(dists)
            centroids[closest][2] += f
        end
    end
    
    println("Found ", length(centroids), " centroids")
    println("Writing to ", out)
    write_fasta(out, [s[1] for s in centroids],
        names=["seq"*string(i)*"_"*string(Int(round(f[2] * numseqs))) 
            for (i, f) in enumerate(centroids)])
end

"""Denoising algorithm designed for Illumina data"""
function menoise3(fl::String; error_rate = 0.005, out = fl[1:end-6]*"-menoiseed.fasta", fastmode = false)
    println("Reading from ", fl)
    if fl[end-5:end] == ".fasta"
        seqs = read_fasta(fl, seqtype=String)
    elseif fl[end-5:end] == ".fastq"
        seqs, scores, _ = read_fastq(fl, seqtype=String)
    else
        error("bad file type: $fl -- must be type .fasta or .fastq")
    end
    println("Read ", length(seqs), " seqs")
    
    #change domain of scores
    scores = [phred_to_p(x) for x in scores]

    # determine actual cutoff
    numseqs = length(seqs)
    seqfreqs = sorted_freqs(seqs)
    toplen = length(seqfreqs[1][2])
    error_rates= [scores[i] for (i,x) in enumerate(seqs) if x == seqfreqs[1][2]]
    errorProbs = [mean([x[i] for x in error_rates]) for i in 1:toplen]
    numstay = seqfreqs[1][1] *numseqs
    numchange = numstay / ((1 - error_rate) ^ toplen) - numstay
    println("numchange: ", numchange)

    totProbs = sum(errorProbs);
    binSize = numstay/e^(-totProbs)
    println("binsize: ", binSize)

    current = maximum(errorProbs);
    sumOthers = totProbs - current;
    probOfNoOtherMutations = e^(-sumOthers);
    probOfOnlyCurrentMut = (current/3)*probOfNoOtherMutations;
    #0.577 is euler macaroni. formula is some black box by ben
    SDs = log(numseqs^2 / (2 * pi * log(numseqs^2/(2 * pi)))) ^ 0.5 * (1+ 0.577/log(numseqs))
    upperThresh = (binSize*probOfOnlyCurrentMut) + 
      SDs*((binSize*probOfOnlyCurrentMut*(1 - probOfOnlyCurrentMut))^0.5)

    cutoff_freq = upperThresh
    #SDs = 6

    println("Using cutoff of $cutoff_freq occurrences")

    # account for weird observed distributions
    cutoff_freq = min(numstay, cutoff_freq)
    cutoff_freq /= numseqs
    k = 5
    centroids = Array{Any, 1}[]
    not_centroids = Array{Any, 1}[]
    freqs = proportionmap(seqs)
    for (s, f) in freqs
        if f >= cutoff_freq
            push!(centroids, [s, f])
        elseif fastmode == false
            push!(not_centroids, [s, f])
        end
    end
    if fastmode == false
        centrkmers = [kmer_count(q[1], k) for q in centroids]
        for (s, f) in not_centroids
            skc = kmer_count(s, k)
            dists = map(q -> corrected_kmer_dist(skc, q, k=k), centrkmers)
            #dists = map(q -> kmer_seeded_edit_dist(s, q[1]), centroids)
            closest = indmin(dists)
            centroids[closest][2] += f
        end
    end
    
    println("Found ", length(centroids), " centroids")
    println("Writing to ", out)
    write_fasta(out, [s[1] for s in centroids],
        names=["seq"*string(i)*"_"*string(Int(round(f[2] * numseqs))) 
            for (i, f) in enumerate(centroids)])
        
end

