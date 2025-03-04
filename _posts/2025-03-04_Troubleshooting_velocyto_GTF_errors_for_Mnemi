---
layout: post
title: Troubleshooting velocyto GTF errors for Mnemi
date: '2025-03-04'
categories: Analysis
tags: scRNAseq, orthofinder
---

This post is to document the steps I took to troubleshoot the velocyto GTF errors for Mnemi. For some reason, this GTF file is not being read in correctly and when I modify i, I am getting unusually high unpliced percentages. 

# Conda setup 

```bash
# Create new conda environment for velocyto
source anaconda3/bin/activate 
conda create -n velocyto python=3.8 -y

# Activate environment
conda activate velocyto
# Install velocyto and dependencies
conda install -c bioconda velocyto.py -y
conda install -c conda-forge numpy scipy cython numba matplotlib scikit-learn h5py click -y
conda install -c bioconda samtools

# Install additional R dependencies
conda install -c conda-forge r-base r-seurat r-matrix -y

# Verify installation
velocyto --help
```

# Concatenated GTF file with mitochondrial genome

```bash
cat Mle_F31_T2T_BrakerAnnotationFinal.gtf 2023_mitogenome_mle.gtf > Mle_v3_cat2.gtf
```


# Original GTF file
```bash
(velocyto) -bash-4.2$ cut -f3 Mle_v3_cat2.gtf | sort | uniq -c

 205250 CDS
 205284 exon
  14360 five_prime_UTR
 183200 intron
   1521 mRNA
  22042 start_codon
  22049 stop_codon
  12910 three_prime_UTR

(velocyto) -bash-4.2$ head -n 20 Mle_v3_cat2.gtf
chromosome_1    GeneMark.hmm3   start_codon     41718   41720   .       +       0       transcript_id "anno2.3794_t"; gene_id "g_16414";
chromosome_1    GeneMark.hmm3   mRNA    41718   42227   .       +       .       transcript_id "anno2.3794_t"; gene_id "g_16414";
chromosome_1    GeneMark.hmm3   CDS     41718   42227   .       +       0       transcript_id "anno2.3794_t"; gene_id "g_16414";
chromosome_1    GeneMark.hmm3   exon    41718   42227   .       +       0       transcript_id "anno2.3794_t"; gene_id "g_16414";
chromosome_1    GeneMark.hmm3   stop_codon      42225   42227   .       +       0       transcript_id "anno2.3794_t"; gene_id "g_16414";
chromosome_1    AUGUSTUS        start_codon     149260  149262  .       +       0       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        CDS     149260  149335  0.67    +       0       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        exon    149260  149335  .       +       .       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        intron  149336  149431  1       +       .       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        CDS     149432  149524  1       +       2       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        exon    149432  149524  .       +       .       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        intron  149525  149740  1       +       .       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        CDS     149741  150193  0.87    +       2       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        exon    149741  150193  .       +       .       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        intron  150194  166802  0.15    +       .       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        CDS     166803  166807  0.15    +       2       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        exon    166803  166807  .       +       .       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    AUGUSTUS        stop_codon      166805  166807  .       +       0       transcript_id "anno1.g689.t1"; gene_id "g_666";
chromosome_1    stringtie2utr   three_prime_UTR 166706  166999  1000    -       .       transcript_id "anno1.g690.t1"; gene_id "g_667";
chromosome_1    AUGUSTUS        stop_codon      167261  167263  .       -       0       transcript_id "anno1.g690.t1"; gene_id "g_667";


(velocyto) -bash-4.2$ head -n 20 Mle_v3_cat2.gtf | grep "exon" | cut -f9 | head -n 2
transcript_id "anno2.3794_t"; gene_id "g_16414";
transcript_id "anno1.g689.t1"; gene_id "g_666";
```

```bash 
# Check feature ordering within transcripts
awk -F"\t" '
{
    if($3 == "exon" || $3 == "intron") {
        # Store feature type and coordinates
        if(match($9, /transcript_id "([^"]+)"/, t)) {
            tid = t[1]
            features[tid] = features[tid] sprintf("%s:%d-%d;", $3, $4, $5)
        }
    }
}
END {
    # Print first 5 transcripts feature structure
    count = 0
    for(tid in features) {
        if(++count <= 5) {
            print "Transcript", tid
            print features[tid]
            print "---"
        }
    }
}' Mle_v3_cat2.gtf

# Check for overlapping features
awk -F"\t" '
$3 == "exon" || $3 == "intron" {
    if(match($9, /transcript_id "([^"]+)"/, t)) {
        tid = t[1]
        # Check overlap with previous feature
        if(last_end[tid] && $4 <= last_end[tid]) {
            print "Overlap in transcript " tid ":"
            print "Previous feature ends at " last_end[tid]
            print "Current " $3 " starts at " $4
        }
        last_end[tid] = $5
    }
}' Mle_v3_cat2.gtf | head -n 15


Overlap in transcript MSTRG.15670.11:
Previous feature ends at 12372240
Current exon starts at 12372205
```

There are overlapping exons in the GTF file. For example, in transcript MSTRG.15670.11:
* One exon ends at position 12372240
* The next exon starts at position 12372205
* This means there's a 35bp overlap (12372240 - 12372205 = 35)
* This overlap could cause velocyto to:
    * Count the same reads multiple times
    * Incorrectly classify reads as unspliced
    * Misinterpret the transcript structure

Let's modify our script to fix these overlaps:

```bash
# Create velocyto-compatible GTF with strict ordering no overlap checking
# Create velocyto-compatible GTF with strict coordinate ordering
# Create velocyto-compatible GTF with strict filtering
awk -F"\t" 'BEGIN{OFS="\t"} 
{
    # Fix chromosome names
    if($1 ~ /^chromosome_/) {
        $1 = "ML_" $1
    } else if($1 ~ /^PG104_mtDNA/) {
        $1 = "ML_mtDNA"
    }
    
    # Store exon info by transcript
    if($3 == "exon") {
        if(match($9, /transcript_id "([^"]+)"; gene_id "([^"]+)";/, arr)) {
            tid = arr[1]
            gid = arr[2]
            chr[tid] = $1
            strand[tid] = $7
            gene[tid] = gid
            
            # Store unique exon coordinates
            key = tid "_" $4 "_" $5
            if(!(key in seen)) {
                seen[key] = 1
                exon_count[tid]++
                start_pos[tid,exon_count[tid]] = $4
                end_pos[tid,exon_count[tid]] = $5
                
                # Update transcript bounds
                if(!min_pos[tid] || $4 < min_pos[tid]) min_pos[tid] = $4
                if(!max_pos[tid] || $5 > max_pos[tid]) max_pos[tid] = $5
            }
        }
    }
}
END {
    # Process transcripts
    for(tid in exon_count) {
        if(exon_count[tid] > 1) {  # Skip single-exon transcripts
            # Sort coordinates
            for(i=1; i<exon_count[tid]; i++) {
                for(j=i+1; j<=exon_count[tid]; j++) {
                    if(start_pos[tid,i] > start_pos[tid,j]) {
                        # Swap starts
                        temp = start_pos[tid,i]
                        start_pos[tid,i] = start_pos[tid,j]
                        start_pos[tid,j] = temp
                        # Swap ends
                        temp = end_pos[tid,i]
                        end_pos[tid,i] = end_pos[tid,j]
                        end_pos[tid,j] = temp
                    }
                }
            }
            
            # Check for overlaps
            valid = 1
            for(i=1; i<exon_count[tid]; i++) {
                if(end_pos[tid,i] >= start_pos[tid,i+1]) {
                    valid = 0
                    break
                }
            }
            
            if(valid) {
                # Print transcript
                print chr[tid], "velocyto", "transcript", min_pos[tid], max_pos[tid], \
                      ".", strand[tid], ".", "gene_id \"" gene[tid] "\"; transcript_id \"" tid "\";"
                
                # Print first exon
                print chr[tid], "velocyto", "exon", start_pos[tid,1], end_pos[tid,1], \
                      ".", strand[tid], ".", "gene_id \"" gene[tid] "\"; transcript_id \"" tid "\";"
                
                # Print alternating introns and exons
                for(i=2; i<=exon_count[tid]; i++) {
                    # Print intron
                    print chr[tid], "velocyto", "intron", end_pos[tid,i-1] + 1, start_pos[tid,i] - 1, \
                          ".", strand[tid], ".", "gene_id \"" gene[tid] "\"; transcript_id \"" tid "\";"
                    # Print exon
                    print chr[tid], "velocyto", "exon", start_pos[tid,i], end_pos[tid,i], \
                          ".", strand[tid], ".", "gene_id \"" gene[tid] "\"; transcript_id \"" tid "\";"
                }
            }
        }
    }
}' Mle_v3_cat2.gtf > Mle_v3_cat2.velocyto.strict.gtf
```
```bash
# Check the output
echo "Feature counts in new GTF:"
cut -f3 Mle_v3_cat2.velocyto.strict.gtf | sort | uniq -c

 200020 exon
 183201 intron
  16819 transcript
```

```bash
echo -e "\nFirst few entries of new GTF:"
head -n 15 Mle_v3_cat2.velocyto.strict.gtf
```

```bash
(velocyto) -bash-4.2$ head -n 15 Mle_v3_cat2.velocyto.strict.gtf
ML_chromosome_2 velocyto        transcript      3504200 3507988 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        exon    3504200 3504228 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        intron  3504229 3504348 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        exon    3504349 3504382 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        intron  3504383 3504457 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        exon    3504458 3504568 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        intron  3504569 3504976 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        exon    3504977 3505964 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        intron  3505965 3506359 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        exon    3506360 3506547 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        intron  3506548 3506671 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        exon    3506672 3506928 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        intron  3506929 3507328 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        exon    3507329 3507467 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
ML_chromosome_2 velocyto        intron  3507468 3507608 .       -       .       gene_id "g_11506"; transcript_id "anno2.MSTRG.7124.4";
```

Lets do more more checks

```bash
# 1. Check for proper exon-intron alternation
awk -F"\t" '
$3 == "exon" || $3 == "intron" {
    if(match($9, /transcript_id "([^"]+)"/, t)) {
        tid = t[1]
        if(last_type[tid] == $3) {
            print "ERROR: Consecutive " $3 "s in transcript " tid
        }
        last_type[tid] = $3
    }
}' Mle_v3_cat2.velocyto.strict.gtf

# 2. Check for coordinate consistency
awk -F"\t" '
{
    if(match($9, /transcript_id "([^"]+)"/, t)) {
        tid = t[1]
        if(last_end[tid] && $4 <= last_end[tid]) {
            print "ERROR: Overlapping/incorrect order in " tid ":"
            print "Previous feature ends at " last_end[tid]
            print "Current " $3 " starts at " $4
        }
        last_end[tid] = $5
    }
}' Mle_v3_cat2.velocyto.strict.gtf

# 3. Check intron validity
awk -F"\t" '
$3 == "exon" || $3 == "intron" {
    if(match($9, /transcript_id "([^"]+)"/, t)) {
        tid = t[1]
        if($3 == "intron") {
            if(!prev_exon_end[tid]) {
                print "ERROR: Intron without preceding exon in " tid
            } else if($4 != prev_exon_end[tid] + 1) {
                print "ERROR: Gap between exon and intron in " tid
                print "Exon ends at " prev_exon_end[tid]
                print "Intron starts at " $4
            }
        } else { # exon
            prev_exon_end[tid] = $5
        }
    }
}' Mle_v3_cat2.velocyto.no_overlaps.gtf

# 4. Check transcript boundaries
awk -F"\t" '
{
    if(match($9, /transcript_id "([^"]+)"/, t)) {
        tid = t[1]
        if($3 == "transcript") {
            trans_start[tid] = $4
            trans_end[tid] = $5
        } else if($4 < trans_start[tid] || $5 > trans_end[tid]) {
            print "ERROR: Feature outside transcript bounds in " tid
            print "Transcript: " trans_start[tid] "-" trans_end[tid]
            print "Feature: " $4 "-" $5
        }
    }
}' Mle_v3_cat2.velocyto.no_overlaps.gtf

# 5. Check for multi-exon transcripts
awk -F"\t" '
{
    if(match($9, /transcript_id "([^"]+)"/, t)) {
        tid = t[1]
        if($3 == "exon") exon_count[tid]++
    }
}
END {
    for(tid in exon_count) {
        if(exon_count[tid] == 1) {
            print "WARNING: Single-exon transcript found: " tid
        }
    }
}' Mle_v3_cat2.velocyto.no_overlaps.gtf

# 6. Verify gene_id and transcript_id presence
awk -F"\t" '
{
    if(!match($9, /gene_id "[^"]+"/) || !match($9, /transcript_id "[^"]+"/)) {
        print "ERROR: Missing gene_id or transcript_id in line:"
        print $0
    }
}' Mle_v3_cat2.velocyto.no_overlaps.gtf
```


