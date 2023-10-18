#!/bin/bash
###########################################################################################
###########################################################################################
#
# Author : Patrick Cichowicz
# 
# Masters Thesis project script for aDNA and UDG (USER) enzyme reduction treatment analysis
# University of Copenhagen, GLOBE Institute, SUND
# 
###########################################################################################
###########################################################################################

    # Script to generate statistics from bioinformatic pipeline
    
    # number of total read pairs sequenced,
    # number of reads trimmed with fastp
    # number of reads mapped to reference genome
    # number of cluster duplicates
    # number of PCR duplicates
    # number of unique mapped reads
    # proportion of reads after trimming (total reads pairs / reads after trimming)
    # proportion reads mapping ( number of mapped reads / total)
    # Clonality -  proportion of duplicate reads
    # Cluster Duplicates - proportion of reads that are cluster duplicates
    # PCR duplicates - proportion of reads that are PCR duplicates
    # Endogenous - proportion of humans reads (including duplicates) after trimming reads
    # Mapped Unigue - proportion of human reads (excluding duplicates) after trimming reads
    # Efficiency - Proportion of human reads (excluding duplicates) out of total number of reads
    # Read Length - Mean read length of mapped reads, excluding duplicate reads
    # xDepth=${depth[5]}
    # yDepth=${depth[8]}
    # mDepth=${depth[11]}
    # autDepth=${depth[14]}

    project_dir="/Users/Patrick/aDNA/Project/${1}"
    source /Users/Patrick/aDNA/temp/colors.txt
    source ${project_dir}/${1}.config

while read samples; do
    name=$(basename $samples)
        
        # Header
    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
    "Sample ID" \
    "Sample Age" \
    "Library" \
    "Treatment" \
    "Total Reads" \
    "Trimmed Reads" \
    "Mapped Reads" \
    "Cluster Duplicates" \
    "PCR Duplicates" \
    "Unique Mapped" \
    "pTrimmed Reads" \
    "pMapped Reads" \
    "Clonality" \
    "pCluster Duplicates" \
    "pPCR duplicates" \
    "Endogenous" \
    "pUnique Mapped" \
    "Efficiency" \
    "Read Length" \
    "xDepth" \
    "yDepth" \
    "mDepth" \
    "autDepth" \
    > stat.txt

    sample_id="$(echo $name | cut -f1 -d_)"
    sample_age="$(echo $name | cut -f2 -d_)"
    library="$(echo $name | cut -f4 -d_)"
    treatment="$(echo $name | cut -f6 -d_)"
    
        # total genome length
    genome_length=$(samtools view -H ${bam_path}/${bam}/${name}_mapped-md.bam | awk -vFS=: '/^@SQ/ {sum+=$3} END {print sum}')

    depth=($(samtools depth -q30 -Q20 -a ${bam_path}/${name}_mapped-md.bam | head -1000 | awk -f sexDetermination.awk))

    total_reads=$(jq '.summary.before_filtering.total_reads' fastp.json)
    trimmed_reads=$(samtools view -c ${bam_path}/${name}.bam)
    mapped_reads=$(samtools view -c ${bam_path}/${name}_mapped.bam)
    cluster_duplicates=5324
    pcr_duplicates=90000
    unique_mapping=$(samtools view -c -F 1024 ${bam_path}/${name}_mapped-md.bam) # after removing duplicate reads ?? or the _mapped.bam
    pTrimmed_reads=$(echo "scale=5; $trimmed_reads / $total_reads" | bc)
    pMapping=$(echo "scale=5; $mapped_reads/$total_reads" | bc )
    clonality=$(echo "scale=5; 1 - $unique_mapping/$mapped_reads" | bc )
    pCluster_duplicates=$(echo "scale=5; $cluster_duplicates / $mapped_reads" | bc ) # cluster_duplicates / mapped_reads
    pPCR_duplicates=$(echo "scale=5; $pcr_duplicates / $mapped_reads" | bc)
    endodgenous=$(echo "scale=5; $mapped_reads / $trimmed_reads" | bc )
    pMapped_unique=$(echo "scale=5; $unique_mapping / $trimmed_reads" | bc )
    efficiency=$(echo "scale=5; $unique_mapping / $total_reads" | bc )
    read_length=$(samtools view -F 1024 ${bam_path}/${name}_mapped-md.bam  | awk '{sum+=length($10)} END {print sum/NR}')
    xDepth=${depth[5]}
    yDepth=${depth[8]}
    mDepth=${depth[11]}
    autDepth=${depth[14]}

    printf "%s,%s,%s,%s,%i,%i,%i,%i,%i,%i,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.2f,%.5f,%.5f,%.5f,%.5f\n" \
    $sample_id \
    $sample_age \
    $library \
    $treatment \
    $total_reads \
    $trimmed_reads \
    $mapped_reads \
    $cluster_duplicates \
    $pcr_duplicates \
    $unique_mapping \
    $pTrimmed_reads \
    $pMapping \
    $clonality \
    $pCluster_duplicates \
    $pPCR_duplicates \
    $endodgenous \
    $pMapped_unique \
    $efficiency \
    $read_length \
    $xDepth \
    $yDepth \
    $mDepth \
    $autDepth \
    >> stat.txt

done < ${project_dir}/${1}_fastqs.list

# samtools view -H 012345P_ia_LV2002787650_LV3003058645_mkri16_U_mapped-md.bam | awk -vFS=: '/^@SQ/ {sum+=$3} END {print sum}'

        # SAM/BAM alignment section; mandatory fields
        ##############################################

        # 1     | 2    | 3     | 4   | 5    | 6     | 7     | 8     | 9    | 10  | 11
        # QNAME | FLAG | RNAME | POS | MAPQ | CIGAR | RNEXT | PNEXT | TLEN | SEQ | QUAL