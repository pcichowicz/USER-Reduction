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
project_dir="/Users/Patrick/aDNA/Project/${1}"
source ${project_dir}/${1}.config
source /Users/Patrick/aDNA/temp/colors.txt

# bwa_path="/Users/Patrick/aDNA/Project/${1}"
# source ${bwa_path}/${1}.config
# source /Users/Patrick/aDNA/temp/colors.txt

echo_time(){
    base=(${IBlue}`basename $0 ".sh"`${Color_Off})
    timestamp="${IGreen}`date "+%a %H:%M:%S"`${Color_Off}"
    echo -e "$timestamp $base $*"
}

echo_time "Running alignment script - BWA aln & samtools\n" | tee ${logs_path}/${1}_bwa-aln.log ${logs_path}/${1}_bam.log

echo
aln_counter=1

while read trimmed; do
    name="$(basename $trimmed ".fastq.gz")"
    #input_ref="${bwa_path}/fastp"
    #output_bwa="${bwa_path}/bwa"
    #output_bam="${bwa_path}/bam"    
    merged="${fastp_path}/${name}_trimmed_merged.fastq.gz"
    R1="${fastp_path}/${name}_trimmed_out_1.fastq.gz"
    R2="${fastp_path}/${name}_trimmed_out_2.fastq.gz"
    S1="${bwa_path}/${name}_trimmed_out_1.fastq.gz.sai"
    S2="${bwa_path}/${name}_trimmed_out_2.fastq.gz.sai"
    M1="${bwa_path}/${name}_trimmed_merged.fastq.gz.sai"

        # @RG tags for bam header, editing is based on how sample is named. Change depending on your naming convention
    ID="${name}_merged"
    ID2="${name}_paired"
    SM="$(echo $name | cut -c1-7)"
    CN="CGG" # Can hardcode string if all sequencing is done from same place etc.
    PL="$(echo $sequencing_platform | tr "[:lower:]" "[:upper:]")"
    LB="$(echo $name | cut -f4 -d_)"
    DS="aDNA USER Reduction project" # Description of sample/project, etc

        # merged reads alignment
    echo_time ">Sample $aln_counter out of $number_samples; merged reads - $name" | tee -a ${logs_path}/${1}_bwa-aln.log
    echo_time "@bwa aln ${reference_genome} $merged > ${M1}" | tee -a ${logs_path}/${1}_bwa-aln.log
    # bwa aln ${reference_genome} $merged > ${M1}

        # paired-end reads alignment (out 1 & 2)
    echo_time ">Sample $aln_counter out of $number_samples; paired-end reads 1 & 2 - $name" | tee -a ${logs_path}/${1}_bwa-aln.log
    echo_time "@bwa aln ${reference_genome} $R1 > ${S1}" | tee -a ${logs_path}/${1}_bwa-aln.log
    # bwa aln ${reference_genome} $R1 > ${S1}
    echo_time "@bwa aln ${reference_genome} $R2 > ${S2}" | tee -a ${logs_path}/${1}_bwa-aln.log
    # bwa aln ${reference_genome} $R2 > ${S2}

        # bwa samse command
    echo_time ">Sample $aln_counter out of $number_samples; samse - $name" | tee -a ${logs_path}/${1}_bwa-aln.log
    echo_time "@bwa samse -r '@RG\tID:${ID}\tSM:${SM}\tCN:${CN}\tPL:${PL}\tLB:${LB}\tDS:${DS}' ${reference_genome} $M1 $merged | samtools sort -o ${bam_path}/${name}_merged.bam" | tee -a ${logs_path}/${1}_bwa-aln.log
    # bwa samse -r "@RG\tID:${ID}\tSM:${SM}\tCN:${CN}\tPL:${PL}\tLB:${LB}\tDS:seqcenter@sund.ku.dk" ${reference_genome} $M1 $merged | samtools sort -o ${bam_path}/${name}_merged.bam 
    
        # bwa sampe command
    echo_time ">Sample $aln_counter out of $number_samples; sampe - $name" | tee -a ${logs_path}/${1}_bwa-aln.log
    echo_time "@bwa sampe -r '@RG\tID:${ID2}\tSM:${SM}\tCN:${CN}\tPL:${PL}\tLB:${LB}\tDS:${DS}' ${reference_genome} $S1 $S2 $R1 $R2 | samtools sort -o ${bam_path}/${name}_paired-reads.bam" | tee -a ${logs_path}/${1}_bwa-aln.log
    # bwa sampe -r "@RG\tID:${ID2}\tSM:${SM}\tCN:CGG\tPL:ILLUMINA\tLB:${LB}\tDS:seqcenter@sund.ku.dk" ${reference_genome} $S1 $S2 $R1 $R2 | samtools sort -o ${bam_path}/${name}_paired-reads.bam

        # Merge both bam files into one containing merged and paired end reads, extract mapped reads
    echo_time ">Sample $aln_counter out of $number_samples; samtools merge - $name" | tee -a ${logs_path}/${1}_bam.log
    echo_time "@samtools merge ${bam_path}/${name}.bam ${bam_path}/${name}*.bam" | tee -a ${logs_path}/${1}_bam.log
    # samtools merge ${bam_path}/${name}.bam ${bam_path}/${name}*.bam

        # Extract mapped reads
    echo_time "@samtools view -b -F 4 -o ${bam_path}/${name}_mapped.bam ${bam_path}/${name}.bam" | tee -a ${logs_path}/${1}_bam.log
    # samtools view -b -F 4 -o ${bam_path}/${name}_mapped.bam ${bam_path}/${name}.bam
    
        # index final bam file
    echo_time "@samtools index -b ${bam_path}/${name}_mapped.bam" | tee -a ${logs_path}/${1}_bam.log
    # samtools index -b ${bam_path}/${name}_mapped.bam
    
    ((aln_counter++))
done < ${project_dir}/${1}_fastqs.list

echo_time "Finished read alignments and creating bam files" | tee -a ${1}_bwa-aln.log ${logs_path}/${1}_bam.log