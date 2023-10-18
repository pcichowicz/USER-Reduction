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

echo_time(){
    base=(${IBlue}`basename $0 ".sh"`${Color_Off})
    timestamp="${IGreen}`date "+%a %H:%M:%S"`${Color_Off}"
    echo -e "$timestamp $base $*"
}

echo_time "Running Picard MarkDuplicates" | tee ${logs_path}/${1}_mark_duplicates.log

markduplicates_counter=1

while read mapped; do
    name="$(basename $mapped)"
    I="${bam_path}/${name}_mapped.bam"
    O="${bam_path}/${name}_mapped-md.bam" 
    opt_duplicates="-OPTICAL_DUPLICATE_PIXEL_DISTANCE 12000"
    remove_duplicates="-REMOVE_DUPLICATES false" # set to true if you want to keep duplicated reads. Some tools required to have duplicates not present
    metrics_files="-METRICS_FILE ${bam_path}/${name}_mapped-md.bam.metrics"
    tag_policy="-TAGGING_POLICY All"
    validation="-VALIDATION_STRINGENCY LENIENT"

    echo_time ">Sample $mark_duplicates out of $number_samples - $name" | tee -a ${logs_path}/${1}_mark_duplicates.log
    echo_time "picard MarkDuplicates $I $O $opt_duplicates $remove_duplicates $metrics_files $tag_policy $validation" | tee -a ${logs_path}/${1}_mark_duplicates.log
    java -jar /usr/local/Cellar/picard/picard.jar MarkDuplicates -I $I -O $O $opt_duplicates $remove_duplicates $metrics_files $tag_policy $validation

    # index mark dupliated bam file
    echo_time "@samtools index -b $O" | tee -a ${logs_path}/${1}_mark_duplicates.log
    samtools index -b $O

    ((markduplicates_counter++))
done < ${project_dir}/${1}_fastqs.list

echo_time "Finished with marking and removing duplicate reads" | tee -a ${logs_path}/${1}_mark_duplicates.log