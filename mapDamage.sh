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
source /Users/Patrick/aDNA/temp/colors.txt
source ${project_dir}/${1}.config

echo_time(){
    base=(${IBlue}`basename $0 ".sh"`${Color_Off})
    timestamp="${IGreen}`date "+%a %H:%M:%S"`${Color_Off}"
    echo -e "$timestamp $base $*"
}

echo_time "Running mapDamage script ..."

mapDam_counter=1
while read samples; do
    name=$(basename $samples)
    bam_input=${bam_path}/${name}_mapped-md.bam
    folder=${mapDamage_path}/${name}_mapdamage
    
    echo_time ">Sample $mapDam_counter out of $number_samples; samse - $name" | tee -a ${logs_path}/${1}_mapDamage.log
    echo_time "@mapDamage -Q 30 -d $folder -i $bam_input  -r ${reference_genome}" | tee -a ${logs_path}/${1}_mapDamage.log
    mapDamage -Q 30 -d $folder -i $bam_input  -r ${reference_genome}

((mapDam_counter++))
done < ${project_dir}/${1}_fastqs.list