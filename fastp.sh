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

# fastp options to be used
# --in1
# --in2
# --out1
# --out2
# --unpaired1
# --unpaired2
# --merge
# --merged_out
# --failed_out
# --detect_adapter_for_pe
# --dont_eval_duplication
# --length_required
# --qualified_quality_phred
# --report_title
# --thread

    # Import name of config from run_scripts.sh 
project_dir="/Users/Patrick/aDNA/Project/${1}"
source /Users/Patrick/aDNA/temp/colors.txt
source ${project_dir}/${1}.config

echo_time(){
    base=(${IBlue}`basename $0 ".sh"`${Color_Off})
    timestamp="${IGreen}`date "+%a %H:%M:%S"`${Color_Off}"
    echo -e "$timestamp $base $*"
}

echo_time "Running adapter trimming and quality control script for project $1" | tee ${logs_path}/${1}_fastp.log

fastq_counter=1

while read samples; do
    name="$(basename $samples ".fastq.gz")"
    
    echo_time ">Sample $fastq_counter out of $number_samples - $name" | tee -a ${logs_path}/${1}_fastp.log

        # variables to be used for fastp command, SE or PE
    in1="--in1 ${samples}_S1_L004_R1_001.fastq.gz"
    in2="--in2 ${samples}_S1_L004_R2_001.fastq.gz"
    out1="--out1 ${fastp_path}/${name}_trimmed_out_1.fastq.gz"
    out2="--out2 ${fastp_path}/${name}_trimmed_out_2.fastq.gz"
    unpair_1="--unpaired1 ${fastp_path}/${name}_trimmed_unpaired.fastq.gz" # unpaired 1/2 can be same file
    unpair_2="--unpaired2 ${fastp_path}/${name}_trimmed_unpaired.fastq.gz" # unpaired 1/2 can be same file
    merge="--merge"
    merge_out="--merged_out ${fastp_path}/${name}_trimmed_merged.fastq.gz"
    failed_out="--failed_out ${fastp_path}/${name}_trimmed_failed_out.fastq.gz"
    pe_adapter="--detect_adapter_for_pe"
    dont_eval="--dont_eval_duplication"
    l="--length_required 30"
    q="--qualified_quality_phred 30"
    report_title="--report_title '$name'"
    html="--html ${fastp_path}/${name}.html"
    json="--json ${fastp_path}/${name}.json"
    w="--thread"
    
        # config file stores the variable, indicating if PE or SE
    if [ "$sequencing_type" == "SE" ]; then
    
        echo_time "@$fastp $l $q $in1 $out1 $failed_out $dont_eval $report_title $html $json\n" | tee -a ${logs_path}/${1}_fastp.log
        fastp $l $q $in1 $out1 $failed_out $dont_eval $report_title $html $json
    else
        
        echo_time "@fastp $l $q $in1 $in2 $out1 $out2 $unpair_1 $unpair_2 $merge $merge_out $failed_out $pe_adapter $dont_eval $report_title $html $json\n" | tee -a ${logs_path}/${1}_fastp.log
        fastp $l $q $in1 $in2 $out1 $out2 $unpair_1 $unpair_2 $merge $merge_out $failed_out $pe_adapter $dont_eval $report_title $json $html
    fi
    
    ((fastq_counter++))
done < ${project_dir}/${1}_fastqs.list # sourced/imported from the config file that has this information

echo_time "fastp trimming and quality control finished." | tee -a ${logs_path}/${1}_fastp.log