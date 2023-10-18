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
source /Users/Patrick/aDNA/temp/colors.txt

echo_time(){
    base=(${IBlue}`basename $0 ".sh"`${Color_Off})
    timestamp="${IGreen}`date "+%a %H:%M:%S"`${Color_Off}"
    echo -e "$timestamp $base $*"
}

# Check if at least one parameter is provided
if [ $# -eq 0 ]; then
  echo "Please provide at least one script name as a parameter."
  echo "fastp | bwa_aln | mark_duplicates | mapdamage"
  exit 1
fi

echo_time "Enter project name for config file:"
read config_file

# # Initialize flags for each script
fastp_flag=false
bwa_aln_flag=false
mark_duplicates_flag=false
mapdamage_flag=false

  # Process the parameters, set to true
for param in "$@"; do
  case $param in
    fastp)
      fastp_flag=true
      ;;
    bwa_aln)
      bwa_aln_flag=true
      ;;
    mark_duplicates)
      mark_duplicates_flag=true
      ;;
    mapdamage)
      mapdamage_flag=true
      ;;
    *)
      echo "Invalid script name: $param"
      ;;
  esac
done

# Execute the scripts in the specified order
if [ "$fastp_flag" = true ]; then
  echo_time "Running fastp.sh script..."
  ./fastp.sh $config_file
fi

if [ "$bwa_aln_flag" = true ]; then
  echo_time "Running bwa_aln.sh script..."
  ./bwa_aln.sh $config_file
fi

if [ "$mark_duplicates_flag" = true ]; then
  echo_time "Running remove_duplicates.sh script..."
  ./mark_duplicates.sh $config_file
fi

if [ "$mapdamage_flag" = true ]; then
  echo_time "Running mapdamage.sh script..."
  ./mapdamage.sh $config_file
fi

echo_time "Running of all script(s) complete, log files are located in project name directory."
