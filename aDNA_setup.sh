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

###########################################################################################
###########################################################################################
# 
# ONLY EDIT Root paths to project, fastq files and optional "naming_convention" to set
# platform used for sequencing.
#
# ex. Illumina Platform sequencing raw fq file naming format
#     sample_name_S1_L001_R1_001.fastq.gz , script will remove the "_S1_L001_R1_001.fastq.gz"
#     tag in order to get unique sample names to populate fastq list file
#

set -e

    ############################################################################################
    #                                    Variables/paths                                       #
    ############################################################################################ 

    # Path to Scripts directory for color text output used for scripts
source /Users/Patrick/aDNA/temp/colors.txt

    # Root path to project directory where all output and inputs will be used/placed
project_path="/Users/Patrick/aDNA/Project"

    # Change this to your fastq directory where fastq files are downloaded
rawfqs_path="/Users/Patrick/aDNA/raw_fq"

    # Set this to type of platform used, must edit platform_type function if not illumina
naming_convention="illumina"

    # Paths to relevant directories for reference genomes
    # GRCh38/hg38
    ref_path_38="/Users/Patrick/aDNA/REFERENCE/GRCh38/Homo38"
    reference_38="genome38.fa"

    #GRCh37/hg19
    ref_path_19="/Users/Patrick/aDNA/REFERENCE/UCSC/Homo"
    reference_19="genome.fa"

    ############################################################################################
    #                                          Functions                                       #
    ############################################################################################ 


    # Function for displaing usage
display_usage() {
    echo -e "\nUsage $0 -n <Project Name> -s <Single-end or Paired-end reads> -r <Reference genome build>\n"
    echo -e "This script requires three arguments; Project name, sequencing type, reference genome build\n"

    echo -e "Example"
    echo -e "$0 -n SEA -seq PE -r hg18\n"  
    echo -e "-h \t| --help: \t\t\t Prints help page.\n"
    echo -e " ----- Required ----- "
    echo -e "-n \t| --name: \t\t Name of project (Project ID)."
    echo -e "-s \t| --sequencing: \t Sequencing type"
    echo -e "\t    <SE || se || single || single-end>\tsingle-end"
    echo -e "\t    <PE || pe || paired || paired-end>\tpaired-end"
    echo -e "-r \t| --reference: \t\t Reference genome build." 
    echo -e "\t    <GRCh38 || hg38>\t2013 human reference build"
    echo -e "\t    <GRCh37 || hg19>\t2009 human reference build"
}

    # Function for checking if directory exists
check_directory() {
    if [ ! -d "$1" ]; then
        echo -e "$1 directory not found. Make sure directory exists or correct path is specified."
        exit 1
    else
        echo
    fi
}

    # Function to check for genome reference build index files
check_ref_index(){
    local index_files=("amb" "ann" "bwt" "pac" "sa")
    local missing_index=()
    local valid_index=()
        #If directory does not exists, script terminates
    check_directory $1
    
    for index in "${index_files[@]}"; do
    name="${ref_path}/${reference}.${index}"

        if [[ (! -f "$name") || (! -s "$name") ]]; then
            missing_index+=($index)
        else
            valid_index+=($index)
        fi 
    done
    
        #Print out missing/empty index files to fix
    if [[ ${#missing_index[@]} -gt 0 ]]; then
        echo -e "${#missing_index[@]} index file(s) missing:${BRed} ${missing_index[*]}.${Color_Off}\n
        Please index reference build with \"bwa index <reference genome build>\" to fix issue. Aborting ..."
        exit 1
    
    fi
        # If all index are valid, continue
    if [[ ${#valid_index[@]} -eq 5 ]]; then
        echo -e "All ${#valid_index[@]} reference genome index files are valid:${BGreen} ${valid_index[*]}${Color_Off}."
    fi
}
        # Function to set what kind of sequencing platform was used
platform_type(){
    case "$1" in
        illumina) 
            naming_convention="S1_L004_R[1,2]_001.fastq.gz"
            platform="ILLUMINA"
            ;;
        nanopore | nano) 
            naming_convention="nano"
            platform="L002_R[1,2]_001.fastq"
            ;;
        *) 
            naming_convention="unknown"
            echo -e "\nSequencing platform not recognized"
            echo -e "Accepted arguments are: Illumina or Nanopore / nano\n"
            ;;
    esac
}

        # Function to print date infront of every line/command
echo_time(){
    base=(${IBlue}`basename $0 ".sh"`${Color_Off})
    timestamp="${IGreen}`date "+%a %H:%M:%S"`${Color_Off}"
    echo -e "$timestamp $base $*"
}

#
#
#
###########################################################################################
###########################################################################################

############################################################################################
#                                Setup script starts here                                  #
############################################################################################

echo_time "Starting aDNA analysis pipeline setup script\n"

    # Check if correct number of arguments are provided
if [ ! "$#" -eq 6 ]; then
    display_usage
    exit 1
fi

    # Parse script command-line arguments, allows for any order of arguments
    # Also sets the reference genome path and file name based on reference build

while [ "$#" -gt 0 ]; do
    case "$1" in
        -n | --name) project_name="$2"; 
        shift 2 ;;
        -s | --sequencing) data_type="$2"; 
            case "$data_type" in
                 "SE" | "se" | "single" | "SINGLE" | "single-end")
                 data_type="$2"
                 ;;
                 "PE" | "pe" | "paired" | "PAIRED" | "paired-end")
                 data_type="$2"
                 ;;
                 *)
                 echo -e "Sequencing type not recognized"
                 display_usage
                 exit 1
            esac
        shift 2 ;;
        -r | --reference) reference_build="$2";
            case "$reference_build" in
                 "GRCh38" | "hg38")
                 ref_path=$ref_path_38
                 reference=$reference_38
                 ;;
                 "GRCh37" | "hg19")
                 ref_path=$ref_path_19
                 reference=$reference_19
                 ;;
                 *)
                 echo -e "Build reference not recognized"
                 display_usage
                 exit 1 
            esac
        shift 2 ;;
        *) echo "Unknown option: $1"; 
        display_usage; 
        exit 1 ;;
    esac
done

    #Check if values are present and correct
if [ -z "$project_name" ] || [ -z "$data_type" ] || [ -z "$reference_build" ]; then
    echo "Missing required arguments"
    display_usage
    exit 1
fi

    # Create Project name/id/ directories for output files

if [ ! -d "$project_path/$project_name" ]; then
    echo_time "Directories do not exist for project name, creating directories now ..."
    mkdir -p $project_path/$project_name/{fastp,bwa,bam,mapDamage,logs,statistics}
    find $project_path/$project_name -type d -exec sh -c "echo '${IGreen}$(date +"%a %H:%M:%S")${Color_Off} $base ${Yellow}{}${Color_Off}'" \;
    echo
    sleep 0.5
else 
    echo_time "Directories already exist for project $project_name ..."
    find $project_path/$project_name -type d -exec sh -c "echo '${IGreen}$(date +"%a %H:%M:%S")${Color_Off} $base ${Yellow}{}${Color_Off}'" \;
    sleep 0.5
fi

    # Check/confirm genome index files are in good status
ref_index="$(check_ref_index "$ref_path" "$reference")"
echo_time $ref_index
    
    # Read in directory name of rawfq files to be processed
    # MUST BE LOCATED IN THE ROOT DIRECTORY OF $rawfqs_path

echo_time "Enter directory name of raw fq files:"
read raw_fq

    # Check if correct fastq directory was provided and create txt file of all sample names in directory
if [ ! -d "$rawfqs_path/$raw_fq" ]; then
    echo_time "$raw_fq directory not found. Make sure the directory exists or correct name was entered."
    exit 1
else
    echo_time "Creating list of fastqs to be processed ..."
    sleep 0.5
    
    if [ -z $naming_convention ]; then
        echo_time "Sequencing convention variable not set."
        naming_convention="unknown"
        sleep 0.5

        while [ "$naming_convention" == "unknown" ]; do

            echo_time "What type of platform was used?"
            read -p "" input_platform 
            platform=$(echo $input_platform | tr "[:upper:]" "[:lower:]")
            platform_type "$platform"
            sleep 0.5
        done
        
        fastq_samples=($(find ${rawfqs_path}/${raw_fq} -maxdepth 1 -type f -name "*.fastq.gz" | sed "s/_${naming_convention}//g" | sort | uniq))
        sleep 0.5

    else
        
        echo_time "Seqencing platform set to: ${Green}$naming_convention${Color_Off}"
        platform_type $naming_convention
        
        fastq_samples=($(find ${rawfqs_path}/${raw_fq} -maxdepth 1 -type f -name "*.fastq.gz" | sed "s/_${naming_convention}//g" | sort | uniq))
    fi
    
    # Number of unique fastq samples to be processed. Note each sample will have R1 and R2 if PE, but list will not indicate this.
    fastq_n=(${#fastq_samples[@]})
        
    echo_time "Number of samples; $fastq_n"
    

    if [ -e "${project_name}_fastqs.list" ]; then
        echo_time "Creating new samples list file for project: ${Green}${project_name}_fastqs.list${Color_Off}"
        > "${project_path}/${project_name}/${project_name}_fastqs.list"
    else
        echo_time "Creating fastq sample list: ${Green}${project_name}_fastqs.list${Color_Off}"
        > "${project_path}/${project_name}/${project_name}_fastqs.list"
    fi
    echo -e "${Yellow}${fastq_samples[*]}${Color_Off}" | tr " " "\n" | pr -to31 #| tee ${project_path}/${project_name}/${project_name}_fastqs.list
    echo -e "${Yellow}${fastq_samples[*]}${Color_Off}" | tr " " "\n" >> ${project_path}/${project_name}/${project_name}_fastqs.list
fi

    # Create config file where all meta-information of automation script holds and stores information for next script.
    # If this script is run several times, config metadata will not be duplicated and appended to the existing file,
    # new file will be created and metadata will be written again.

if [ -e "${project_path}/${project_name}.config" ]; then
    echo_time "Creating new config file: ${Green}${project_name}.config${Color_Off}\n"
    > "${project_path}/${project_name}/${project_name}.config"

else
    echo_time "Creating config file: ${Green}${project_name}.config${Color_Off}"
    > "${project_path}/${project_name}/${project_name}.config"
fi

config_file="${project_path}/${project_name}/${project_name}.config"

    # Indexed array with configuration variables
config_variables=(
    "project_name=\"$project_name\""
    "project_dir=\"$project_path/$project_name\""
    "reference_build=\"$reference_build\""
    "sequencing_type=\"$data_type\""
    "sequencing_platform=\"$platform\""
    "rawfqs_path=\"$rawfqs_path/$raw_fq\""
    "number_samples=\"$fastq_n\""
    "reference_genome=\"${ref_path}/${reference}\""
    "fastp_path=\"${project_path}/${project_name}/fastp\""
    "bwa_path=\"${project_path}/${project_name}/bwa\""
    "bam_path=\"${project_path}/${project_name}/bam\""
    "mapDamage_path=\"${project_path}/${project_name}/mapDamage\""
    "logs_path=\"${project_path}/${project_name}/logs\""
    "statistics_path=\"${project_path}/${project_name}/statistics\""
)

    # Loop through the indexed array and write each variable to the config file
for var in "${config_variables[@]}"; do
    echo -e "${Yellow}${var}" | pr -to31
    echo -e "${Yellow}${var}" >> $config_file
done

echo_time "Setup script complete, run ${IGreen}./aDNA_run_scripts${Color_Off} with options for analysis."

    
    
    # Working progress to add in more relevant config settings/metadata

# ID - read group identifier
# BC - sample/library barcode
# CN - Sequencing center
# DT - Date run was produced
# LB - Library
# PL - Platform/technology used (Must be all Capital?)
# DS - Description may be used
# SM - Sample, use pool name if pool is being sequenced


