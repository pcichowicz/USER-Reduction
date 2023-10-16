# USER-Reduction
MSc Thesis project analysis pipeline of aDNA and USER Treatment . Supervised by Thorfinn Korneliussen (GLOBE Institute, Section for Geogenetics - Assistant Professor) and Hugh McColl (GLOBE Institute - Section for Geogenetics - Postdoc at Willerslev Group)

Abstract :
High through sequencing data has revolutionized many different research fields in DNA studies and wetlab and drylab protocols has been developed for modern DNA. The vast amounts of sequencing data that can be generated with recent sequencing machines revolutionized entire research fields and has pushed the envelope for not just answering existing questions but formed the basis of a paradigm shift in fields previously unrelated to genetics, such as archeology and linguistics. However due to the specific problematic nature of ancient DNA these methods might not be optimal for degraded DNA material, and require additional preparations in laboratories or new statistical models to better handle and prepare ancient DNA samples for analysis. In this thesis we will investigate  the performance of recently developed wetlab protocols that specifically aim to solve some of the idiosyncrasies of aDNA data. Specifically these are UDG treatment (USER enzyme treatment) which is responsible for removing damage signals in ancient DNA and "cleaning" the DNA in order for better generation of aDNA data. Library preparation is both costly and time consuming to prepare samples for sequencing. The work in this thesis aims to evaluate if the enzymatic removal of uracil residues from ancient DNA extracts prior to the generation of NGS libraries shows the same efficency and data generation by using 1/4 of the original volume of enzyme. The reduction of the USER enzyme would result in a reduction of cost per 96 library preparations up to 75\%, from 9,722 DKK to 2,443 DKK per plate of reactions.

# Bioinformatic Pipeline

**(WORK IN PROGRESS TO UPLOAD - Currently improving scripts for more robustness)**

Here are the steps that were taken to analyze NGS data from aDNA samples that have been treated with USER enzyme to remove DNA damage and clean the DNA for analysis. The pipeline is as follows;

- Fastp - Adapter and trimming quality control of PE data;
  quality and lenght >= 30

- BWA aln - Alignment and mapping
  bwa aln is used to favor the short length of aDNA

- SAMtools - Utility for human readable output. SAM/BAM files
  Sort, merge, extracting mapped reads only, 

- Picard - Marking and removing duplicate reads

- mapDamage - Nucleotide misincorporations patterns

- ContaMix - Mitochondrial contamination

- ANGSD - X chromosome contaminaion

- Summary statistics
  - Total reads - Number of total reads from sequencing;
  - Total trimmed reads - Number of reads after fastp trimming, removing adapter sequences, length ≥ 30, quality score ≥ 30;
  - Read length - Mean read length of mapped reads;
  - Total # of mapped reads - Number of mapped reads to human reference genome build GRCh37 (hg19);
  - Number of duplicate reads - Number of reads that are PCR & cluster duplicates;
  - Number of unique reads - Number of reads that uniquely map to the human reference genome;
  - Mapping - Proportion of reads mapping;
  - Clonality - Proportion of reads that are duplicates ( 1 - unique mapped reads / number of mapping reads);
  - Endogeous - Proportion of human reads including duplicates, after trimming adapters and short fragments. (mapping Reads / Number of reads trimmed);
  - Unique mapping - Proportion of human reads excluding duplicates, after trimming adapters and short fragments (unique mapping reads / number of reads trimmed);
  - Effeciency - Proportion of human reads excluding duplicates, out of total reads sequenced (unique reads / total number of reads)l;
  
## Set up
This pipeline has specific input/output directory structure. Below shows the instructions of how to setup the directories for this pipeline.

```
|-- Project
|-- Reference
|-- raw_fqs
|-- Scripts
```

### Required inputs
The aDNA pipeline requires the users to have the following to be located in the right directories:

* Raw fastq files to be downloaded and stored in the `raw_fqs` directory with project name/ID.
* Indexed reference genome in the `Reference` directory that will be used. Can have multiple genome builds/species, however other than hg19 will have to be manually added in the 'setup_script.sh'
* All script files to be placed in the `Scripts` directory
* The `Project` directory should be empty (can have project name/ID directories if pipeline has been used before)

#### Running the `setup_script.sh`
Script requires 3 input arguments, -n <Project name/ID>, -s <Single/Paired-end reads>, and -r <reference build>;
'eg. ./setup_script.sh -n SEA -s PE -r hg19'

This script will require two user inputs, creating a 'project name/ID.config' (stores all the metadata needed) file that will be used by the following `run_scripts.sh` where depending on which tools selected will use the 'project name/ID.config' metadata. Input prompts are as follows:

* Directory name of fastq files (relative name), 'eg . SEA'
* Type of sequencing platform used (can be set in script if all reads are from the same platform). 'eg. Illumina'

Once the script is finished, the '.config' file is produced and contains the following information;

* project_name="SEA"
* project_dir="/Users/Patrick/aDNA/Project/SEA"
* reference_build="hg19"
* sequencing_type="PE"
* sequencing_platform="illumina"
* rawfqs_path="/path_to_main_directory/raw_fq/SEA"
* number_samples="<number of samples to be processed"
* reference_genome="/path_to_main_directory/REFERENCE/UCSC/Homo/genome.fa"
* fastp_path="/path_to_main_directory/Project/single/fastp"
* bwa_path="/path_to_main_directory/Project/single/bwa"
* bam_path="/path_to_main_directory/Project/single/bam"
* final_path="/path_to_main_directory/Project/single/final"
* mapDamage_path="/path_to_main_directory/Project/single/mapDamage"
* logs_path="/path_to_main_directory/Project/single/logs"
* statistics_path="/path_to_main_directory/Project/single/statistics"
