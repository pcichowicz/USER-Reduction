# USER-Reduction
MSc Thesis project analysis pipeline of aDNA and USER Treatment . Supervised by Thorfinn Korneliussen (GLOBE Institute, Section for Geogenetics - Assistant Professor) and Hugh McColl (GLOBE Institute - Section for Geogenetics - Postdoc at Willerslev Group)

Title - Comparative analysis of USER enzyme reduction for aDNA next generation sequencing

Abstract - 
High through sequencing data has revolutionized many different research fields in DNA studies and wetlab and drylab protocols has been developed for modern DNA. The vast amounts of sequencing data that can be generated with recent sequencing machines revolutionized entire research fields and has pushed the envelope for not just answering existing questions but formed the basis of a paradigm shift in fields previously unrelated to genetics, such as archeology and linguistics. However due to the specific problematic nature of ancient DNA these methods might not be optimal for degraded DNA material, and require additional preparations in laboratories or new statistical models to better handle and prepare ancient DNA samples for analysis. In this thesis we will investigate  the performance of recently developed wetlab protocols that specifically aim to solve some of the idiosyncrasies of aDNA data. Specifically these are UDG treatment (USER enzyme treatment) which is responsible for removing damage signals in ancient DNA and "cleaning" the DNA in order for better generation of aDNA data. Library preparation is both costly and time consuming to prepare samples for sequencing. The work in this thesis aims to evaluate if the enzymatic removal of uracil residues from ancient DNA extracts prior to the generation of NGS libraries shows the same efficency and data generation by using 1/4 of the original volume of enzyme. The reduction of the USER enzyme would result in a reduction of cost per 96 library preparations up to 75\%, from 9,722 DKK to 2,443 DKK per plate of reactions.

# Bioinformatic Pipeline

**(Currently improving scripts for more robustness, mainly plotting and statistics files)**

Here are the steps that were taken to analyze NGS data from aDNA samples that have been treated with USER enzyme to remove DNA damage and clean the DNA for analysis. The pipeline is as follows;

- Fastp - Adapter and trimming quality control of PE data;
  quality and lenght >= 30
- BWA aln - Alignment and mapping
  bwa aln is used to favor the short length of aDNA (higher precision relative to BWA-mem)
- SAMtools - Utility for human readable output. SAM/BAM files
  Sort, merge, extracting mapped reads only, index ...
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
  - Effeciency - Proportion of human reads excluding duplicates, out of total reads sequenced (unique reads / total number of reads);
  - xDepth - Depth of coverage of x Chromosome
  - yDepth - Depth of coverage of y Chromosome
  - mDepth  - Depth of coverage of mitochonrial Chromosome
  - authDepth - Depth of automsome Chromosomes
  - Sex - male or female based (SKOGLUND x / y)
  - HaploType - Assignment of mitochondrial haplotype
  - HaploProb - probability of the given haplotype
  - Contamix 1X estimation (upper/lower bound, mapAutentic values)
  - Contamix 5X estimation (upper/lower bound, mapAutentic values)
  - ANGSD contamination values (Method 1 and 2)


  
## Set up
This pipeline has specific input/output directory structure. Below shows the instructions of how to setup the directories for this pipeline. The base directory is the directory that houses the following;

```ShellSession
|-- Project/
|-- Reference/
   |-- hg19/
      |-- genome.fa
      |-- genome.fa.amb
      |-- genome.fa.ann
      |-- genome.fa.bwt
      |-- genome.fa.oac
      |-- genome.fa.sa
|-- raw_fqs/
   |-- SEA/
      |-- sample1.fastq.gz
      |-- sample2.fastq.gz
|-- Scripts/
   |-- aDNA_setup.sh
   |-- aDNA_run.sh
   |-- fastp.sh
   |-- bwa_aln.sh
   |-- mark_duplicates.sh
   |-- mapDamage.sh
   |-- sexDetermination.awk
   |-- generate_statistics.sh
   |-- colors.txt
```

### Required inputs
The aDNA pipeline requires the users to have the following to be located in the right directories:

* Raw fastq files to be downloaded and stored in the `raw_fqs` directory with project name/ID.
* Indexed reference genome in the `Reference` directory that will be used. Can have multiple genome builds/species, however other than hg19 will have to be manually added in the **setup_script.sh**
* All script files to be placed in the `Scripts` directory
* The `Project` directory should be empty (can have project name/ID directories if pipeline has been used before)

#### Running the `setup_script.sh`
Script requires 3 input arguments, -n <Project name/ID>, -s <Single/Paired-end reads>, and -r <reference> ;

**eg. ./setup_script.sh -n SEA -s PE -r hg19**

This script will require two user inputs, creating a directory name provided by argument, **project name/ID.config** (stores all the metadata needed) and **project_name/ID_fastq.list (contains all unique fastq reads) files that will be used by the following `run_scripts.sh`. Input prompts are as follows:

* Directory name of fastq files (relative name), **eg. SEA**
* Type of sequencing platform used (can be set in script if all reads are from the same platform). **eg. Illumina**

Once the script is finished, the **.config** file is produced and contains the following information;

```ShellSession
project_name="SEA"
project_dir="/path_to_main_directory/Project/SEA"
reference_build="hg19"
sequencing_type="PE"
sequencing_platform="illumina"
rawfqs_path="/path_to_main_directory/raw_fq/SEA"
number_samples="<number of samples to be processed"
reference_genome="/path_to_main_directory/REFERENCE/UCSC/Homo/genome.fa"
fastp_path="/path_to_main_directory/Project/single/fastp"
bwa_path="/path_to_main_directory/Project/single/bwa"
bam_path="/path_to_main_directory/Project/single/bam"
mapDamage_path="/path_to_main_directory/Project/single/mapDamage"
logs_path="/path_to_main_directory/Project/single/logs"
statistics_path="/path_to_main_directory/Project/single/statistics"
```

and the **project name/ID_fastqs.list** file containing all samples (absolute path) to be processed, one sample per line;

```ShellSession
/path_to_main_directory/raw_fq/Project_name/012345P_ia_LV2002787650_LV3003058645_mkri16_U
/path_to_main_directory/raw_fq/Project_name/012345P_ia_LV2002787650_LV3003058650_mkri16_U
/path_to_main_directory/raw_fq/Project_name/012345P_ia_LV2002787650_LV3003058655_mkri16_U
/path_to_main_directory/raw_fq/Project_name/067890T_cwc_LV2002787650_LV3003058675_mkri16_U
/path_to_main_directory/raw_fq/Project_name/067890T_cwc_LV2002787650_LV3003058680_mkri16_U
/path_to_main_directory/raw_fq/Project_name/067890T_cwc_LV2002787650_LV3003058685_mkri16_U
```

Directories now should be organized as follows;

```ShellSession
|-- Project/
   |-- SEA/
      |-- SEA.config
      |-- SEA_fastqs.list
      |-- fastp/
      |-- bwa/
      |-- bam/
      |-- mapDamage/
      |-- logs/
      |-- statistics/
|-- Reference/
   |-- hg19/
      |-- genome.fa
      |-- genome.fa.amb
      |-- genome.fa.ann
      |-- genome.fa.bwt
      |-- genome.fa.oac
      |-- genome.fa.sa
|-- raw_fqs/
   |-- SEA/
      |-- sample1.fastq.gz
      |-- sample2.fastq.gz
|-- Scripts/
   |-- aDNA_setup.sh
   |-- aDNA_run_scripts.sh
   |-- fastp.sh
   |-- bwa_aln.sh
   |-- mark_duplicates.sh
   |-- mapDamage.sh
   |-- sexDetermination.awk
   |-- generate_statistics.sh
   |-- colors.txt
```

#### Running the `run_scripts.sh`
This script requires the arguments of tools that the USER would like to perform, such as fastp, bwa_aln, mark_duplicates, etc. Order of arguments does not matter as script is designed to follow logical order of tools. If only mark duplicates is selected, the prior files and steps are required to be already done, fastp, alingment, etc.

eg.

./run_scripts.sh fastp mark_duplicates bwa_aln mapdamage

Script will then only require one user input, the project name (this is the name given to the -n <Project name> from the **setup_script.sh**). Both **.config** and **fastq.list** files will be sourced and use metadata (input/parameters/samples) for tools selected. Example of fastp, bwa aln and picard tool being used below.

```ShellSession
Wed 13:05:52 run_scripts Enter project name for config file:
single
Wed 13:05:57 run_scripts Running fastp.sh script...
Wed 13:05:57 fastp Running adapter trimming and quality control script for project single
Wed 13:05:57 fastp >Sample 1 out of 1 - 012345P_ia_LV2002787650_LV3003058645_mkri16_U
Wed 13:05:57 fastp fastp trimming and quality control finished.
Wed 13:05:57 fastp @fastp --length_required 30 --qualified_quality_phred 30 --in1 /path_to_main_directory/raw_fq/single/012345P_ia_LV2002787650_LV3003058645_mkri16_U_S1_L004_R1_001.fastq.gz --in2 /path_to_main_directory/raw_fq/single/012345P_ia_LV2002787650_LV3003058645_mkri16_U_S1_L004_R2_001.fastq.gz --out1 /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_out_1.fastq.gz --out2 /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_out_2.fastq.gz --unpaired1 /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_unpaired.fastq.gz --unpaired2 /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_unpaired.fastq.gz --merge --merged_out /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_merged.fastq.gz --failed_out /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_failed_out.fastq.gz --detect_adapter_for_pe --dont_eval_duplication --report_title '012345P_ia_LV2002787650_LV3003058645_mkri16_U' --html /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U.html --json /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U.json
Wed 13:31:52 run_scripts Running bwa_aln.sh script...
Wed 13:31:52 bwa_aln Running alignment script - BWA aln & samtools
Wed 13:31:52 bwa_aln >Sample 1 out of 1; merged reads - 012345P_ia_LV2002787650_LV3003058645_mkri16_U
Wed 13:31:52 bwa_aln @bwa aln /path_to_main_directory/REFERENCE/UCSC/Homo/genome.fa /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_merged.fastq.gz > /path_to_main_directory/Project/single/bwa/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_merged.fastq.gz.sai
Wed 13:31:52 bwa_aln >Sample 1 out of 1; paired-end reads 1 & 2 - 012345P_ia_LV2002787650_LV3003058645_mkri16_U
Wed 13:31:52 bwa_aln @bwa aln /path_to_main_directory/REFERENCE/UCSC/Homo/genome.fa /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_out_1.fastq.gz > /path_to_main_directory/Project/single/bwa/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_out_1.fastq.gz.sai
Wed 13:31:52 bwa_aln @bwa aln /path_to_main_directory/REFERENCE/UCSC/Homo/genome.fa /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_out_2.fastq.gz > /path_to_main_directory/Project/single/bwa/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_out_2.fastq.gz.sai
Wed 13:31:52 bwa_aln >Sample 1 out of 1; samse - 012345P_ia_LV2002787650_LV3003058645_mkri16_U
Wed 13:31:52 bwa_aln @bwa samse -r '@RG	ID:012345P_ia_LV2002787650_LV3003058645_mkri16_U_merged	SM:012345P	CN:CGG	PL:ILLUMINA	LB:LV3003058645	DS:aDNA USER Reduction project' /path_to_main_directory/REFERENCE/UCSC/Homo/genome.fa /path_to_main_directory/Project/single/bwa/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_merged.fastq.gz.sai /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_merged.fastq.gz | samtools sort -o /path_to_main_directory/Project/single/bam/012345P_ia_LV2002787650_LV3003058645_mkri16_U_merged.bam
Wed 13:31:52 bwa_aln >Sample 1 out of 1; sampe - 012345P_ia_LV2002787650_LV3003058645_mkri16_U
Wed 13:31:52 bwa_aln @bwa sampe -r '@RG	ID:012345P_ia_LV2002787650_LV3003058645_mkri16_U_paired	SM:012345P	CN:CGG	PL:ILLUMINA	LB:LV3003058645	DS:aDNA USER Reduction project' /path_to_main_directory/REFERENCE/UCSC/Homo/genome.fa /path_to_main_directory/Project/single/bwa/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_out_1.fastq.gz.sai /path_to_main_directory/Project/single/bwa/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_out_2.fastq.gz.sai /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_out_1.fastq.gz /path_to_main_directory/Project/single/fastp/012345P_ia_LV2002787650_LV3003058645_mkri16_U_trimmed_out_2.fastq.gz | samtools sort -o /path_to_main_directory/Project/single/bam/012345P_ia_LV2002787650_LV3003058645_mkri16_U_paired-reads.bam
Wed 13:31:52 bwa_aln >Sample 1 out of 1; samtools merge - 012345P_ia_LV2002787650_LV3003058645_mkri16_U
Wed 13:31:52 bwa_aln @samtools merge /path_to_main_directory/Project/single/bam/012345P_ia_LV2002787650_LV3003058645_mkri16_U.bam /path_to_main_directory/Project/single/bam/012345P_ia_LV2002787650_LV3003058645_mkri16_U*.bam
Wed 13:31:52 bwa_aln @samtools view -b -F 4 -o /path_to_main_directory/Project/single/bam/012345P_ia_LV2002787650_LV3003058645_mkri16_U_mapped.bam /path_to_main_directory/Project/single/bam/012345P_ia_LV2002787650_LV3003058645_mkri16_U.bam
Wed 13:31:52 bwa_aln @samtools index -b /path_to_main_directory/Project/single/bam/012345P_ia_LV2002787650_LV3003058645_mkri16_U_mapped.bam
Wed 13:31:52 bwa_aln Finished read alignments and creating bam files
Wed 13:31:52 run_scripts Running remove_duplicates.sh script...
Wed 13:31:52 mark_duplicates Running Picard MarkDuplicates
Wed 13:31:54 mark_duplicates picard MarkDuplicates /path_to_main_directory/Project/single/bam/012345P_ia_LV2002787650_LV3003058645_mkri16_U_mapped.bam /path_to_main_directory/Project/single/bam/012345P_ia_LV2002787650_LV3003058645_mkri16_U_mapped-md.bam -OPTICAL_DUPLICATE_PIXEL_DISTANCE 12000 -REMOVE_DUPLICATES false -METRICS_FILE /path_to_main_directory/Project/single/bam/012345P_ia_LV2002787650_LV3003058645_mkri16_U_mapped-md.bam.metrics -TAGGING_POLICY All -VALIDATION_STRINGENCY LENIENT
Wed 13:31:54 mark_duplicates Finished with marking and removing duplicate reads
```

Upon finishing, the Project directory should have these files;

```ShellSession
|-- Project/
   |-- SEA/
      |-- fastp/
         |-- sample1_trimmed_failed_out.fastq.gz
         |-- sample1_trimmed_merged.fastq.gz
         |-- sample1_trimmed_out1.fastq.gz
         |-- sample1_trimmed_out2.fastq.gz
         |-- sample1_trimmed_unpaired.fastq.gz
      |-- bwa/
         |-- sample1_trimmed_merged.fastq.gz.sai
         |-- sample1_trimmed_out1.fastq.gz.sai
         |-- sample1_trimmed_out2.fastq.gz.sai
      |-- bam/
         |-- sample1_paired-reads.bam
         |-- sample1_merged.bam
         |-- sample1.bam
         |-- sample1_mapped.bam
         |-- sample1_mapped-md.bam
         |-- sample1_mapped-md.bam.bai
         |-- sample1_mapped-md.bam.metrics
      |-- mapDamage/
      |-- logs/
         |-- SEA_fast.log
         |-- SEA_bwa-aln.log
         |-- SEA_bam.log
         |-- SEA_mark-duplicates.log
      |-- statistics/
         |-- SEA_stats.txt
```
Statistics file is generated as csv, here I show it in tsv for clarity. Since there are more than 20 headers, only few are shown

| Sample ID | Sample Age | Library | Treatment | Total Reads | Trimmed Reads | Mapped Reads | PCR Duplicates | Clonality |
| --------- | ---------- | ------- | --------- | ----------- | ------------- | ------------ | -------------- | --------- |
| 015683P | ia | LV3003058645 | U | 30847755 | 27592750 | 20872074 | 5936949 | 0.284444613 |
| 016054P | cwc | LV3003058689 | U | 48117377 | 41375776 | 33889887 | 13130966 | 0.38745977 |
| 022376B | bamed | LV3003058656 | E | 34732559 | 29481065 | 23598566 | 15221736 | 0.35497199 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |

# R plots and data visualizations

Box plots;

![Clonality Boxplot](https://github.com/pcichowicz/USER-Reduction/assets/81156946/6f7856f4-03d6-443a-83f4-46f03155b2f5)

mapDamage misincorporation frequency plot (mapDamage creates plot, but this plot is custom made to my liking);

![Custom R script to create mapDamage misincorporation at 5p' (C > T) and 3p' (G > A)](/021130P_mapDamage_misincorpotation.png)

USER treated vs non-USER treated, C > T misincorporation at the first 3 positions;

![all_EU](https://github.com/pcichowicz/USER-Reduction/assets/81156946/8b0667e9-aa48-435e-a3ed-710709656f53)

USER treated vs non-USER treated, using other South East Asia samples for comparison at the first position

![one_EU_sea](https://github.com/pcichowicz/USER-Reduction/assets/81156946/a80cbdd6-5872-4e11-b3e1-aef33108066d)

Ridge plots showing the 3 treatments, mapDamage lambda parameter;

![lambda_para](https://github.com/pcichowicz/USER-Reduction/assets/81156946/963a393b-34d8-430b-ae2a-d9f673fe5112)

Finally a scatter plot, relationship of the DeltaD mapDamage parameter and internal C > T mean damage frequency;

![internal_deltad](https://github.com/pcichowicz/USER-Reduction/assets/81156946/ed35da45-ada7-4a8e-97a0-52206c4b55fd)





