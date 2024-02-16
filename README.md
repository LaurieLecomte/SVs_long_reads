# SV calling pipeline from long-read sequencing data

## Pipeline Overview

1. **Call SVs** : the 3 tools may be used independently in any order or at the same time.
* 1.1. Sniffles : scripts `01.1` to `01.4`
* 2.2. SVIM : scripts `02.1` to `02.4`
* 3.2. NanoVar : scripts `03.1` to `03.4`
4. **Merge SV calls** across callers : `04_merge_callers.sh`
5. Format merged output : `05_format_merged.sh`
6. **Filter** merged output : `06_filter_merged.sh` 


### Additional scripts

Other scripts targeting a specific step or operation conducted in one of the main scripts or allowing additional analyses are provided in the `01_scripts/utils` subdirectory.

* `01_scripts/utils/format_add_ALTseq_LR.R` : adds an explicit alternate sequence (when possible) to the merged SVs. Called by the `01_scripts/utils/format_merged.R` script featured in the `05_format_merged.sh` main script.
* `01_scripts/utils/format_merged_sample_names.R` : add unique sample names to the merged VCF, also called by the `05_format_merged.sh` main script.
* `01_scripts/utils/nanovar_add_rnames.R` : add unique sample names to the NanoVar VCF. Called by the `03.1_nanovar_call.sh` script.
* `01_scripts/utils/combined_plot_by_caller.R` : used for plotting filtered short-read SVs, called by the `01_scripts/utils/summarize_plot.sh` script (Supp. Fig. 4 from the paper [Investigating structural variant, indel and single nucleotide polymorphism differentiation between locally adapted Atlantic salmon populations using whole genome sequencing and a hybrid genomic polymorphism detection approach](https://www.biorxiv.org/content/10.1101/2023.09.12.557169v1))

Older scripts used for development or debugging purposes are stored in the `01_scripts/archive` folder for future reference if needed. These are not meant to be used in their current state and may be obsolete.

## Prerequisites

### Files

* A reference genome named `genome.fasta` and its index (.fai) in `03_genome`
* Bam files for all samples and their index. These can be soft-linked in the 04_bam folder for easier handling : if `$BAM_PATH` is the remote path to bam files, use `for file in $(ls -1 $BAM_PATH/*); do ln -s $file ./04_bam; done`. These should be named as `SAMPLEID.bam` (see sample ID list below).
* A bam files list in `02_infos`. This list can be generated with the following command, where `$BAM_DIR` is the path of the directory where bam files are located : `ls -1 $BAM_DIR/*.bam > 02_infos/bam_list.txt`
* A chromosomes list (or contigs, or sites) in `02_infos`. This list is used for parallelizing the SV calling step. It can be produced from the indexed genome file (`"$GENOME".fai`) : `less "$GENOME".fai | cut -f1 > 02_infos/chr_list.txt`
* If some chromosomes are to be excluded from the SV calling step, such as unplaced contigs, these must be listed in `02_infos/excl_chrs.txt`, which needs to be encoded in linux format AND have a newline at the end.
* A chromosomes list in the form of a bed file (`02_infos/chrs.bed`). We use the one produced by `00_prepare_regions.sh` from the [SVs_SR_pipeline](https://github.com/LaurieLecomte/SVs_short_reads/blob/main/01_scripts/00_prepare_regions.sh)
* A sample IDs list in `02_infos`, one ID per line. This list can be used for renaming bam files symlinks in `$BAM_DIR`, adjust `grep` command as required (warning : use carefully): `less 02_infos/ind_ONT.txt | while read ID; do BAM_NAME=$(ls $BAM_DIR/*.bam | grep "$ID"); mv $BAM_NAME $BAM_DIR/"$ID".bam; done` and `less 02_infos/ind_ONT.txt | while read ID; do BAM_NAME=$(ls $BAM_DIR/*.bai | grep "$ID"); mv $BAM_NAME $BAM_DIR/"$ID".bam.bai; done`


### Software
#### For Manitou users
Custom conda environments are required for running `NanoVar`, `SVIM`, `sniffles2` and `jasmine`, as these programs are not available on Manitou; See the [Conda environment preparation](#conda-environment-preparation) section below. 

#### For users working with other computing clusters and servers
The program versions specified in this pipeline refer to the versions available on IBIS' bioinformatics servers when this pipeline was built in 2021-2022, and are likely not available on all other servers. 
Please add a '#' at the beginning of each line in the `#LOAD REQUIRED MODULES` section in each script (or remove these lines), and follow the [Conda environment preparation](#conda-environment-preparation) to create custom conda environments with correct program versions and dependencies.
A R installation is also required.


## Detailed Walkthrough

For running each script, copy the `srun` command from the script's header to the terminal and adjust parameters (memory, partition, time limit) if necessary.  
The header also features a brief description of the script's contents. 


### Conda environment preparation

#### SV calling environments (`SVs_LR` + `NanoVar`)
From the main directory, run `conda create --name SVs_LR --file SVs_LR_env.txt` and `conda create --name NanoVar --file NanoVar_env.txt`

These environments are used for calling SVs and contain the following callers:
* SVIM 2.0.0
* Sniffles 2.0.7
* NanoVar 1.3.8
* bcftools 1.13

#### SV merging environment (`jasmine_1.1.5`)
From the main directory, run `conda create --name jasmine_1.1.5 --file jasmine_1.1.5_env.txt`

This environment is used for merging SVs across callers, and contains [jasmine 1.1.5](https://github.com/mkirsche/Jasmine) and bcftools 1.13.


### Main pipeline

#### 1. Prepare region files (`00_prepare_regions.sh`)

This script prepares the bed files required for specifying the regions in which SVs must be called or must not be called. It first produces a bed file from the reference fasta in order to yield : 

* A text and a bed file of excluded chromosomes or contigs


#### 2. Call SVs using 3 seperate tools


##### Sniffles (scripts 01.1 to 01.4)

Before running each script for Sniffles, activate the `SVs_LR` env: `conda activate SVs_LR`

* `01.1_sniffles_call.sh`
* `01.2_sniffles_refine.sh`
* `01.3_sniffles_merge.sh`
* `01.4_sniffles_filter_format.sh`


##### SVIM (scripts 02.1 to 02.4)
Before running each script for SVIM, activate the `SVs_LR` env: `conda activate SVs_LR`

* `02.1_svim_call.sh`
* `02.2_svim_refine.sh`
* `02.3_svim_merge.sh`
* `02.4_svim_filter_format.sh`


##### NanoVar (scripts 03.1 to 03.4)
Before running each script for NanoVar, activate the `NanoVar` env: `conda activate NanoVar`

* `03.1_nanovar_call.sh`: **warning :**NanoVar indexes are picky and will crash if other indexes are present in the genome directory, so we need to provide a genome directory in which NanoVar can add its own indexes. Genome indexing steps prior to SV calling are very LONG, so we run run these steps only once for the first sample, then we run others
* `03.2_nanovar_refine.sh`
* `03.3_nanovar_merge.sh`
* `03.4_nanovar_filter_format.sh`


#### 3. Merge SV calls across callers (`04_merge_callers.sh`)
Before running this script, activate the `jasmine_1.1.5` env (even if you are working on Manitou): `conda activate jasmine_1.1.5`

#### 4. Format merged output (`05_format_merged.sh`)

#### 5. Filter merged SVs (`06_filter_merged.sh`)
Keep SVs supported by at least 2/3 tools and larger than 50 bp.
