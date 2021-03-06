# E2EAssembler 1.0 User Manual

End to end chromosome assembly for long read sequencing data.

----

## Contents ##

* [Requirements](#requirements)
* [Installation](#initial-installation)
* [How to run](#how-to-run)
* [Output files](#output-files)

----

## Requirements ##

1. [seqkit](https://bioinf.shenwei.me/seqkit/) (version >= 2.1.0)
2. [quickmerge](https://github.com/mahulchak/quickmerge) (version >= 0.3)
3. [bbmap](https://sourceforge.net/projects/bbmap/)
4. [samtools](http://www.htslib.org/) (version >= 1.15)
5. [minimap2](https://github.com/lh3/minimap2)
6. [bedtools](https://github.com/arq5x/bedtools2/releases) (version >= 2.30.0)
7. [canu](https://github.com/marbl/canu/releases) (version >= v2.3-development)
8. [clustalo](http://www.clustal.org/omega/) (version >= 1.2.4)
9. [hmmer](http://hmmer.org/) (version >= 3.3)
10. [UCSC tools](https://hgdownload.soe.ucsc.edu/admin/exe/)
11. [mummer](http://mummer.sourceforge.net/)

Please install the required software in a location of your choice.

Please note that MUMMER and quinckmerge executables must also be in your path.
```
export PATH=/path/to/mummer:$PATH
export PATH=/path/to/quickmerge:$PATH
```

Before using E2EAssembler, make sure you edit the configuration file (see install section).

----

## Initial Installation ##

To install E2EAssembler you need to:

* Create a clone of the repository:

    ``$ git clone https://github.com/jflucier/E2EAssembler.git ``

    Note: Creating a clone of the repository requires [Github](https://github.com/) to be installed.

----

## How to run ##

E2EAssembler was developped by execution steps.

First you need to define software variables that will be used by pipeline. I recommend copying and modifying example file from the E2EAssembler install folder to your working directory and edit configuration files based on your environment.

```

cd /path/to/working_dir
cp /path/to/E2EAssembler/my.example.config .

```

The first step is to split long read fastq based on whole genome coverage. In the configuration file, a parameter DATASET_SPLIT_COVERAGE
is defined to 60X by default. The algoithm will try to optimise split for requested coverage. It will try to maximise coverage and read
number in split fastq based on original coverage found in provided fastq file.

```

bash ${E2EAssembler}/01_split_by_coverage.sh my.example.config

# loading and validating env
# ################################################################################################################
# ## Checking all software dependencies
# ## checking if all E2EAssembler variables properly defined
# ## NANOPORE datapath: /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/data/211029_bar8_sir3.fastq.gz
# ## GENOME file: /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/analysis/211029_bar7_sir1/GCA_000146045.2_R64_genomic.fna
# ## GENOMESIZE file: /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/analysis/211029_bar7_sir1/# GCA_000146045.2_R64_genomic.fna.genomesize
# ## GENOME size: 12222226
# ## CHROMSIZES file: /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/analysis/211029_bar7_sir1/# GCA_000146045.2_R64_genomic.fna.chromsizes
# ## CANU_OUTPATH: /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/analysis/211029_bar8_sir3/canu_assembly
# ## DATASET_SPLIT_COVERAGE: 60
# ## TELOTAG: ACAGAGAATATGTGTAGACTG
# ## TELOMOTIF: TGTGGGTGTGGTG
# ## SEED_LEN: 6
# ## LOCAL_THREAD: 4
# ## LOCAL_MEMORY: 8
# ## SLURM_CANU_THREAD: 48
# ## SLURM_CANU_MEMORY: 251
# ## SLURM_ALLOCATION: def-mundy7
# ## SLURM_WALLTIME: 72:00:00
# ################################################################################################################
# /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/data/211029_bar8_sir3.fastq already found. No unzipping required.
# FASTQ to FASTA using seqkit
# /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/data/211029_bar8_sir3.fasta already found.
# Requested COVERAGE=60X
# Total FASTQ COVERAGE is 415X
# Based on FASTQ coverage, will generate 6 x fastq with 69.1747638278003X coverage
# outputting /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/data/211029_bar8_sir3.0.fastq
# outputting /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/data/211029_bar8_sir3.1.fastq
# outputting /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/data/211029_bar8_sir3.2.fastq
# outputting /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/data/211029_bar8_sir3.3.fastq
# outputting /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/data/211029_bar8_sir3.4.fastq
# outputting /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/data/211029_bar8_sir3.5.fastq
# Next step is to run canu denovo assembly for each 69.1747638278003X fastq files
# Generating shell script for canu assembly
# Generating SLURM script for canu assembly.
# WARNING: Make sure you EDIT slurm script prior to using.
# To execute locally: bash /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/analysis/211029_bar8_sir3/canu_assembly/02_exec_canu.slurm.sh
# -- OR --
# To submit to slurm: sbatch --array=1-7 /nfs3_ib/ip29-ib/ip29/jflucier/service/externe/wellinger/analysis/20211116_subtelomere_analysis/analysis/211029_bar8_sir3/canu_assembly/# 02_exec_canu.slurm.sh
# ** DONE **


```

The following step will execute Canu assembly software on each of the splitted fastq files. You can either run it locally on your computer (long!) or submit on a slurm HPC cluster.

```
# Script are outputted in folder specified in vairable CANU_OUTPATH in configuration file
# to run canu assembly locally:
bash ${CANU_OUTPATH}/02_exec_canu.slurm.sh

# OR

# you can submit the canu assembly step on your slurm cluster using sbatch command as outputted by previous step execution log
# As mentionned, first edit 02_exec_canu.slurm.sh with correct #SBATCH parameters based on your compute allocation
sbatch ${CANU_OUTPATH}/02_exec_canu.slurm.sh

```

Next step is to merge the produced assemblies together

```

bash ${E2EAssembler}/03_merge_assemblies.sh my.example.config

```

Next step is to annotate assembly with telomere repeats, subtelomeres and telotag.

Make sure the configuration variable ANNOTATION_FASTA is well defined. A fasta for each subtelomeres family of sequences to be identified must be provided.
The bash associative array must bu in the format: label:/path/to/fasta

```
## for example, in config file:

export ANNOTATION_FASTA=(
    "Y:/home/jflucier/Documents/service/externe/wellinger/20211116_subtelomere_analysis/test/yeast_subtelo_Y.fa"
    "XC:/home/jflucier/Documents/service/externe/wellinger/20211116_subtelomere_analysis/test/X_XC_subtelo.fa"
    "XCR:/home/jflucier/Documents/service/externe/wellinger/20211116_subtelomere_analysis/test/XCR_subtelo.fa"
)

```

As you can see in the above example,

* 3 fasta are provided in the followinf format name:/path/to/fasta.
* Y subtelomeric elements located /home/jflucier/Documents/service/externe/wellinger/20211116_subtelomere_analysis/test/yeast_subtelo_Y.fa
* XC subtelomeric elements located /home/jflucier/Documents/service/externe/wellinger/20211116_subtelomere_analysis/test/X_XC_subtelo.fa
* XCR subtelomeric elements located /home/jflucier/Documents/service/externe/wellinger/20211116_subtelomere_analysis/test/XCR_subtelo.fa

To start this analysis

```

bash ${E2EAssembler}/04_annotate_asembly.sh my.example.config

```

Next available step is used to calculte telomere length based on reads.

```

bash ${E2EAssembler}/05_reads_telomere_length_distribution.split.sh my.example.config

```

Finally, you can optionally generate trackhubs files

```

bash ${E2EAssembler}/99_generate_hub_files.sh my.example.config

```

----

## Output files ##

When all steps of E2EAssembler are run, main output files will be created in these folders:

### 1. Merged assembly files (step 3 analysis) ###

* Final assembly files are found in the merged_assembly folder: <<my_sample>>.tsv and <<my_sample>>.fasta.
* The fasta file contains all assembled chromosome where telomeric lotifs where found at start and end of sequence.
* The tsv file is tab seperated file with the following columns: chromosome name and assembled chromosome sequence

### 2. Merged assembly annotation files (step 4 analysis) ###

* A file named *.annotation.bed contains all indentified features in bed format is generated in annotation_assembly folder
* A tab seperated file named *.annotation.report.tsv is generated in the report folder. This file contains the following columns:
  * chr: chromosome name based on reference genome
  * *_cnt: List of subtelomeric element count. The name of column is based on the provided fasta list in config file parameter variable ANNOTATION_FASTA
  * telomere_repeats_cnt: # of telomeric repeats identified
  * telotag_cnt: # of telotag identified
  * organisation_str: A string representing chromosomal organisation in the format subtelo_name1-subtelo_name2-subtelo_name...
  * organisation_coords: A string representing chromosomal organisation start positions in the format (subtelo_name1.start:subtelo_name1.end)-(subtelo_name2.start:subtelo_name2.end)-(subtelo_name3.start:subtelo_name3.end)...

### 3. Telomere length from reads files (step 5 analysis) ###

* All files generated by this step are located in the report folder
* *.pdf: A pdf file showing telomeric length distribution by chromosome
* *.telomotif.tsv: A tsv file is tab seperated file with the following columns:
  * telomere: telomere tag representing a chromosome arm. For example, 1L means chromosome 1 left arm (5') and 1R means chromosome 1 right arm (3')
  * all_length_median: A median length of all telomere identified reads (all chromosome arms)
  * length_median: Median telomere length for chromosome arm
  * length_average: Average telomere length for chromosome arm
  * length_stdev: Standar deviation of telomere length for chromosome arm
  * datanbr: Number of read identified with telomere
  * length_data: Length value extracted for each reads

### 4. Hub files (step 99 analysis) ###

* All hub files are locaed in hub folder.
* See UCSC Genome browser documentation to setup hub [here](https://genome.ucsc.edu/goldenPath/help/hubQuickStart.html)


## Run analysis directly on original fastq ( NO fastq splitting -- EXPERIMENTAL) ##

```

### no splitting analysis
bash ${E2EAssembler}/01_build_assembly_script.sh my.example.config
bash ${E2EAssembler}/03_generate_assemblies.sh my.example.config
bash ${E2EAssembler}/04_annotate_asembly.sh my.example.config
bash ${E2EAssembler}/05_reads_telomere_length_distribution.sh my.example.config
bash ${E2EAssembler}/99_generate_hub_files.sh my.example.config


```
