#!/bin/bash

set -e

echo "################################################################################################################"

if [[ -z "${E2EAssembler}" ]]; then
    echo "## E2EAssembler install path variable must be defined: export E2EAssembler=/path/to/E2EAssembler"
    exit 1
fi

echo "## Checking all software dependencies"

if ! command -v "${SEQTK}" &> /dev/null
then
    echo "##**** seqtk could not be found ****"
    echo "## Please install seqtk and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export SEQTK=/path/to/seqtk/seqtk"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${SEQKIT}" &> /dev/null
then
    echo "##**** SEQKIT could not be found ****"
    echo "## Please install seqkit (https://bioinf.shenwei.me/seqkit/download/) and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export SEQKIT=/path/to/seqkit"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "nucmer" &> /dev/null
then
    echo "##**** MUMMER could not be found ****"
    echo "## Please install MUMMER and put in PATH"
    echo "## export PATH=/path/to/mummer:$PATH"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v quickmerge &> /dev/null
then
    echo "##**** quickmerge could not be found ****"
    echo "## Please install quickmerge and put in PATH"
    echo "## export PATH=/path/to/quickmerge:$PATH"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${QUICKMERGE}" &> /dev/null
then
    echo "##**** quickmerge could not be found ****"
    echo "## Please install quickmerge and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export QUICKMERGE=/path/to/quickmerge/merge_wrapper.py"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${BBMAP_STATS}" &> /dev/null
then
    echo "##**** BBMAP could not be found ****"
    echo "## Please install BBMAP and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export BBMAP_STATS=/path/to/BBMAP/stats.sh"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${SAMTOOLS}" &> /dev/null
then
    echo "##**** SAMTOOLS could not be found ****"
    echo "## Please install SAMTOOLS and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export SAMTOOLS=/path/to/samtools-x.xx/samtools"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${MINIMAP2}" &> /dev/null
then
    echo "##**** MINIMAP2 could not be found ****"
    echo "## Please install MINIMAP2 and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export MINIMAP2=/path/to/minimap2/minimap2"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${BEDTOOLS}" &> /dev/null
then
    echo "##**** BEDTOOLS could not be found ****"
    echo "## Please install BEDTOOLS and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export BEDTOOLS=/path/to/bedtools/bin/bedtools"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v java &> /dev/null
then
    echo "##**** JAVA could not be found ****"
    echo "## Please install java and include in PATH"
    echo "## export PATH=/path/to/java"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${CANU}" &> /dev/null
then
    echo "##**** CANU could not be found ****"
    echo "## Please install CANU and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export CANU=/path/to/canu/build/bin/canu"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v Rscript &> /dev/null
then
    echo "##**** R could not be found ****"
    echo "## Please install R and include in PATH"
    echo "## export PATH=/path/to/R"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${CLUSTALO}" &> /dev/null
then
    echo "##**** CLUSTALO could not be found ****"
    echo "## Please install CLUSTALO and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export CLUSTALO=/path/to/bin/clustalo"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${NHMMER}" &> /dev/null
then
    echo "##**** NHMMER could not be found ****"
    echo "## Please install NHMMER and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export NHMMER=/path/to/bin/nhmmer"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${BED2BB}" &> /dev/null
then
    echo "##**** bedToBigBed could not be found ****"
    echo "## Please install UCSC userapps and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export BED2BB=/path/to/bin/bedToBigBed"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${FA2BIT}" &> /dev/null
then
    echo "##**** faToTwoBit could not be found ****"
    echo "## Please install UCSC userapps and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export FA2BIT=/path/to/bin/faToTwoBit"
    echo "##**********************************"
    echo "##"
    exit 1
fi

if ! command -v "${BAM2BED}" &> /dev/null
then
    echo "##**** bamToBed could not be found ****"
    echo "## Please install UCSC userapps and edit config file ${E2EAssembler}/E2EAssembler.config"
    echo "## Modify this line: export BAM2BED=/path/to/bin/bamToBed"
    echo "##**********************************"
    echo "##"
    exit 1
fi


echo "## checking if all E2EAssembler variables properly defined"

if [[ -z "${NANOPORE_FASTQ}" ]]; then
    echo "## FATAL: NANOPORE_FASTQ variable must be defined. To set, edit config file: export NANOPORE_FASTQ=/path/to/nanopore_reads.fastq"
    exit 1
elif [ ! -f "$NANOPORE_FASTQ" ]; then
    echo "## FATAL: $NANOPORE_FASTQ file does not exist. Please specifiy a valid path. To set, edit config file: export NANOPORE_FASTQ=/path/to/nanopore_reads.fastq"
    exit 1
fi
echo "## NANOPORE datapath: $NANOPORE_FASTQ"

if [[ -z "${GENOME}" ]]; then
    echo "## FATAL: GENOME variable must be defined. To set, edit config file: export GENOME=/path/to/genome.fa"
    exit 1
elif [ ! -f "$GENOME" ]; then
    echo "## FATAL: $GENOME file does not exist. Please specifiy a valid path. To set, edit config file: export GENOME=/path/to/genome.fa"
    exit 1
fi
echo "## GENOME file: $GENOME"
if [ ! -f "$GENOME.genomesize" ]; then
    echo "## generating $GENOME.genomesize"
    read g_header_lines g_header_words g_header_chars <<< $(grep -e '>' $GENOME | wc )
    read g_total_lines total_words g_total_chars file <<< $(wc $GENOME )
    GENOME_SIZE=$(( $g_total_chars-$g_header_chars ))
    echo "$GENOME_SIZE" > $GENOME.genomesize
else
    GENOME_SIZE=$(cat $GENOME.genomesize)
fi
echo "## GENOMESIZE file: $GENOME.genomesize"
echo "## GENOME size: $GENOME_SIZE"

if [ ! -f "$GENOME.chromsizes" ]; then
    echo "## generating $GENOME.chromsizes using samtools"
    ${SAMTOOLS} faidx $GENOME
    cut -f1,2 $GENOME.fai > $GENOME.chromsizes
fi
echo "## CHROMSIZES file: $GENOME.chromsizes"

if [[ -z "${CANU_OUTPATH}" ]]; then
    echo "##**********************************"
    echo "## WARNING: CANU_OUTPATH is not defined. To set, edit config file: export CANU_OUTPATH=<<path>>"
    echo "## Will set CANU_OUTPATH to default $PWD/canu_assembly/"
    echo "##**********************************"
    echo "##"
    export CANU_OUTPATH=$PWD/canu_assembly/
else
    echo "## CANU_OUTPATH: $CANU_OUTPATH"
fi

if [[ -z "${DATASET_SPLIT_COVERAGE}" ]]; then
    echo "##**********************************"
    echo "## WARNING: DATASET_SPLIT_COVERAGE is not defined. To set, edit config file: export DATASET_SPLIT_COVERAGE=<<int>>"
    echo "## Will set DATASET_SPLIT_COVERAGE to default 60X COVERAGE"
    echo "##**********************************"
    echo "##"
    export DATASET_SPLIT_COVERAGE=60
else
    echo "## DATASET_SPLIT_COVERAGE: $DATASET_SPLIT_COVERAGE"
fi

if [[ -z "${TELOTAG}" ]]; then
    echo "##**********************************"
    echo "## WARNING: TELOTAG is not defined. To set, edit config file: export TELOTAG=<<sequence>>"
    echo "## Will set TELOTAG to default CAGTCTACACATATTCTCTGT"
    echo "##**********************************"
    echo "##"
    export TELOTAG=ACAGAGAATATGTGTAGACTG
else
    echo "## TELOTAG: $TELOTAG"
fi

if [[ -z "${TELOMOTIF}" ]]; then
    echo "##**********************************"
    echo "## WARNING: TELOMOTIF is not defined. To set, edit config file: export TELOMOTIF=<<sequence>>"
    echo "## Will set TELOMOTIF to default yeast telomere repeat TGTGGGTGTGGTG"
    echo "##**********************************"
    echo "##"
    export TELOMOTIF=TGTGGGTGTGGTG
else
    echo "## TELOMOTIF: $TELOMOTIF"
fi

if [[ -z "${SEED_LEN}" ]]; then
    echo "##**********************************"
    echo "## WARNING: SEED_LEN is not defined. To set, edit config file: export SEED_LEN=<<sequence>>"
    echo "## Will set SEED_LEN to default 8"
    echo "##**********************************"
    echo "##"
    export SEED_LEN=8
else
    echo "## SEED_LEN: $SEED_LEN"
fi


if [[ -z "${LOCAL_THREAD}" ]]; then
    echo "##**********************************"
    echo "## WARNING: LOCAL_THREAD is not defined. To set, edit config file: export LOCAL_THREAD=<<int>>"
    echo "## Will set LOCAL_THREAD to default LOCAL_THREAD=4"
    echo "##**********************************"
    echo "##"
    export LOCAL_THREAD=4
else
    echo "## LOCAL_THREAD: $LOCAL_THREAD"
fi

if [[ -z "${LOCAL_MEMORY}" ]]; then
    echo "##**********************************"
    echo "## WARNING: LOCAL_MEMORY is not defined. To set, edit config file: export LOCAL_MEMORY=<<mem in G>>"
    echo "## Will set LOCAL_MEMORY to default LOCAL_MEMORY=8"
    echo "##**********************************"
    echo "##"
    export LOCAL_MEMORY=8
else
    echo "## LOCAL_MEMORY: $LOCAL_MEMORY"
fi

if [[ -z "${SLURM_CANU_THREAD}" ]]; then
    echo "##**********************************"
    echo "## WARNING: SLURM_CANU_THREAD is not defined. To set, edit config file: export SLURM_CANU_THREAD=<<mem in G>>"
    echo "## Will set SLURM_CANU_THREAD to default SLURM_CANU_THREAD=24"
    echo "##**********************************"
    echo "##"
    export SLURM_CANU_THREAD=24
else
    echo "## SLURM_CANU_THREAD: $SLURM_CANU_THREAD"
fi

if [[ -z "${SLURM_CANU_MEMORY}" ]]; then
    echo "##**********************************"
    echo "## WARNING: SLURM_CANU_MEMORY is not defined. To set, edit config file: export SLURM_CANU_MEMORY=<<mem in G>>"
    echo "## Will set SLURM_CANU_MEMORY to default SLURM_CANU_MEMORY=30"
    echo "##**********************************"
    echo "##"
    export SLURM_CANU_MEMORY=30
else
    echo "## SLURM_CANU_MEMORY: $SLURM_CANU_MEMORY"
fi

if [[ -z "${SLURM_ALLOCATION}" ]]; then
    echo "##**********************************"
    echo "## WARNING: SLURM_ALLOCATION is not defined. To set, edit config file: export SLURM_ALLOCATION=<<slurm_account_name>>"
    echo "## Will set SLURM_ALLOCATION to to empty string. Make sure you modify cnau slurm script if you plan to use!"
    echo "##**********************************"
    echo "##"
    export SLURM_ALLOCATION=""
else
    echo "## SLURM_ALLOCATION: $SLURM_ALLOCATION"
fi

if [[ -z "${SLURM_WALLTIME}" ]]; then
    echo "##**********************************"
    echo "## WARNING: SLURM_WALLTIME is not defined. To set, edit config file: export SLURM_WALLTIME=<<HH:MM:SS>>"
    echo "## Will set SLURM_WALLTIME to default SLURM_WALLTIME=24:00:00"
    echo "##**********************************"
    echo "##"
    export SLURM_WALLTIME=24:00:00
else
    echo "## SLURM_WALLTIME: $SLURM_WALLTIME"
fi

echo "################################################################################################################"
