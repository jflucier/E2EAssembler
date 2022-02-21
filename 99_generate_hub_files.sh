#!/bin/bash

set -e

if ! command -v minimap2 &> /dev/null
then
    echo "minimap2 could not be found. Please install and put in PATH. export PATH=/path/to/minimap2:\$PATH"
    exit 1
fi

if ! command -v samtools &> /dev/null
then
    echo "Samtools could not be found. Please install Samtools and put in PATH. NANOPORE_NAMEort PATH=/path/to/samtools:\$PATH"
    exit 1
fi


if ! command -v bamToBed &> /dev/null
then
    echo "bedtools could not be found. Please install bedtools and put in PATH. NANOPORE_NAMEort PATH=/path/to/bedtools:\$PATH"
    exit 1
fi

if ! command -v bedSort &> /dev/null
then
    echo "bedSort and bedToBigBed could not be found. Please install UCSC command line tools (https://hgdownload.soe.ucsc.edu/admin/exe/) and put in PATH. export PATH=/path/to/ucsctools:\$PATH"
    exit 1
fi


echo "** analysing $NANOPORE_FASTQ **"

if [[ $NANOPORE_FASTQ == *.fastq ]]; then
    NANOPORE_BASE=${NANOPORE_FASTQ%.fastq}
elif [[ $NANOPORE_FASTQ == *.gz ]]; then
    file_no_gz=${NANOPORE_FASTQ%.gz}
    if [ ! -f "$file_no_gz" ]; then
        echo "unzipping fastq $NANOPORE_FASTQ"
        NANOPORE_BASE=${NANOPORE_FASTQ%.fastq.gz}
        zcat $NANOPORE_FASTQ > ${NANOPORE_BASE}.fastq
    else
        # file .fastq exists
        NANOPORE_BASE=${file_no_gz%.fastq}
    fi
fi

# check if genome file exists
if [[ -z "${GENOME}" ]]; then
    echo "GENOME variable must be defined: NANOPORE_NAMEort GENOME=/path/to/genome.fa"
    exit 1
fi

if [ ! -f "$GENOME" ]; then
    echo "Genome file $GENOME file does not exist."
    exit 1
else
    read g_header_lines g_header_words g_header_chars <<< $(grep -e '>' $GENOME | wc )
    read g_total_lines total_words g_total_chars file <<< $(wc $GENOME )
    GENOME_SIZE=$(( $g_total_chars-$g_header_chars ))
    echo "genome file is $GENOME"
    echo "genome size is $GENOME_SIZE nt"
    if [ ! -f "$GENOME.chromsizes" ]; then
        echo "generating $GENOME.chromsizes using samtools"
        samtools faidx $GENOME
        cut -f1,2 $GENOME.fai > $GENOME.chromsizes
    fi
    CHROMSIZES=$GENOME.chromsizes
    echo "CHROMSIZES file is $CHROMSIZES"
fi


NANOPORE_NAME=$(basename ${NANOPORE_BASE})

echo "cleaning previous hub files"
rm -fr $PWD/hub/${NANOPORE_NAME} 2 > /dev/null
echo "generating tracks for hub for $NANOPORE_NAME"
mkdir $PWD/hub/${NANOPORE_NAME}
minimap2 -a -t 32 $GENOME ${PWD}/${NANOPORE_NAME}.complete_contigs.fasta > $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.sam
samtools view --reference $GENOME -bS $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.sam > $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.bam
bamToBed -i $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.bam > $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.bed
bedSort $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.bed $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.sorted.bed
bedToBigBed $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.sorted.bed $CHROMSIZES $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.sorted.bb

echo "done generating hub files in $PWD/hub/${NANOPORE_NAME}/"
