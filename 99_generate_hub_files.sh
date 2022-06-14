#!/bin/bash

set -e

echo "loading and validating env"
export E2EAssembler=$(dirname "$0")
if [ -z ${1+x} ]; then
    echo "Please provide a configuration file. See ${E2EAssembler}/my.example.config for an example."
    exit 1
fi

# load and valdiate env
source $1
${E2EAssembler}/00_check_environment.sh

if [[ $NANOPORE_FASTQ == *.fastq ]]; then
    export NANOPORE_BASE=${NANOPORE_FASTQ%.fastq}
elif [[ $NANOPORE_FASTQ == *.fastq.gz ]]; then
    export NANOPORE_BASE=${NANOPORE_FASTQ%.fastq.gz}
else
    export NANOPORE_BASE=${NANOPORE_FASTQ}
fi

export REF_ASSEMBLY_NAME=$(basename ${NANOPORE_BASE})

mkdir -p ${PWD}/hub

#create hub file
echo "hub E2EAssembler
shortLabel E2EAssembler
longLabel E2EAssembler hub
genomesFile genomes.txt
email myEmail@address
" > ${PWD}/hub/hub.txt

mkdir -p ${PWD}/hub/${REF_ASSEMBLY_NAME}
cp ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.fasta ${PWD}/hub/${REF_ASSEMBLY_NAME}/
echo "generating 2 bit file of assembly for hub"
${FA2BIT} ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.fasta ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.2bit
${SAMTOOLS} faidx ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.fasta
cut -f1,2 ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.fasta.fai > ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.chromsizes

#groups $REF_ASSEMBLY_NAME/groups.txt
echo "
genome $REF_ASSEMBLY_NAME
trackDb $REF_ASSEMBLY_NAME/trackDb.txt
description $REF_ASSEMBLY_NAME E2Eassembly results hub
twoBitPath $REF_ASSEMBLY_NAME/$REF_ASSEMBLY_NAME.2bit
organism $REF_ASSEMBLY_NAME
defaultPos 1:1-5000
orderKey 4700
" > ${PWD}/hub/genomes.txt

GENOME_NAME=$(basename $GENOME)
echo "generating ref chromosome track on assembly"
$MINIMAP2 -t $LOCAL_THREAD -ax map-ont ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.fasta $GENOME \
| $SAMTOOLS view -Sh -q 20 -F 2048 -F 256 \
| $SAMTOOLS sort -o ${PWD}/hub/${REF_ASSEMBLY_NAME}/${GENOME_NAME}.sam
$SAMTOOLS view --reference ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.fasta -bS ${PWD}/hub/${REF_ASSEMBLY_NAME}/${GENOME_NAME}.sam > ${PWD}/hub/${REF_ASSEMBLY_NAME}/${GENOME_NAME}.bam
$BAM2BED -i ${PWD}/hub/${REF_ASSEMBLY_NAME}/${GENOME_NAME}.bam > ${PWD}/hub/${REF_ASSEMBLY_NAME}/${GENOME_NAME}.bed
${BEDTOOLS} sort -i ${PWD}/hub/${REF_ASSEMBLY_NAME}/${GENOME_NAME}.bed > ${PWD}/hub/${REF_ASSEMBLY_NAME}/${GENOME_NAME}.sort.bed
${BED2BB} ${PWD}/hub/${REF_ASSEMBLY_NAME}/${GENOME_NAME}.sort.bed ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.chromsizes ${PWD}/hub/${REF_ASSEMBLY_NAME}/${GENOME_NAME}.sort.bb

echo "generating annotation track on assembly"
cp ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.annotation.bed ${PWD}/hub/${REF_ASSEMBLY_NAME}/
${BEDTOOLS} sort -i ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.annotation.bed > ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.annotation.sorted.bed
${BED2BB} ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.annotation.sorted.bed ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.chromsizes ${PWD}/hub/${REF_ASSEMBLY_NAME}/${REF_ASSEMBLY_NAME}.annotation.sorted.bb
echo "
track ${REF_ASSEMBLY_NAME}_annotation
longLabel ${REF_ASSEMBLY_NAME} telomere annotation
shortLabel ${REF_ASSEMBLY_NAME}_annot
bigDataUrl ${REF_ASSEMBLY_NAME}.annotation.sorted.bb
type bigBed 6
html alignment

track $GENOME_NAME
longLabel ${GENOME_NAME} reference genome alignment
shortLabel ${GENOME_NAME}
bigDataUrl ${GENOME_NAME}.sort.bb
type bigBed 6
html alignment
" > ${PWD}/hub/${REF_ASSEMBLY_NAME}/trackDb.txt


echo "done generating hub files in ${PWD}/hub/"
