#!/bin/bash

set -e

if ! command -v seqtk &> /dev/null
then
    echo "seqtk could not be found. Please install and put in PATH. export PATH=/path/to/seqtk:\$PATH"
    exit 1
fi

if ! command -v merge_wrapper.py &> /dev/null
then
    echo "quickmerge could not be found. Please install quickmerge (https://github.com/mahulchak/quickmerge) and put in PATH. export PATH=/path/to/quickmerge:\$PATH"
    exit 1
fi

if ! command -v stats.sh &> /dev/null
then
    echo "BBMap could not be found. Please install BBMap and put in PATH. export PATH=/path/to/bbmap:\$PATH"
    exit 1
fi

if ! command -v samtools &> /dev/null
then
    echo "Samtools could not be found. Please install Samtools and put in PATH. export PATH=/path/to/bbmap:\$PATH"
    exit 1
fi

if ! command -v minimap2 &> /dev/null
then
    echo "minimap2 could not be found. Please install and put in PATH. export PATH=/path/to/minimap2:\$PATH"
    exit 1
fi

if ! command -v bamToBed &> /dev/null
then
    echo "bedtools could not be found. Please install bedtools and put in PATH. NANOPORE_NAMEort PATH=/path/to/bedtools:\$PATH"
    exit 1
fi

if ! command -v bedtools &> /dev/null
then
    echo "bedtools could not be found. Please install bedtools and put in PATH. NANOPORE_NAMEort PATH=/path/to/bedtools:\$PATH"
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
    echo "GENOME variable must be defined: REF_ASSEMBLY_NAMEort GENOME=/path/to/genome.fa"
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

if [[ -z "${CANU_OUTPATH}" ]]; then
    echo "CANU_OUTPATH is not defined. Will set CANU_OUTPATH to $PWD/canu_assembly/"
    CANU_OUTPATH=$PWD/canu_assembly/
fi

if [[ -z "${TELOTAG5}" ]]; then
    echo "TELOTAG5 is not defined. Will set TELOTAG5 to CAGTCTACACATATTCTCTGT"
    TELOTAG5=CAGTCTACACATATTCTCTGT
fi

if [[ -z "${TELOTAG3}" ]]; then
    echo "TELOTAG5 is not defined. Will set TELOTAG5 to ACAGAGAATATGTGTAGACTG"
    TELOTAG3=ACAGAGAATATGTGTAGACTG
fi


echo "removing previous run results in merged_assembly"
rm -fr merged_assembly/*  2>/dev/null
mkdir -p merged_assembly

NANOPORE_NAME=$(basename ${NANOPORE_BASE})
REF_ASSEMBLY_DIR=${CANU_OUTPATH}/${NANOPORE_NAME}.0
REF_ASSEMBLY_NAME_PART=$(basename $REF_ASSEMBLY_DIR)
REF_ASSEMBLY_NAME=${REF_ASSEMBLY_NAME_PART%.0}
echo "cleaning previous runs temp files for $REF_ASSEMBLY_NAME"
rm $REF_ASSEMBLY_NAME.tmp.incomplete_contigs.fasta $REF_ASSEMBLY_NAME.complete_contigs.fasta $REF_ASSEMBLY_NAME.incomplete_contigs.fasta

echo "running assembly merge on $REF_ASSEMBLY_NAME"
for f in $CANU_OUTPATH/${REF_ASSEMBLY_NAME}.*
do
    echo "****** running $f ******"
    ASSEMBLY_NAME=$(basename $f)
    mkdir -p merged_assembly/$ASSEMBLY_NAME

    if [ "$f" = "$REF_ASSEMBLY_DIR" ]; then
        echo "Using $f as reference for 1st pass"
        MERGED_ASSEMBLY_FA=$REF_ASSEMBLY_DIR/${ASSEMBLY_NAME}.contigs.fasta

    else
        echo "Merging $ASSEMBLY_NAME to previous incomplete contings"
        HYBRID_ASS=$PWD/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta
        SELF_ASS=$f/${ASSEMBLY_NAME}.contigs.fasta

        N50=$(stats.sh in=$SELF_ASS format=6 | perl -ne '
        chomp($_);
        if($_ !~ /^\#/){
            my @t = split("\t",$_);
            print $t[6] . "\n";
        }
        ')

        echo "Assembly N50=$N50"
        cd merged_assembly/$ASSEMBLY_NAME
        /ip29/jflucier/service/externe/wellinger/program/quickmerge/merge_wrapper.py -l $N50 --prefix $ASSEMBLY_NAME ../../$HYBRID_ASS ../../$SELF_ASS
        cd ../../

        MERGED_ASSEMBLY_FA=merged_assembly/${ASSEMBLY_NAME}/merged_${ASSEMBLY_NAME}.fasta

    fi

    perl -e '
    use Bio::SeqIO;

    my $fa_in = Bio::SeqIO->new(
        -file => "<'$MERGED_ASSEMBLY_FA'",
        -format => "fasta"
    );

    my $good_seq = Bio::SeqIO->new(
        -file   => ">>'${PWD}'/'${REF_ASSEMBLY_NAME}'.complete_contigs.fasta",
        -format => "fasta"
    );

    my $bad_seq = Bio::SeqIO->new(
        -file   => ">'${PWD}'/'${REF_ASSEMBLY_NAME}'.tmp.incomplete_contigs.fasta",
        -format => "fasta"
    );

    while (my $seq = $fa_in->next_seq) {
        my $seq_str = $seq->seq();
        if($seq_str =~ /'${TELOTAG5}'.+'${TELOTAG3}'/ ){
            $good_seq->write_seq($seq);
        }
        else{
            $bad_seq->write_seq($seq);
        }
    }

    '

    #cat tmp.incomplete_contigs.fasta $__UNNASS_REF > incomplete_contigs.fasta
    rm  ${PWD}/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta
    mv  ${PWD}/${REF_ASSEMBLY_NAME}.tmp.incomplete_contigs.fasta ${PWD}/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta

    ## backup curretn assembly state
    cp ${PWD}/${REF_ASSEMBLY_NAME}.complete_contigs.fasta merged_assembly/$ASSEMBLY_NAME/
    cp ${PWD}/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta merged_assembly/${ASSEMBLY_NAME}/

    read good_contigs g_words g_chars <<< $(grep -e '>' ${REF_ASSEMBLY_NAME}.complete_contigs.fasta | wc)
    read bad_contigs g_words g_chars <<< $(grep -e '>' ${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta | wc)
    echo "valid end to end scaffols with telomeres = $good_contigs"
    echo "incomplete scaffolds = $bad_contigs"

    # echo "select telomeric sequences for $f"
    # seqkit grep -s -R -200:-1 -r -p ${TELOTAG3} canu_assembly/${ASSEMBLY_NAME}/${ASSEMBLY_NAME}.correctedReads.fasta.gz > telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.filtered.GTGTG.fastq
    # seqkit grep -s -R 1:200 -r -p ${TELOTAG5} canu_assembly/${ASSEMBLY_NAME}/${ASSEMBLY_NAME}.correctedReads.fasta.gz > telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.filtered.CACAC.fastq
    # ##Combine the reads that are tailed
    # cat telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.filtered.GTGTG.fastq telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.filtered.CACAC.fastq > telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.filtered.telomotif.fastq

done

read good_contigs g_words g_chars <<< $(grep -e '>' ${PWD}/${REF_ASSEMBLY_NAME}.complete_contigs.fasta | wc)
read bad_contigs g_words g_chars <<< $(grep -e '>' ${PWD}/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta | wc)
echo "*************** ${REF_ASSEMBLY_NAME} *******************"
echo "valid end to end scaffols with telomeres = $good_contigs"
echo "incomplete scaffolds = $bad_contigs"
echo "******************************************"

echo "cleaning previous hub files"
rm -fr $PWD/hub/${NANOPORE_NAME} 2 > /dev/null
echo "generating tracks for hub for $NANOPORE_NAME"
mkdir $PWD/hub/${NANOPORE_NAME}
minimap2 -a -t 32 $GENOME ${PWD}/${NANOPORE_NAME}.complete_contigs.fasta > $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.sam
samtools view --reference $GENOME -bS $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.sam > $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.bam
bedtools bamtobed -i $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.bam > $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.bed
bedtools sort $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.fasta.bed $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.sorted.bed
bedToBigBed $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.sorted.bed $CHROMSIZES $PWD/hub/${NANOPORE_NAME}/${NANOPORE_NAME}.complete_contigs.sorted.bb

echo "done generating hub files in $PWD/hub/${NANOPORE_NAME}/"

echo "done assembly for $REF_ASSEMBLY_NAME"
