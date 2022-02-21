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


if ! command -v bedtools &> /dev/null
then
    echo "bedtools could not be found. Please install bedtools and put in PATH. NANOPORE_NAMEort PATH=/path/to/bedtools:\$PATH"
    exit 1
fi

if ! command -v seqkit &> /dev/null
then
    echo "seqkit could not be found. Please install seqkit and put in PATH. export PATH=/path/to/seqkit:\$PATH"
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

if [[ -z "${TELOTAG5}" ]]; then
    echo "TELOTAG5 is not defined. Will set TELOTAG5 to CAGTCTACACATATTCTCTGT"
    TELOTAG5=CAGTCTACACATATTCTCTGT
fi

if [[ -z "${TELOTAG3}" ]]; then
    echo "TELOTAG5 is not defined. Will set TELOTAG5 to ACAGAGAATATGTGTAGACTG"
    TELOTAG3=ACAGAGAATATGTGTAGACTG
fi



NANOPORE_NAME=$(basename ${NANOPORE_BASE})
REF_ASSEMBLY_DIR=${CANU_OUTPATH}/${NANOPORE_NAME}.0
REF_ASSEMBLY_NAME_PART=$(basename $REF_ASSEMBLY_DIR)
REF_ASSEMBLY_NAME=${REF_ASSEMBLY_NAME_PART%.0}

echo "removing previous run results in telomere"
rm -fr telomotif/$REF_ASSEMBLY_NAME  2>/dev/null
mkdir -p telomotif/$REF_ASSEMBLY_NAME

echo "running assembly merge on $REF_ASSEMBLY_NAME"
for f in $CANU_OUTPATH/${REF_ASSEMBLY_NAME}.*
do
    echo "****** running $f ******"
    ASSEMBLY_NAME=$(basename $f)

    # chr start
    # TACTTCGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTCAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTCCACACCCACACCACACCCACACACCAC
    #                                                  CAGTCTACACATATTCTCTGT
    #
    # chr end
    # TGGTGTGGGTGTGGTGTGTGGGTGTGGAAAAAAAAAAAAAACAGAGAATATGTGTAGACTGAGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACGTAACTGAACG
    #                                         ACAGAGAATATGTGTAGACTG

    # fetch reads with telotag motif
    mkdir -p telomotif/
    echo "select telomeric sequences for $f"
    seqkit grep -s -R -200:-1 -r -p ${TELOTAG3} canu_assembly/${ASSEMBLY_NAME}/${ASSEMBLY_NAME}.correctedReads.fasta.gz > telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.${TELOTAG3}.fastq
    seqkit grep -s -R 1:200 -r -p ${TELOTAG5} canu_assembly/${ASSEMBLY_NAME}/${ASSEMBLY_NAME}.correctedReads.fasta.gz > telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.${TELOTAG5}.fastq
    ##Combine the reads that are tailed
    cat telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.${TELOTAG3}.fastq telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.${TELOTAG5}.fastq > telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.telomotif.fastq

done

# calculate telomeric sewuquence length in each read (distance between telotag and subtelomeric sequence)
# shodes method
minimap2 -ax map-ont /atium/Data/projects/sam_telo_AG/Num2kbtelo11RsacCer3.fa chopped.fastq | samtools view -Sbh -bq 20 -F 2048 -F 256 | samtools sort -o chopped.primary.bam
samtools index chopped.primary.bam
bedtools bamtobed -i chopped.primary.bam > chopped.primary.bed
./Telo_length_AG.R -c 2 -s telo_start_position_WT_2kb.txt -i chopped.primary.bed -o telomere_lengths.tsv


### testing
telo=../telomotif/21029_bar6_wt/all.filtered.telomotif.fastq
seqkit grep -s -R -200:-1 -r -p ${TELOTAG3} canu_assembly/${ASSEMBLY_NAME}/${ASSEMBLY_NAME}.correctedReads.fasta.gz > telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.${TELOTAG3}.fastq
seqkit grep -s -R 1:200 -r -p ${TELOTAG5} canu_assembly/${ASSEMBLY_NAME}/${ASSEMBLY_NAME}.correctedReads.fasta.gz > telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.${TELOTAG5}.fastq

echo "FASTQ to FASTA using seqtk"
seqtk seq -a ${telo} > ${telo}.fasta
seqkit grep -s -R -200:-1 -r -p ${TELOTAG3} ../${telo}.fasta > telo3.fasta
seqkit grep -s -R 1:200 -r -p ${TELOTAG5} ../${telo}.fasta > telo5.fasta

perl -e '
open(my $FH, "<telo5.fasta");
my @l = <$FH>;
my $c=0;
my $printseq=0;
foreach my $li (@l){
    if($c == 1 and $li =~/^\>/){
        print $li;
        $printseq = 1;
        $c++;
    }
    elsif($c == 0 and $li =~/^\>/){
        $c++;
    }
    elsif($c > 1 and $li =~/^\>/){
        exit;
    }
    elsif($printseq){
        print $li;
    }
}

' > telo5.2.fasta



perl -e '
open(my $FH, "<telo5.2.fasta");
my @l = <$FH>;
chomp(@l);
my $c=0;
my $seq="";
foreach my $li (@l){
    if($li =~/^\>/){
        print $li . "\n";
        $c++;
    }
    elsif($c > 1 and $li =~/^\>/){
        last;
    }
    else{
        $seq = $seq . $li;
    }
}
print "$seq\n";

my @nts

'

echo "done"
