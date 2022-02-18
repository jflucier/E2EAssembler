#!/bin/bash

set -e

# You need seqtk in your path
if ! command -v seqtk &> /dev/null
then
    echo "seqtk could not be found. Please install and put in PATH. export PATH=/path/to/seqtk:\$PATH"
    exit 1
fi

# check if genome file exists
if [[ -z "${GENOME}" ]]; then
    echo "GENOME variable must be defined: export GENOME=/path/to/genome.fa"
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
fi

# check if NANOPORE_FASTQ exists
if [[ -z "${NANOPORE_FASTQ}" ]]; then
    echo "NANOPORE_FASTQ variable must be defined: export NANOPORE_FASTQ=/path/to/nanopore_reads.fastq"
    exit 1
elif [ ! -f "$NANOPORE_FASTQ" ]; then
    echo "NANOPORE_FASTQ file $NANOPORE_FASTQ file does not exist."
    exit 1
fi

if [[ -z "${DATASET_SPLIT_COVERAGE}" ]]; then
    echo "DATASET_SPLIT_COVERAGE is not defined. Will set DATASET_SPLIT_COVERAGE to default 60X COVERAGE"
    DATASET_SPLIT_COVERAGE=60
fi

if [[ -z "${CANU_THREAD}" ]]; then
    echo "CANU_THREAD is not defined. Will set CANU_THREAD to 4 cpu"
    CANU_THREAD=4
fi

if [[ -z "${CANU_MEMORY}" ]]; then
    echo "CANU_MEMORY is not defined. Will set CANU_MEMORY to 10G"
    CANU_MEMORY=10
fi

if [[ -z "${CANU_OUTPATH}" ]]; then
    echo "CANU_OUTPATH is not defined. Will set CANU_OUTPATH to $PWD/canu_assembly/"
    CANU_OUTPATH=$PWD/canu_assembly/
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
        NANOPORE_BASE=$file_no_gz
    fi
fi

echo "FASTQ to FASTA using seqtk"
seqtk seq -a ${NANOPORE_BASE}.fastq > ${NANOPORE_BASE}.fasta

echo "Requested COVERAGE=${DATASET_SPLIT_COVERAGE}X"
NTS_PER_FILE=$(( $GENOME_SIZE * $DATASET_SPLIT_COVERAGE ))
echo "For ${DATASET_SPLIT_COVERAGE}X COVERAGE we need $NTS_PER_FILE nts per FASTA file"

echo "Extracting converage info from ${NANOPORE_BASE}.fasta"
read header_lines header_words header_chars <<< $(grep -e '>' ${NANOPORE_BASE}.fasta | wc )
read total_lines total_words total_chars file <<< $(wc ${NANOPORE_BASE}.fasta )
COVERAGE=$(( ($total_chars-$header_chars) / $GENOME_SIZE ))  #= 287X
echo "Total COVERAGE is ${COVERAGE}X"
NBR_FILES=$(( (($total_chars-$header_chars) / $NTS_PER_FILE) + 1 ))
echo "Based on total sample coverage, will generate $NBR_FILES x ${DATASET_SPLIT_COVERAGE}X fastq files"

echo "generating ${DATASET_SPLIT_COVERAGE}X dataset fastq files"
perl -e '
open(my $FH,"<'${NANOPORE_BASE}'.fastq");
my $count=0;
my $out_cnt=0;
my $out_lbl="'${NANOPORE_BASE}'.$out_cnt.fastq";
print "outputting '${NANOPORE_BASE}'.$out_cnt.fastq\n";
open(my $OUT,">$out_lbl");
while( my $h = <$FH>)  {
    my $s = <$FH>;
    my $x = <$FH>;
    my $q = <$FH>;

    $count+=length($s);
    print $OUT "$h\n";
    print $OUT "$s\n";
    print $OUT "$x\n";
    print $OUT "$q\n";

    if($count > '$NTS_PER_FILE'){
        close($OUT);

        # print "nt count = $count\n";
        # print "genome size = '$GENOME_SIZE'\n";
        # print "max nt / file = '$NTS_PER_FILE'\n";
        # my $cov = $count / '$GENOME_SIZE';
        # print "COVERAGE = $cov\n";

        $out_cnt++;
        $out_lbl="'${NANOPORE_BASE}'.$out_cnt.fastq";
        print "outputting '${NANOPORE_BASE}'.$out_cnt.fastq\n";
        open($OUT,">$out_lbl");
        $count = 0;
    }

}
close($OUT);
# print "nt count = $count\n";
# print "genome size = '$GENOME_SIZE'\n";
# print "max nt / file = '$NTS_PER_FILE'\n";
# my $cov = $count / '$GENOME_SIZE';
# print "COVERAGE = $cov\n";
'

echo "Next step is to run canu denovo assembly for each ${DATASET_SPLIT_COVERAGE}X fastq files"
CANU_GENOME_SIZE=$(bc <<<"scale=2; $GENOME_SIZE/1000000")

echo "Generating shell script for canu assembly"
echo '#!/bin/bash' > 02_exec_canu.sh
echo '
set -e

if ! command -v java \&> /dev/null
then
    echo "java could not be found. Please install and put in PATH."
    exit 1
fi

#check if canu in path or die
if ! command -v canu \&> /dev/null
then
    echo "canu could not be found. Please install and put in PATH. export PATH=/path/to/canu:\$PATH"
    exit 1
fi

mkdir -p '$CANU_OUTPATH'
for f in '${NANOPORE_BASE}'.*.fastq
do
    __fastq_file=$(basename $f)
    __fastq_group=${__fastq_file%.fastq}
    mkdir -p '$CANU_OUTPATH'/${__fastq_group}
    rm -fr  '$CANU_OUTPATH'/${__fastq_group}/*
    echo "running assembly on $f"
    echo "ouptutting resulting assembly in '$PWD'/canu_assembly/${__fastq_group}"
    canu useGrid=false \
    -p ${__fastq_group} -d '$CANU_OUTPATH'/${__fastq_group} \
    genomeSize='$CANU_GENOME_SIZE'm maxThreads='$CANU_THREAD' maxMemory='$CANU_MEMORY' \
    -nanopore $f
    echo "done assembly of $f"
done

echo "done running all assemblies"
' > 02_exec_canu.sh

echo "Generating SLURM script for canu assembly. Open and edit this script prior to using."
echo '#!/bin/bash' > 02_exec_canu.slurm.sh
echo '
#SBATCH --mail-type=END,FAIL
#SBATCH -D '$CANU_OUTPATH'
#SBATCH -o '$CANU_OUTPATH'/canu-%A.out
#SBATCH --time=requested_exec_time
#SBATCH --mem='$CANU_MEMORY'G
#SBATCH -N 1
#SBATCH -n '$CANU_THREAD'
#SBATCH -A requested_allocation
#SBATCH --mail-user=your_email@yourdomain.ca
#SBATCH -J canu

if ! command -v java \&> /dev/null
then
    echo "java could not be found. Please install and put in PATH."
    exit 1
fi

if ! command -v canu &> /dev/null
then
    echo "canu could not be found. Please install and put in PATH. export PATH=/path/to/canu:$PATH"
    exit 1
fi

export __fastq=$(ls '${NANOPORE_BASE}'.*.fastq | awk "NR==$SLURM_ARRAY_TASK_ID")

export __fastq_file=$(basename $__fastq)
export __fastq_group=${__fastq_file%.fastq}

echo "copying fastqs"
cp $__fastq $SLURM_TMPDIR/${__fastq_file}

mkdir $SLURM_TMPDIR/${__fastq_group}
echo "exec canu"
canu useGrid=false \
-p ${__fastq_group} -d $SLURM_TMPDIR/${__fastq_group}/ \
genomeSize='$CANU_GENOME_SIZE'm maxThreads='$CANU_THREAD' maxMemory='$CANU_MEMORY' \
-nanopore $SLURM_TMPDIR/${__fastq_file}

echo "copying results to ip29"
cp -r $SLURM_TMPDIR/${__fastq_group} '$CANU_OUTPATH'/

echo "done"

' > 02_exec_canu.slurm.sh

echo "** DONE **"
