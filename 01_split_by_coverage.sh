#!/bin/bash

set -e

# load and valdiate env
source ${E2EAssembler}/E2EAssembler.config
${E2EAssembler}/00_check_environment.sh

if [[ $NANOPORE_FASTQ == *.fastq ]]; then
    export NANOPORE_BASE=${NANOPORE_FASTQ%.fastq}
elif [[ $NANOPORE_FASTQ == *.fastq.gz ]]; then
    export NANOPORE_BASE=${NANOPORE_FASTQ%.fastq.gz}
else
    export NANOPORE_BASE=${NANOPORE_FASTQ}
fi

export NANOPORE_NAME=$(basename ${NANOPORE_BASE})
export GENOME_SIZE=$(cat $GENOME.genomesize)

echo "Unzipping FASTQ"
if [ ! -f "${NANOPORE_BASE}.fastq" ]; then
    echo "unzipping FASTQ"
    # zcat $NANOPORE_FASTQ > ${NANOPORE_BASE}.fastq
    pigz -p $LOCAL_THREAD -dk $NANOPORE_FASTQ
else
    echo "${NANOPORE_BASE}.fastq already found. No unzipping required."
fi

echo "FASTQ to FASTA using seqkit"
if [ ! -f "${NANOPORE_BASE}.fasta" ]; then
    ${SEQKIT} fq2fa --threads $LOCAL_THREAD $NANOPORE_FASTQ > ${NANOPORE_BASE}.fasta
else
    echo "${NANOPORE_BASE}.fasta already found."
fi

echo "Requested COVERAGE=${DATASET_SPLIT_COVERAGE}X"
NTS_PER_FILE=$(( $GENOME_SIZE * $DATASET_SPLIT_COVERAGE ))
echo "For ${DATASET_SPLIT_COVERAGE}X COVERAGE we need $NTS_PER_FILE nts per FASTA file"

read header_lines header_words header_chars <<< $(grep -e '>' ${NANOPORE_BASE}.fasta | wc )
read total_lines total_words total_chars file <<< $(wc ${NANOPORE_BASE}.fasta )
COVERAGE=$(( ($total_chars-$header_chars) / $GENOME_SIZE ))  #= 287X
echo "Total FASTQ COVERAGE is ${COVERAGE}X"
NBR_FILES=$(( (($total_chars-$header_chars) / $NTS_PER_FILE) + 1 ))
echo "Based on FASTQ coverage, will generate $NBR_FILES x fastq with ${DATASET_SPLIT_COVERAGE}X coverage"

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
        $out_cnt++;
        $out_lbl="'${NANOPORE_BASE}'.$out_cnt.fastq";
        print "outputting '${NANOPORE_BASE}'.$out_cnt.fastq\n";
        open($OUT,">$out_lbl");
        $count = 0;
    }

}
close($OUT);
'

echo "Next step is to run canu denovo assembly for each ${DATASET_SPLIT_COVERAGE}X fastq files"
CANU_GENOME_SIZE=$(bc <<<"scale=2; $GENOME_SIZE/1000000")

echo "Generating shell script for canu assembly"
mkdir -p $CANU_OUTPATH

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
    '$CANU'/canu useGrid=false \
    -p ${__fastq_group} -d '$CANU_OUTPATH'/${__fastq_group} \
    genomeSize='$CANU_GENOME_SIZE'm maxThreads='$LOCAL_THREAD' maxMemory='$LOCAL_MEMORY' \
    -nanopore $f
    echo "done assembly of $f"
done

echo "done running all assemblies"
' >> 02_exec_canu.sh

echo "Generating SLURM script for canu assembly."
echo '#!/bin/bash' > 02_exec_canu.slurm.sh
echo '
#SBATCH --mail-type=END,FAIL
#SBATCH -D '$CANU_OUTPATH'
#SBATCH -o '$CANU_OUTPATH'/canu-%A_%a.out
#SBATCH --time='$SLURM_WALLTIME'
#SBATCH --mem='$SLURM_CANU_MEMORY'G
#SBATCH -N 1
#SBATCH -n '$SLURM_CANU_THREAD'
#SBATCH -A '$SLURM_ALLOCATION'
#SBATCH -J canu

export __fastq=$(ls '${NANOPORE_BASE}'.*.fastq | awk "NR==$SLURM_ARRAY_TASK_ID")

export __fastq_file=$(basename $__fastq)
export __fastq_group=${__fastq_file%.fastq}

echo "copying fastqs"
cp $__fastq $SLURM_TMPDIR/${__fastq_file}

mkdir $SLURM_TMPDIR/${__fastq_group}
echo "exec canu"
'$CANU'/canu useGrid=false \
-p ${__fastq_group} -d $SLURM_TMPDIR/${__fastq_group}/ \
genomeSize='$CANU_GENOME_SIZE'm maxThreads='$SLURM_CANU_THREAD' maxMemory='$SLURM_CANU_MEMORY' \
-nanopore $SLURM_TMPDIR/${__fastq_file}

echo "copying results to ip29"
cp -r $SLURM_TMPDIR/${__fastq_group} '$CANU_OUTPATH'/

echo "done"

' >> 02_exec_canu.slurm.sh

echo "WARNING: Make sure you EDIT slurm script prior to using."
read sbatch_array_max header_words header_chars <<< $(ls ${NANOPORE_BASE}.*.fastq | wc )
echo "To submit to slurm: sbatch --array=1-${sbatch_array_max} 02_exec_canu.slurm.sh"

echo "** DONE **"
