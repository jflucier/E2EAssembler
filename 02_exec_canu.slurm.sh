#!/bin/bash

#SBATCH --mail-type=END,FAIL
#SBATCH -D /home/jflucier/localhost/projet/E2EAssembler/canu_assembly
#SBATCH -o /home/jflucier/localhost/projet/E2EAssembler/canu_assembly/canu-%A.out
#SBATCH --time=24:00:00
#SBATCH --mem=251G
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -A def-mundy7
#SBATCH -J canu

export __fastq=$(ls /home/jflucier/Documents/service/externe/wellinger/20211116_subtelomere_analysis/data/21029_bar6_wt.*.fastq | awk "NR==$SLURM_ARRAY_TASK_ID")

export __fastq_file=$(basename $__fastq)
export __fastq_group=${__fastq_file%.fastq}

echo "copying fastqs"
cp $__fastq $SLURM_TMPDIR/${__fastq_file}

mkdir $SLURM_TMPDIR/${__fastq_group}
echo "exec canu"
/home/jflucier/app/canu/build/bin/canu/canu useGrid=false \
-p ${__fastq_group} -d $SLURM_TMPDIR/${__fastq_group}/ \
genomeSize=12.21m maxThreads=48 maxMemory=251 \
-nanopore $SLURM_TMPDIR/${__fastq_file}

echo "copying results to ip29"
cp -r $SLURM_TMPDIR/${__fastq_group} /home/jflucier/localhost/projet/E2EAssembler/canu_assembly/

echo "done"


