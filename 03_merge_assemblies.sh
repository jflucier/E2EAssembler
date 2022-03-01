#!/bin/bash

set -e

# load and valdiate env
source ${E2EAssembler}/E2EAssembler.config
${E2EAssembler}/00_check_environment.sh

export TELOMOTIF_RC=$(perl -e '
my $seq = reverse("'$TELOMOTIF'");
$seq =~ tr/ACGTUacgtu/TGCAAtgcaa/;
print $seq;
')

if [[ $NANOPORE_FASTQ == *.fastq ]]; then
    export NANOPORE_BASE=${NANOPORE_FASTQ%.fastq}
elif [[ $NANOPORE_FASTQ == *.fastq.gz ]]; then
    export NANOPORE_BASE=${NANOPORE_FASTQ%.fastq.gz}
else
    export NANOPORE_BASE=${NANOPORE_FASTQ}
fi

export REF_ASSEMBLY_NAME=$(basename ${NANOPORE_BASE})

echo "removing previous run results in merged_assembly"
rm -fr $PWD/merged_assembly/*  2>/dev/null
mkdir -p $PWD/merged_assembly

echo "running assembly merge on $REF_ASSEMBLY_NAME"
echo "We will look for assembled contigs having telomere motif at each end."
echo "5' motif is $TELOMOTIF_RC and 3' motif $TELOMOTIF"

for f in $CANU_OUTPATH/${REF_ASSEMBLY_NAME}.*
do
    NEW_ASSEMBLY=$(basename $f)
    echo "****** running ${NEW_ASSEMBLY} ******"

    mkdir -p $PWD/merged_assembly/$NEW_ASSEMBLY

    if [ "$f" = "$CANU_OUTPATH/${REF_ASSEMBLY_NAME}.0" ]; then
        echo "Using $NEW_ASSEMBLY as reference for 1st pass"
        cp $f/${NEW_ASSEMBLY}.contigs.fasta $PWD/merged_assembly/${NEW_ASSEMBLY}/merged_${NEW_ASSEMBLY}.fasta
        MERGED_ASSEMBLY_FA=$PWD/merged_assembly/${NEW_ASSEMBLY}/merged_${NEW_ASSEMBLY}.fasta
    else
        echo "Merging $NEW_ASSEMBLY to previous incomplete contigs"
        HYBRID_ASS=$PWD/merged_assembly/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta
        SELF_ASS=$f/${NEW_ASSEMBLY}.contigs.fasta

        N50=$(${BBMAP_STATS} in=$SELF_ASS format=6 | perl -ne '
        chomp($_);
        if($_ !~ /^\#/){
            my @t = split("\t",$_);
            print $t[6] . "\n";
        }
        ')

        echo "Assembly N50=$N50"
        cd $PWD/merged_assembly/$NEW_ASSEMBLY
        ${QUICKMERGE} -l $N50 --prefix $NEW_ASSEMBLY $HYBRID_ASS $SELF_ASS
        cd ../../

        MERGED_ASSEMBLY_FA=$PWD/merged_assembly/${NEW_ASSEMBLY}/merged_${NEW_ASSEMBLY}.fasta

    fi

    perl -e '
    use Bio::SeqIO;

    my $fa_in = Bio::SeqIO->new(
        -file => "<'$MERGED_ASSEMBLY_FA'",
        -format => "fasta"
    );

    my $good_seq = Bio::SeqIO->new(
        -file   => ">>'${PWD}'/merged_assembly/'${REF_ASSEMBLY_NAME}'.complete_contigs.fasta",
        -format => "fasta"
    );

    my $bad_seq = Bio::SeqIO->new(
        -file   => ">'${PWD}'/merged_assembly/'${REF_ASSEMBLY_NAME}'.tmp.incomplete_contigs.fasta",
        -format => "fasta"
    );

    while (my $seq = $fa_in->next_seq) {
        my $seq_str = $seq->seq();
        if($seq_str =~ /'$TELOMOTIF_RC'.+'$TELOMOTIF'/ ){
            $good_seq->write_seq($seq);
        }
        else{
            $bad_seq->write_seq($seq);
        }
    }

    '

    #cat tmp.incomplete_contigs.fasta $__UNNASS_REF > incomplete_contigs.fasta
    rm -f $PWD/merged_assembly/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta
    mv  $PWD/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.incomplete_contigs.fasta ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta

    ## backup curretn assembly state
    cp ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.fasta ${PWD}/merged_assembly/${NEW_ASSEMBLY}/${REF_ASSEMBLY_NAME}.complete_contigs.fasta
    cp ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta ${PWD}/merged_assembly/${NEW_ASSEMBLY}/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta

    read good_contigs g_words g_chars <<< $(grep -e '>' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.fasta | wc )
    read bad_contigs g_words g_chars <<< $(grep -e '>' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta | wc )
    echo "${NEW_ASSEMBLY}: valid end to end scaffols with telomeres = $good_contigs"
    echo "${NEW_ASSEMBLY}: incomplete scaffolds = $bad_contigs"

done

read good_contigs g_words g_chars <<< $(grep -e '>' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.fasta | wc )
read bad_contigs g_words g_chars <<< $(grep -e '>' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta | wc )

echo "*************** ${REF_ASSEMBLY_NAME} *******************"
echo "FINAL: valid end to end contigs with telomeres = $good_contigs"
echo "FINAL: incomplete contigs = $bad_contigs"
echo "******************************************"

echo "done assembly for $REF_ASSEMBLY_NAME"
