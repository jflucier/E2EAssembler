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

# need mummer in path for quickmerge soft
MUMMER_PATH=$(whereis nucmer | perl -ane '
chomp($_);
my @a=split(":");
my $s = substr($a[1],1);
print $s . "\n";'
)
export PATH=${MUMMER_PATH}:$PATH

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

# N50=$(${BBMAP_STATS} in=$GENOME format=6 | perl -ne '
# chomp($_);
# if($_ !~ /^\#/){
#     my @t = split("\t",$_);
#     print $t[6] . "\n";
# }
# ')
# echo "Genome N50=$N50"

for f in $CANU_OUTPATH/${REF_ASSEMBLY_NAME}.*
do
    NEW_ASSEMBLY=$(basename $f)
    echo "****** running ${NEW_ASSEMBLY} ******"

    mkdir -p $PWD/merged_assembly/$NEW_ASSEMBLY

    perl -e '
    use Bio::SeqIO;

    my $fa_in = Bio::SeqIO->new(
        -file => "<'$f'/'${NEW_ASSEMBLY}'.contigs.fasta",
        -format => "fasta"
    );

    my $good_seq = Bio::SeqIO->new(
        -file   => ">'$f'/'${NEW_ASSEMBLY}'.complete_contigs.fasta",
        -format => "fasta"
    );

    while (my $seq = $fa_in->next_seq) {
        my $seq_str = $seq->seq();
        if($seq_str =~ /'$TELOMOTIF_RC'/ ){
            $good_seq->write_seq($seq);
        }
        elsif($seq_str =~ /'$TELOMOTIF'/ ){
            $good_seq->write_seq($seq);
        }
    }
    '

    echo "## BEFORE MERGE"
    read good_contigs g_words g_chars <<< $(grep -e '>' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.fasta | wc )
    read bad_contigs g_words g_chars <<< $(grep -e '>' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta | wc )
    echo "${NEW_ASSEMBLY}: valid end to end scaffols with telomeres = $good_contigs"
    echo "${NEW_ASSEMBLY}: incomplete scaffolds = $bad_contigs"


    if [ "$f" = "$CANU_OUTPATH/${REF_ASSEMBLY_NAME}.0" ]; then
        echo "Using $NEW_ASSEMBLY as reference for 1st pass"
        #cp $f/${NEW_ASSEMBLY}.contigs.fasta $PWD/merged_assembly/${NEW_ASSEMBLY}/merged_${NEW_ASSEMBLY}.fasta
        cp $f/${NEW_ASSEMBLY}.complete_contigs.fasta $PWD/merged_assembly/${NEW_ASSEMBLY}/merged_${NEW_ASSEMBLY}.fasta
        MERGED_ASSEMBLY_FA=$PWD/merged_assembly/${NEW_ASSEMBLY}/merged_${NEW_ASSEMBLY}.fasta
    else
        echo "Merging $NEW_ASSEMBLY to previous incomplete contigs"
        HYBRID_ASS=$PWD/merged_assembly/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta
        SELF_ASS=$f/${NEW_ASSEMBLY}.complete_contigs.fasta

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

    ## backup current assembly state
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

cat ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.incomplete_contigs.fasta >> ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.fasta

echo "initialising annotation in sqlite db"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite "
DROP TABLE IF EXISTS assembly_mapping;
create table assembly_mapping (
    assembly_chr text,
    flag integer,
    ref_chr text,
    align_pos integer,
    mapq integer,
    new_ass_chr_name text
);
"

perl -ne '
chomp($_);
if($_ =~ /^\>/){
  my($id,$len,$r) = $_ =~ /^\>(tig\d+)_?\s?len=(\d+)_?\s?reads=(\d+)/;
  print ">" . $id . "_". $len . "_". $r . "\n";
  #print "|" . $id . "|\n";
}
else{
  print $_ . "\n";
}
' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.fasta > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta

perl -e '
open(my $FA, "<'${PWD}'/merged_assembly/'${REF_ASSEMBLY_NAME}'.complete_contigs.reformat.fasta");
my @lines = <$FA>;
chomp(@lines);
my $c = 1;
foreach my $l (@lines){
    if($l =~ /^\>/){
        print ">tmp_" . $c . "\n";
        $c++;
    }
    else{
        print $l . "\n";
    }
}
' > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat2.fasta

echo "aligning assembly on reference genome using minimap2"
$MINIMAP2 -t $LOCAL_THREAD -ax map-ont $GENOME ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat2.fasta \
| $SAMTOOLS view -Sh -q 20 -F 2048 -F 256 \
| $SAMTOOLS sort -o ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.sam

echo "import SAM alignment information to sqlitedb"
perl -ne '
chomp($_);
if($_ !~ /^\@/){
    # not header line
    my @f = split("\t",$_);
    print
        $f[0] . "\t"
        . $f[1] . "\t"
        . $f[2] . "\t"
        . $f[3] . "\t"
        . $f[4] . "\n";
}
' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.sam > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.sam.tsv

echo "renaming assembly chr based on alignment"
perl -e '
open(my $FH, "<'${PWD}'/merged_assembly/'${REF_ASSEMBLY_NAME}'.complete_contigs.reformat.sam.tsv");
my @lines = <$FH>;
chomp(@lines);
my %struct;
foreach my $l (@lines){
    my($assembly_chr,$flag,$ref_chr,$align_pos,$mapq) = split("\t",$l);
    if(!exists($struct{$ref_chr})){
        $struct{$ref_chr} = [];
    }
    push(@{$struct{$ref_chr}},$l);
}

foreach my $chr (sort(keys(%struct))){
    #print "$chr\n";
    if(scalar(@{$struct{$chr}}) == 1){
        my($assembly_chr,$flag,$ref_chr,$align_pos,$mapq) = split("\t",$struct{$chr}->[0]);
        #print "$ref_chr\t$flag\t$ref_chr\t$align_pos\t$mapq\n";
        print $struct{$chr}->[0] . "\t" . $ref_chr . "\n";
    }
    else{
        my $c = 1;
        foreach my $l (@{$struct{$chr}}){
            my($assembly_chr,$flag,$ref_chr,$align_pos,$mapq) = split("\t",$l);
            #print $ref_chr . "_" . $c . "\t$flag\t$ref_chr\t$align_pos\t$mapq\n";
            print $l . "\t" . $ref_chr . "_" . $c . "\n";
            $c++;
        }
    }
}
' > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.sam.renamed.tsv
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.sam.renamed.tsv assembly_mapping"

echo "generate annotation fasta using same chr name as reference genome"
perl -ne '
chomp($_);
if($_ =~ /^\>/){
  my($id) = $_ =~ /^\>(.*)$/;
  print "\n".$id."\t";
}
else{
  print $_;
}
' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat2.fasta > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat2.fasta.tsv
sed -i '1d' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat2.fasta.tsv

echo "adding reverse complement sequence to assembly"
perl -ne '
chomp($_);
my @t = split("\t",$_);
my $rc_seq = reverse($t[1]);
$rc_seq =~ tr/ACGTUacgtu/TGCAAtgcaa/;
print $t[0] . "\t" . $t[1] . "\t" . $rc_seq . "\n";
' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat2.fasta.tsv > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat2.fasta.rc.tsv

echo "adding new assembly to sqlitedb"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite "
DROP TABLE IF EXISTS assembly_chr;
create table assembly_chr (
    assembly_chr text,
    sequence text,
    sequence_rc text
);
"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat2.fasta.rc.tsv assembly_chr"

echo "output new assembly with new chr names based on reference genome"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' "
SELECT
    am.new_ass_chr_name chr,
    CASE
        WHEN am.flag = 0 THEN ac.sequence
        WHEN am.flag = 16 THEN ac.sequence_rc
    END sequence
FROM
    assembly_mapping am
    join assembly_chr ac on ac.assembly_chr=am.assembly_chr
ORDER BY cast(am.new_ass_chr_name as integer) asc;
" > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tsv

perl -ne '
chomp($_);
my @t = split("\t",$_);
print ">" . $t[0] . "\n";
print $t[1] . "\n";
' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tsv > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.fasta

echo "done assembly for $REF_ASSEMBLY_NAME"
echo "FASTA output: ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.fasta"
