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

f=$CANU_OUTPATH/${REF_ASSEMBLY_NAME}
NEW_ASSEMBLY=$(basename $f)
mkdir -p $PWD/merged_assembly/${NEW_ASSEMBLY}/
cp $f/${NEW_ASSEMBLY}.contigs.fasta $PWD/merged_assembly/${NEW_ASSEMBLY}/merged_${NEW_ASSEMBLY}.fasta
MERGED_ASSEMBLY_FA=$PWD/merged_assembly/${NEW_ASSEMBLY}/merged_${NEW_ASSEMBLY}.fasta

# polish
ln -s $f/${NEW_ASSEMBLY}.contigs.fasta $PWD/merged_assembly/${NEW_ASSEMBLY}/raw_reads.fasta
ln -s $MERGED_ASSEMBLY_FA $PWD/merged_assembly/${NEW_ASSEMBLY}/contigs.fasta

echo "#### polishing assembly (remove duplicate contigs) ####"
MUMMER_PATH=$(whereis mummer | perl -ane '
chomp($_);
my @a=split(":");
my $s = substr($a[1],1);
print $s . "\n";'
)
/usr/bin/python2.7 ${FINISHERSC} \
-par $LOCAL_THREAD $PWD/merged_assembly/${NEW_ASSEMBLY}/ ${MUMMER_PATH}

perl -e '
use Bio::SeqIO;

my $fa_in = Bio::SeqIO->new(
    -file => "<'$PWD'/merged_assembly/'${NEW_ASSEMBLY}'/improved3.fasta",
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
echo "*************** ${REF_ASSEMBLY_NAME} *******************"
echo "FINAL: valid end to end contigs with telomeres = $good_contigs"
echo "FINAL: incomplete contigs = $bad_contigs"
echo "******************************************"


echo "initialising annotaton in sqlite db"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite "
DROP TABLE IF EXISTS assembly_mapping;
create table assembly_mapping (
    assembly_chr text,
    flag integer,
    ref_chr text,
    align_pos integer,
    mapq integer
);
"

perl -e '
open(my $FA, "<'${PWD}'/merged_assembly/'${REF_ASSEMBLY_NAME}'.complete_contigs.fasta");
my @lines = <$FA>;
chomp(@lines);
my $c = 1;
foreach my $l (@lines){
    if($l =~ /^\>/){
        print $l . "_" . $c . "\n";
        $c++;
    }
    else{
        print $l . "\n";
    }
}
' > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta

echo "aligning assembly on reference genome using minimap2"
$MINIMAP2 -t $LOCAL_THREAD -ax map-ont $GENOME ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta \
| $SAMTOOLS view --threads $LOCAL_THREAD -Sh -q 20 -F 2048 -F 256 \
| $SAMTOOLS sort --threads $LOCAL_THREAD -o ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.sam

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
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.sam.tsv assembly_mapping"

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
' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta.tsv
sed -i '1d' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta.tsv

echo "adding reverse complement sequence to assembly"
perl -ne '
chomp($_);
my @t = split("\t",$_);
my $rc_seq = reverse($t[1]);
$rc_seq =~ tr/ACGTUacgtu/TGCAAtgcaa/;
print $t[0] . "\t" . $t[1] . "\t" . $rc_seq . "\n";
' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta.tsv > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta.rc.tsv

echo "adding new assembly to sqlitedb"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite "
DROP TABLE IF EXISTS assembly_chr;
create table assembly_chr (
    assembly_chr text,
    sequence text,
    sequence_rc text
);
"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta.rc.tsv assembly_chr"

echo "output new assembly with new chr names based on reference genome"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' "
SELECT
    am.ref_chr chr,
    CASE
        WHEN am.flag = 0 THEN ac.sequence
        WHEN am.flag = 16 THEN ac.sequence_rc
    END sequence
FROM
    assembly_mapping am
    join assembly_chr ac on ac.assembly_chr=am.assembly_chr
ORDER BY cast(am.ref_chr as integer) asc;
" > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tsv

perl -ne '
chomp($_);
my @t = split("\t",$_);
print ">" . $t[0] . "\n";
print $t[1] . "\n";
' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tsv > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.fasta

echo "done assembly for $REF_ASSEMBLY_NAME"
echo "FASTA output: ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.fasta"
