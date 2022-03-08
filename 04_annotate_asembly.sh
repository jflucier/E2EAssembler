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

export REF_ASSEMBLY_NAME=$(basename ${NANOPORE_BASE})

export TELOTAG_RC=$(perl -e '
my $seq = reverse("'$TELOTAG'");
$seq =~ tr/ACGTUacgtu/TGCAAtgcaa/;
print $seq;
')

export TELOMOTIF_RC=$(perl -e '
my $seq = reverse("'$TELOMOTIF'");
$seq =~ tr/ACGTUacgtu/TGCAAtgcaa/;
print $seq;
')

mkdir -p ${PWD}/annotation_assembly
rm -f ${PWD}/annotation_assembly/*

echo "initialising sqlite db"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite "
DROP TABLE IF EXISTS assembly_annotation;
create table assembly_annotation (
    name integer,
    assembly_chr text,
    start integer,
    end integer,
    strand text
);
CREATE INDEX coord_assembly_annotation_idx on assembly_annotation(start,end);
"


echo "## annotating telotags on assembly"
perl -ne '
chomp($_);
if($_ =~ /^\>/){
  my($id,$len) = $_ =~ /^\>(.*) len=(\d+) reads=/;
  print ">".$id."_".$len."\n";
}
else{
  print $_ . "\n";
}
' ${PWD}/merged_assembly/${NEW_ASSEMBLY}/${REF_ASSEMBLY_NAME}.complete_contigs.fasta > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta


${SEQKIT} locate --threads ${LOCAL_THREAD} --max-mismatch 3 --bed \
-p $TELOTAG ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.telotag.bed

perl -ne '
chomp($_);

my($id_str,$start,$end,$motif,$x,$strand) = split("\t",$_);
my($id,$len) = split("_",$id_str);

if($start < 200){
    print "telotag\t$id\t$start\t$end\t$strand\n";
}

if($start > ($len - 200)){
    print "telotag\t$id\t$start\t$end\t$strand\n";
}
' ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.telotag.bed > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.telotag.filtered.tsv
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.telotag.filtered.tsv assembly_annotation"



echo "## annotate telomere repeats on assembly"
export TELOMERIC_REGEX=$(perl -e '
my $telomotif = "'$TELOMOTIF'";
my $telomotif_str = $telomotif x 50;
my @seeds = $telomotif_str =~ /.{'$SEED_LEN'}.?/g;
my %h = map { $_ => 1 } @seeds;
my $teloregex = join(" -p ",keys(%h));

$telomotif = "'$TELOMOTIF_RC'";
$telomotif_str = $telomotif x 50;
@seeds = $telomotif_str =~ /.{'$SEED_LEN'}.?/g;
%h = map { $_ => 1 } @seeds;
$teloregex .= " -p " . join(" -p ",keys(%h));

print $teloregex . "\n";
')

${SEQKIT} locate --threads ${LOCAL_THREAD} -r --bed \
-p $TELOMERIC_REGEX -V 0 -m 0 -P \
${PWD}/merged_assembly/${NEW_ASSEMBLY}/${REF_ASSEMBLY_NAME}.complete_contigs.fasta > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.telomere_repeats.bed

${E2EAssembler}/extract_telomeric_repeat_assembly.pl ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.telomere_repeats.bed $SEED_LEN > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.telomere_repeats.filtered.tsv
perl -ne 'print "telomere_repeats\t$_";' ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.telomere_repeats.filtered.tsv > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.telomere_repeats.filtered.tagged.tsv
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.telomere_repeats.filtered.tagged.tsv assembly_annotation"

echo "## annotate subtelomere"
cp $SUBTELO_X ${PWD}/annotation_assembly/subtelo_X.fa
cp $SUBTELO_Y ${PWD}/annotation_assembly/subtelo_Y.fa

for SUBTELO in subtelo_X subtelo_Y
do

    perl -ne '
    chomp($_);
    if($_ =~ /^\>/){
      my($id,$len) = $_ =~ /^\>(.*) len=(\d+) reads=/;
      print ">".$id."\n";
    }
    else{
      print $_ . "\n";
    }
    ' ${PWD}/merged_assembly/${NEW_ASSEMBLY}/${REF_ASSEMBLY_NAME}.complete_contigs.fasta > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta

    echo "aligning $SUBTELO subtelomere elements using clustalo"
    ${CLUSTALO} --force --threads=$LOCAL_THREAD --outfmt=st \
    --in ${PWD}/annotation_assembly/${SUBTELO}.fa \
    --out ${PWD}/annotation_assembly/${SUBTELO}.sto

    echo "running nhmmer on alignment"
    nhmmer --cpu $LOCAL_THREAD -E 1e-3 \
    --tblout ${PWD}/annotation_assembly/${SUBTELO}.sto.raw.tbl \
    ${PWD}/annotation_assembly/${SUBTELO}.sto ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta

    perl -ne '
    chomp($_);
    if($_ !~ /^\#/){
        my @t = split(/\s+/,$_);
        print "'$SUBTELO'\t" . join("\t",@t) . "\n";
    }
    ' ${PWD}/annotation_assembly/${SUBTELO}.sto.raw.tbl > ${PWD}/annotation_assembly/${SUBTELO}.sto.raw.tsv

    cut -f 1,2,8,9,13  ${PWD}/annotation_assembly/${SUBTELO}.sto.raw.tsv > ${PWD}/annotation_assembly/${SUBTELO}.sto.tsv
    # reverse start end when pn negative strand
    perl -ne '
    chomp($_);
    my($name,$chr,$start,$end,$strand) = split("\t",$_);
    if($strand eq "-"){
        print "$name\t$chr\t$end\t$start\t$strand\n";
    }
    else{
        print "$name\t$chr\t$start\t$end\t$strand\n";
    }
    ' ${PWD}/annotation_assembly/${SUBTELO}.sto.tsv > ${PWD}/annotation_assembly/${SUBTELO}.sto.reformat.tsv

    sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import ${PWD}/annotation_assembly/${SUBTELO}.sto.reformat.tsv assembly_annotation"
done

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


echo "aligning assembly on reference genome using minimap2"
$MINIMAP2 -t $LOCAL_THREAD -ax map-ont $GENOME ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta \
| $SAMTOOLS view --threads $LOCAL_THREAD -Sh -q 20 -F 2048 -F 256 \
| $SAMTOOLS sort --threads $LOCAL_THREAD -o ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.sam

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
' ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.sam > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.tsv
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.tsv assembly_mapping"

echo "generate annotation bed"
# name integer,
# assembly_chr text,
# start integer,
# end integer,
# strand text
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' "
SELECT
    am.ref_chr chr,
    aa.start,
    aa.end,
    aa.name,
    0,
    aa.strand,
    '',
    '',
    CASE
        WHEN aa.name = 'telotag' THEN '255,0,0'
        WHEN aa.name = 'telomere_repeats' THEN '255,128,0'
        WHEN aa.name = 'telomere_repeats' THEN '255,255,0'
    END
FROM
    assembly_mapping am
    join assembly_annotation aa on aa.assembly_chr=am.assembly_chr
ORDER BY am.ref_chr asc, aa.start asc, aa.end asc
limit 10;
"


# Y nhmer results needs flattenning --> NO!
#${E2EAssembler}/flatten_subtelomere_results.pl ${PWD}/annotation_assembly/${SUBTELO_NAME}.sto.raw.tsv

# echo "## annotate subtelomere elements X"
# SUBTELO_NAME=$(basename $SUBTELO_X)
#
# perl -ne '
# chomp($_);
# if($_ =~ /^\>/){
#   my($id,$len) = $_ =~ /^\>(.*) len=(\d+) reads=/;
#   print ">".$id."\n";
# }
# else{
#   print $_ . "\n";
# }
# ' ${PWD}/merged_assembly/${NEW_ASSEMBLY}/${REF_ASSEMBLY_NAME}.complete_contigs.fasta > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta
#
# echo "aligning X subtelomere elements using clustalo"
# ${CLUSTALO} --force --threads=$LOCAL_THREAD --outfmt=st \
# --in $SUBTELO_X \
# --out ${PWD}/annotation_assembly/${SUBTELO_NAME}.sto
#
# echo "running nhmmer on alignment"
# nhmmer --cpu $LOCAL_THREAD -E 1e-3 \
# --tblout ${PWD}/annotation_assembly/${SUBTELO_NAME}.sto.raw.tbl \
# ${PWD}/annotation_assembly/${SUBTELO_NAME}.sto ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.complete_contigs.reformat.fasta
#
# perl -ne '
# chomp($_);
# if($_ !~ /^\#/){
#     my @t = split(/\s+/,$_);
#     print join("\t",@t) . "\n";
# }
# ' ${PWD}/annotation_assembly/${SUBTELO_NAME}.sto.raw.tbl > ${PWD}/annotation_assembly/${SUBTELO_NAME}.sto.raw.tsv

echo "done"









#
#
# echo "combining all ${DATASET_SPLIT_COVERAGE}X assembly corrected reads with 3' telotag motif"
# cat $PWD/telomotif/${REF_ASSEMBLY_NAME}/*.${TELOTAG}.fasta > $PWD/telomotif/$REF_ASSEMBLY_NAME.${TELOTAG}.fasta
#
# # calculate telomeric sewuquence length in each read (distance between telotag and subtelomeric sequence)
# # shodes method
# # minimap2 -ax map-ont /atium/Data/projects/sam_telo_AG/Num2kbtelo11RsacCer3.fa chopped.fastq | samtools view -Sbh -bq 20 -F 2048 -F 256 | samtools sort -o chopped.primary.bam
# # samtools index chopped.primary.bam
# # bedtools bamtobed -i chopped.primary.bam > chopped.primary.bed
# # ./Telo_length_AG.R -c 2 -s telo_start_position_WT_2kb.txt -i chopped.primary.bed -o telomere_lengths.tsv
#
#
# echo "Generating telomeric seeds"
# export TELOMERIC_REGEX=$(perl -e '
# my $telomotif = "'$TELOMOTIF'";
# my $telomotif_str = $telomotif x 10;
# my @seeds = $telomotif_str =~ /.{6}.?/g;
# my $teloregex = join(" -p ",@seeds);
#
# $telomotif = "'$TELOMOTIF_RC'";
# $telomotif_str = $telomotif x 10;
# @seeds = $telomotif_str =~ /.{6}.?/g;
# $teloregex .= " -p " . join(" -p ",@seeds);
# print $teloregex . "\n";
# ')
#
# echo $TELOMERIC_REGEX
#
# echo "running seqkit locate using telomeric seeds"
# cat $PWD/telomotif/${REF_ASSEMBLY_NAME}.telomotif.fasta | \
# ${SEQKIT} locate --threads ${LOCAL_THREAD} -r --bed \
# -p $TELOMERIC_REGEX -V 0 -m 0 > $PWD/telomotif/${REF_ASSEMBLY_NAME}.telomotif.bed


# inport telo reads in db
# map telo reads to reference
# import read / chr mapping to db
# run seqkit with telo seed
