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

echo "removing previous run results"
rm -fr $PWD/telomotif/${REF_ASSEMBLY_NAME}* 2>/dev/null
mkdir -p $PWD/telomotif/$REF_ASSEMBLY_NAME

echo "extracting reads with telotag on $REF_ASSEMBLY_NAME corrected reads"
for f in $CANU_OUTPATH/${REF_ASSEMBLY_NAME}.*
do
    ASSEMBLY_NAME=$(basename $f)

    # chr start
    # TACTTCGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTCAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTCCACACCCACACCACACCCACACACCAC
    #                                                  CAGTCTACACATATTCTCTGT
    echo "select reads with 5' telotag sequences for $ASSEMBLY_NAME"
    ${SEQKIT} grep --threads ${LOCAL_THREAD} -s -P -m 0 -R 1:200 -p ${TELOTAG_RC} $f/${ASSEMBLY_NAME}.correctedReads.fasta.gz > $PWD/telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.${TELOTAG_RC}.fasta

    # chr end
    # TGGTGTGGGTGTGGTGTGTGGGTGTGGAAAAAAAAAAAAAACAGAGAATATGTGTAGACTGAGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACGTAACTGAACG
    #                                         ACAGAGAATATGTGTAGACTG
    # fetch reads with telotag motif
    echo "select reads with 3' telotag sequences for $ASSEMBLY_NAME"
    ${SEQKIT} grep --threads ${LOCAL_THREAD} -s -P -m 0 -R -200:-1 \
    -p ${TELOTAG} $f/${ASSEMBLY_NAME}.correctedReads.fasta.gz > $PWD/telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.${TELOTAG}.fasta

done

echo "initialising sqlite db"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite "
DROP TABLE IF EXISTS telo_reads;
create table telo_reads (
    extr integer,
    id text,
    seq text,
    telo_trimmed_seq text
);
CREATE INDEX extr_id_telo_reads_idx on telo_reads(extr,id);

DROP TABLE IF EXISTS telo_stretch;
create table telo_stretch (
    extr integer,
    telo_read_id text,
    start integer,
    end nteger,
    len integer
);
CREATE INDEX extr_id_telo_stretch_idx on telo_stretch(extr,telo_read_id);

DROP TABLE IF EXISTS telo_reads_mapping;
create table telo_reads_mapping (
    extr integer,
    telo_read_id text,
    flag integer,
    chr_map text,
    align_pos integer,
    mapq integer
);
CREATE INDEX extr_id_telo_reads_mapping_idx on telo_reads_mapping(extr,telo_read_id);

"


for extr in 5 3
do
    if [ "$extr" = "3" ]; then
        telo=${TELOTAG}
        motif=${TELOMOTIF}
    else
        telo=${TELOTAG_RC}
        motif=${TELOMOTIF_RC}
    fi

    ###### Combine the reads that have 5' telotag ######
    echo "combining all ${DATASET_SPLIT_COVERAGE}X assembly corrected reads with ${extr}' telotag motif"
    cat $PWD/telomotif/${REF_ASSEMBLY_NAME}/*.${telo}.fasta > $PWD/telomotif/$REF_ASSEMBLY_NAME.${telo}.fasta

    echo "convert fasta to tsv"
    # fasta to tsv
    perl -ne '
    chomp($_);
    if($_ =~ /^\>/){
      my($id) = $_ =~ /^\>(.*) id=\d+$/;
      print "\n'$extr'\t" . $id . "\t";
    }
    else{
      print $_;
    }
    ' $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.fasta > $PWD/telomotif/$REF_ASSEMBLY_NAME.${telo}.tsv
    tail -n +2 "$PWD/telomotif/$REF_ASSEMBLY_NAME.${telo}.tsv" > "$PWD/telomotif/$REF_ASSEMBLY_NAME.${telo}.tsv.tmp"
    mv "$PWD/telomotif/$REF_ASSEMBLY_NAME.${telo}.tsv.tmp" "$PWD/telomotif/$REF_ASSEMBLY_NAME.${telo}.tsv"

    echo "trim telotag from all ${extr}' telotag reads sequences"
    perl -e '
    open(my $FH,"<'$PWD'/telomotif/'$REF_ASSEMBLY_NAME'.'${telo}'.tsv");
    while( my $l = <$FH>)  {
        chomp($l);
        my($extr,$id,$seq) = split("\t",$l);
        my($pre,$match,$after) = $seq =~ /^(.*)('${telo}')(.*)$/;
        if(!defined($match)){
            print STDERR "$l\n";
            die("suppose to fing telotag in sequence. Problematic sequence is with id $h");
        }

        if('$extr' == 5){
            print "$extr\t$id\t$seq\t$after\n";
        }
        else{
            print "$extr\t$id\t$seq\t$pre\n";
        }
    }
    ' > $PWD/telomotif/$REF_ASSEMBLY_NAME.${telo}.telotag_trimmed.tsv

    echo "import ${extr}' telotag reads to sqlite db"
    sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import $PWD/telomotif/$REF_ASSEMBLY_NAME.${telo}.telotag_trimmed.tsv telo_reads"

    echo "using $SEED_LEN nt as seed lenght to generate ${extr}' telomeric repeat db"
    export TELOMERIC_REGEX=$(perl -e '
    my $telomotif = "'$motif'";
    my $telomotif_str = $telomotif x 50;
    my @seeds = $telomotif_str =~ /.{'$SEED_LEN'}.?/g;
    my %h = map { $_ => 1 } @seeds;
    $teloregex .= join(" -p ",keys(%h));
    print $teloregex . "\n";
    ')

    echo "generate telomeric reads fasta from sqlitedb for seqkit"
    sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' "
    select
        id || '_' || length(telo_trimmed_seq),
        telo_trimmed_seq
    from telo_reads
    where
        extr = "${extr}";
    " > $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.tsv
    perl -e '
    open(my $FH,"<'$PWD'/telomotif/'$REF_ASSEMBLY_NAME'.'${telo}'.telotag_trimmed.reformat.tsv");
    while( my $l = <$FH>)  {
        chomp($l);
        my($h,$s) = split("\t",$l);
        print ">$h\n$s\n";
    }
    ' > $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.fasta

    echo "running seqkit locate to identify ${extr}' telomeric seeds"
    cat $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.fasta | \
    ${SEQKIT} locate --threads ${LOCAL_THREAD} -r --bed \
    -p $TELOMERIC_REGEX -V 0 -m 0 -P > $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.bed

    echo "extracting telomeric stretch from reads"
    ${E2EAssembler}/extract_telomeric_stretch.pl $extr $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.bed $SEED_LEN > $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.filtered.tsv

    echo "importing telomeric stretch to sqlite db"
    sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.filtered.tsv telo_stretch"

    # map reads to reference to determine chr
    echo "generate telomeric reads fasta from sqlitedb for minimap2"
    sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' "
    select
        id,
        telo_trimmed_seq
    from telo_reads
    where
        extr = "${extr}";
    " > $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.tsv
    perl -e '
    open(my $FH,"<'$PWD'/telomotif/'$REF_ASSEMBLY_NAME'.'${telo}'.telotag_trimmed.reformat.tsv");
    while( my $l = <$FH>)  {
        chomp($l);
        my($h,$s) = split("\t",$l);
        print ">$h\n$s\n";
    }
    ' > $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.fasta

    echo "aligning reads on reference genome using minimap2"
    $MINIMAP2 -t $LOCAL_THREAD -ax map-ont $GENOME $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.fasta \
    | $SAMTOOLS view --threads $LOCAL_THREAD -Sh -q 20 -F 2048 -F 256 \
    | $SAMTOOLS sort --threads $LOCAL_THREAD -o $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.sam

    echo "import SAM alignment information to sqlitedb"
    perl -ne '
    chomp($_);
    if($_ !~ /^\@/){
        # not header line
        my @f = split("\t",$_);
        print
            "'${extr}'\t"
            . $f[0] . "\t"
            . $f[1] . "\t"
            . $f[2] . "\t"
            . $f[3] . "\t"
            . $f[4] . "\n";
    }
    ' $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.sam > $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.tsv
    sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.tsv telo_reads_mapping"

done

# generate report graph using sqlitedb data
echo "generate relomere length report"
mkdir -p $PWD/report
## tsv report by chr: CHR, extr, median_all, avg_all, stdev_all, median_chr, avg_chr, stdev_chr
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' "
SELECT
    case trm.extr
        WHEN '5' THEN trm.chr_map || 'L'
        WHEN '3' THEN trm.chr_map || 'R'
        ELSE trm.chr_map || '?'
    END telomere,
    group_concat(ts.len)
from
    telo_reads_mapping trm
    join telo_stretch ts on ts.telo_read_id=trm.telo_read_id and ts.extr=trm.extr
GROUP BY trm.chr_map,trm.extr
ORDER BY cast(trm.chr_map as integer) ASC, trm.extr DESC
" > $PWD/report/${REF_ASSEMBLY_NAME}.telomotif.data.tsv

ALL_MEDIAN=$(sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' "
SELECT
    len
FROM telo_stretch
ORDER BY len
LIMIT 1
OFFSET (
    SELECT COUNT(*)
    FROM telo_stretch
) / 2;
")

echo -e "telomere\tall_length_median\tlength_median\tlength_average\tlength_stdev\tdatanbr\tlength_data" > $PWD/report/${REF_ASSEMBLY_NAME}.telomotif2.tsv
perl -ne '
use Statistics::Basic qw(:all nofill);
chomp($_);
my($chr_lbl,$data_str) = split("\t",$_);
my @lens = split(",",$data_str);
my @sorted = sort { $a <=> $b } @lens;
my $nbr = scalar(@lens);
my $v1  = vector(@lens);
my $med = median($v1);
my $avg = mean($v1);
my $stddev = stddev($v1);
print "$chr_lbl\t'$ALL_MEDIAN'\t$med\t$avg\t$stddev\t$nbr\t".join(",",@sorted)."\n";
' $PWD/report/${REF_ASSEMBLY_NAME}.telomotif.data.tsv >> $PWD/report/${REF_ASSEMBLY_NAME}.telomotif.tsv
rm $PWD/report/${REF_ASSEMBLY_NAME}.telomotif.data.tsv

sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' '.header on' "
SELECT
    ts.len length,
    case trm.extr
        WHEN '5' THEN trm.chr_map || 'L'
        WHEN '3' THEN trm.chr_map || 'R'
        ELSE trm.chr_map || '?'
    END telomere,
    '"${REF_ASSEMBLY_NAME}"' genotype
from
    telo_reads_mapping trm
    join telo_stretch ts on ts.telo_read_id=trm.telo_read_id and ts.extr=trm.extr
ORDER BY cast(trm.chr_map as integer) ASC, trm.extr DESC
" > $PWD/report/${REF_ASSEMBLY_NAME}.telomotif.rdata.tsv

echo "ouputting distribution plots"
Rscript $PWD/report/${REF_ASSEMBLY_NAME}.telomotif.rdata.tsv

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
