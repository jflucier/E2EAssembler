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

echo "sqlitedb setup to import and analyse telomeric reads"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite "
DROP TABLE IF EXISTS genome_info;
create table genome_info (
    ref_chr text,
    len integer
);
"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import $GENOME.chromsizes genome_info"

echo "removing previous run results"
rm -fr $PWD/telomotif/${REF_ASSEMBLY_NAME}* 2>/dev/null
mkdir -p $PWD/telomotif/$REF_ASSEMBLY_NAME

echo "extracting reads with telotag on $REF_ASSEMBLY_NAME corrected reads"
f=$CANU_OUTPATH/${REF_ASSEMBLY_NAME}


ASSEMBLY_NAME=$(basename $f)
FA_IN=$f/${ASSEMBLY_NAME}.correctedReads.fasta.gz

echo "changing orientation to be same as reference genome"

echo "unzipping fasta $FA_IN"
pigz -dc -p $LOCAL_THREAD $FA_IN > $f/${ASSEMBLY_NAME}.correctedReads.fasta

echo "sqlitedb setup to import selected reads"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite "
DROP TABLE IF EXISTS tmp_telo_reads;
create table tmp_telo_reads (
    id text,
    sequence text,
    sequence_rc text
);

DROP TABLE IF EXISTS tmp_reads_mapping;
create table tmp_reads_mapping (
    id text,
    flag integer,
    ref_chr text,
    align_pos integer,
    mapq integer
);
"

echo "fasta to tsv"
perl -ne '
chomp($_);
if($_ =~ /^\>/){
  my($id) = $_ =~ /^\>(.*) id=/;
  print "\n".$id."\t";
}
else{
  print $_;
}
' $f/${ASSEMBLY_NAME}.correctedReads.fasta > $f/${ASSEMBLY_NAME}.correctedReads.fasta.tsv
sed -i '1d' $f/${ASSEMBLY_NAME}.correctedReads.fasta.tsv

echo "adding reverse complement to reads tsv"
perl -ne '
chomp($_);
my @t = split("\t",$_);
my $rc_seq = reverse($t[1]);
$rc_seq =~ tr/ACGTUacgtu/TGCAAtgcaa/;
print $t[0] . "\t" . $t[1] . "\t" . $rc_seq . "\n";
' $f/${ASSEMBLY_NAME}.correctedReads.fasta.tsv > $f/${ASSEMBLY_NAME}.correctedReads.fasta.rc.tsv

echo "importing reads sequence + rc"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import $f/${ASSEMBLY_NAME}.correctedReads.fasta.rc.tsv tmp_telo_reads"

### align using minimap
echo "aligning reads on reference genome to determine orientation"
$MINIMAP2 -t $LOCAL_THREAD -ax map-ont $GENOME $f/${ASSEMBLY_NAME}.correctedReads.fasta \
| $SAMTOOLS view --threads $LOCAL_THREAD -Sh -q 20 -F 2048 -F 256 \
| $SAMTOOLS sort --threads $LOCAL_THREAD -o $f/${ASSEMBLY_NAME}.correctedReads.fasta.sam

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
' $f/${ASSEMBLY_NAME}.correctedReads.fasta.sam > $f/${ASSEMBLY_NAME}.correctedReads.fasta.sam.tsv

echo "importing reads mapping on ref genome with original sequence"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import $f/${ASSEMBLY_NAME}.correctedReads.fasta.sam.tsv tmp_reads_mapping"

echo "regen read fasta with correct read sequence orientation"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' "
SELECT
    am.id,
    CASE
        WHEN am.flag = 0 THEN ac.sequence
        WHEN am.flag = 16 THEN ac.sequence_rc
    END
FROM
    tmp_reads_mapping am
    join tmp_telo_reads ac on ac.id=am.id;
" > $PWD/telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.correctedReads.ok.tsv

echo "fasta to tsv"
perl -ne '
chomp($_);
my @t = split("\t",$_);
print ">" . $t[0] . "\n";
print $t[1] . "\n";
' $PWD/telomotif/${REF_ASSEMBLY_NAME}/${ASSEMBLY_NAME}.correctedReads.ok.tsv > $PWD/telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.correctedReads.ok.fasta

# chr start
# TACTTCGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTCAGTCTACACATATTCTCTGTTTTTTTTTTTTTTTTTTTTTTTCCACACCCACACCACACCCACACACCAC
#                                                  CAGTCTACACATATTCTCTGT
echo "select reads with 5' telotag sequences for $ASSEMBLY_NAME"
${SEQKIT} grep --threads ${LOCAL_THREAD} -s -P -m 1 -R 1:200 \
-p ${TELOTAG_RC} $PWD/telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.correctedReads.ok.fasta > $PWD/telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.${TELOTAG_RC}.fasta

# chr end
# TGGTGTGGGTGTGGTGTGTGGGTGTGGAAAAAAAAAAAAAACAGAGAATATGTGTAGACTGAGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACGTAACTGAACG
#                                         ACAGAGAATATGTGTAGACTG
# fetch reads with telotag motif
echo "select reads with 3' telotag sequences for $ASSEMBLY_NAME"
${SEQKIT} grep --threads ${LOCAL_THREAD} -s -P -m 1 -R -200:-1 \
-p ${TELOTAG} $PWD/telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.correctedReads.ok.fasta > $PWD/telomotif/$REF_ASSEMBLY_NAME/${ASSEMBLY_NAME}.${TELOTAG}.fasta



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
      my($id) = $_ =~ /^\>(.*)$/;
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
            die("suppose to fing telotag in sequence. Problematic sequence is with id $h");
        }
        elsif('$extr' == 5){
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
    $MINIMAP2 -t $LOCAL_THREAD -ax map-ont ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.fasta $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.fasta \
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

    sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite "
    DROP TABLE IF EXISTS telo_reads_mapping_"$extr";
    create table telo_reads_mapping_"$extr" (
        extr integer,
        telo_read_id text,
        flag integer,
        chr_map text,
        align_pos integer,
        mapq integer
    );

    "
    sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import $PWD/telomotif/${REF_ASSEMBLY_NAME}.${telo}.telotag_trimmed.reformat.tsv telo_reads_mapping_"$extr""

done

${SAMTOOLS} faidx ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.fasta
cut -f1,2 ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.fasta.fai > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.fasta.chromsizes

sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite "
DROP TABLE IF EXISTS assembly_info;
create table assembly_info (
    ref_chr text,
    len integer
);
"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.fasta.chromsizes assembly_info"

sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite "
DROP TABLE IF EXISTS telo_reads_mapping;

CREATE TABLE telo_reads_mapping as
WITH q5 AS (
select
    telo_read_id,
    count(*) c
FROM telo_reads_mapping_5
GROUP BY 1
HAVING c >= 2
),
q3 as (
select
    telo_read_id,
    count(*) c
FROM telo_reads_mapping_3
GROUP BY 1
HAVING c >= 2
)
SELECT
    *
from telo_reads_mapping_5
where
    telo_read_id not in (select telo_read_id from q5)
UNION
SELECT
    *
from telo_reads_mapping_3
where
    telo_read_id not in (select telo_read_id from q3);

CREATE INDEX extr_id_telo_reads_mapping_idx on telo_reads_mapping(extr,telo_read_id);
"

# generate report graph using sqlitedb data
echo "generate telomere length report (reads filtered for minimal telomeric stretch of $TELO_MINLEN and aligning within $READ_TELOTAG_EXTR_CUTOFF nt of chromosome extremity)"
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
    join genome_info g on trm.chr_map=g.ref_chr
    join telo_reads tr on tr.id=trm.telo_read_id and tr.extr=trm.extr
WHERE
    ts.len >= "$TELO_MINLEN"
    AND (
    trm.align_pos <= "$READ_TELOTAG_EXTR_CUTOFF"
    or (trm.align_pos + length(tr.telo_trimmed_seq)) >= (g.len - "$READ_TELOTAG_EXTR_CUTOFF")
    )
GROUP BY trm.chr_map,trm.extr
ORDER BY cast(trm.chr_map as integer) ASC, trm.extr DESC
" > $PWD/report/${REF_ASSEMBLY_NAME}.telomotif.data.tsv

ALL_MEDIAN=$(sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' "
SELECT
    ts.len
FROM
    telo_reads_mapping trm
    join telo_stretch ts on ts.telo_read_id=trm.telo_read_id and ts.extr=trm.extr
    join genome_info g on trm.chr_map=g.ref_chr
    join telo_reads tr on tr.id=trm.telo_read_id and tr.extr=trm.extr
where
    ts.len >= "$TELO_MINLEN"
    AND (
    trm.align_pos <= "$READ_TELOTAG_EXTR_CUTOFF"
    or (trm.align_pos + length(tr.telo_trimmed_seq)) >= (g.len - "$READ_TELOTAG_EXTR_CUTOFF")
    )
ORDER BY ts.len
LIMIT 1
OFFSET (
    SELECT COUNT(*)
    FROM
        telo_reads_mapping trm
        join telo_stretch ts on ts.telo_read_id=trm.telo_read_id and ts.extr=trm.extr
        join genome_info g on trm.chr_map=g.ref_chr
        join telo_reads tr on tr.id=trm.telo_read_id and tr.extr=trm.extr
    where
        ts.len >= "$TELO_MINLEN"
        AND (
        trm.align_pos <= "$READ_TELOTAG_EXTR_CUTOFF"
        or (trm.align_pos + length(tr.telo_trimmed_seq)) >= (g.len - "$READ_TELOTAG_EXTR_CUTOFF")
        )
) / 2;
")

echo -e "telomere\tall_length_median\tlength_median\tlength_average\tlength_stdev\tdatanbr\tlength_data" > $PWD/report/${REF_ASSEMBLY_NAME}.telomotif.tsv
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
    join genome_info g on trm.chr_map=g.ref_chr
    join telo_reads tr on tr.id=trm.telo_read_id and tr.extr=trm.extr
WHERE
    ts.len >= "$TELO_MINLEN"
    AND (
    trm.align_pos <= "$READ_TELOTAG_EXTR_CUTOFF"
    or (trm.align_pos + length(tr.telo_trimmed_seq)) >= (g.len - "$READ_TELOTAG_EXTR_CUTOFF")
    )
ORDER BY cast(trm.chr_map as integer) ASC, trm.extr DESC
" > $PWD/report/${REF_ASSEMBLY_NAME}.telomotif.rdata.tsv

echo "filtering chromosome telomere length where there is only a single data point"
perl -e '
open(my $FH,"<'$PWD'/report/'${REF_ASSEMBLY_NAME}'.telomotif.rdata.tsv");
my @lines = <$FH>;
chomp(@lines);
my %struct;
foreach my $l (@lines){
    if($l =~ /^length\t/){
        next;
    }
    my($len,$chr,$lbl) = split("\t",$l);
    if(!exists($struct{$chr})){
        $struct{$chr} = [];
    }
    push(@{$struct{$chr}},$l);
}

print "length\ttelomere\tgenotype\n";
foreach my $c (keys(%struct)){
    if(scalar(@{$struct{$c}}) > 1){
        foreach my $i (@{$struct{$c}}){
            print $i . "\n";
        }
    }
    else{
        print STDERR "#### CHR $c removed since it contains only a single length data: $i\n"
    }
}
' > $PWD/report/${REF_ASSEMBLY_NAME}.telomotif.rdata.filtered.tsv


echo "ouputting distribution plots"
Rscript ${E2EAssembler}/gen_telo_length.R $PWD/report/${REF_ASSEMBLY_NAME}.telomotif.rdata.filtered.tsv

echo "done"
