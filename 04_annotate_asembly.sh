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

sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' "
SELECT
    am.ref_chr || '_' || length(ac.sequence) chr,
    CASE
        WHEN am.flag = 0 THEN ac.sequence
        WHEN am.flag = 16 THEN ac.sequence_rc
    END sequence
FROM
    assembly_mapping am
    join assembly_chr ac on ac.assembly_chr=am.assembly_chr
ORDER BY cast(am.ref_chr as integer) asc;
" > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.tsv

perl -ne '
chomp($_);
my @t = split("\t",$_);
print ">" . $t[0] . "\n";
print $t[1] . "\n";
' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.tsv > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.fasta

${SEQKIT} locate --threads ${LOCAL_THREAD} --max-mismatch 3 --bed \
-p $TELOTAG ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.fasta > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telotag.bed

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
' ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telotag.bed > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telotag.filtered.tsv
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telotag.filtered.tsv assembly_annotation"

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
" > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.tsv

perl -ne '
chomp($_);
my @t = split("\t",$_);
print ">" . $t[0] . "\n";
print $t[1] . "\n";
' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.tsv > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.fasta

${SEQKIT} locate --threads ${LOCAL_THREAD} -r --bed \
-p $TELOMERIC_REGEX -V 0 -m 0 -P \
${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.fasta > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telomere_repeats.bed

${E2EAssembler}/extract_telomeric_repeat_assembly.pl ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telomere_repeats.bed $SEED_LEN > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telomere_repeats.filtered.tsv
perl -ne 'print "telomere_repeats\t$_";' ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telomere_repeats.filtered.tsv > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telomere_repeats.filtered.tagged.tsv
# filter for TELO_MINLEN > 25nt
perl -ne '
chomp($_);
my @t = split("\t",$_);
my $len = $t[3] - $t[2];
if($len >= '$TELO_MINLEN' ){
    print $_ . "\n";
}
' ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telomere_repeats.filtered.tagged.tsv > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telomere_repeats.filtered.tagged.filtered.tsv
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' ".import ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.telomere_repeats.filtered.tagged.filtered.tsv assembly_annotation"

echo "## annotate subtelomere"

for SUBTELO_ENTRY in "${ANNOTATION_FASTA[@]}" ; do
    SUBTELO=${SUBTELO_ENTRY%%:*}
    SUBTELO_FA_PATH=${SUBTELO_ENTRY#*:}

    cp $SUBTELO_FA_PATH ${PWD}/annotation_assembly/${SUBTELO}.fa
    echo "finding $SUBTELO"

    echo "aligning $SUBTELO subtelomere elements using clustalo"
    ${CLUSTALO} --force --threads=$LOCAL_THREAD --outfmt=st \
    --in ${PWD}/annotation_assembly/${SUBTELO}.fa \
    --out ${PWD}/annotation_assembly/${SUBTELO}.sto

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
    " > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.tsv

    perl -ne '
    chomp($_);
    my @t = split("\t",$_);
    print ">" . $t[0] . "\n";
    print $t[1] . "\n";
    ' ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.tsv > ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.fasta

    echo "running nhmmer on alignment"
    nhmmer --cpu $LOCAL_THREAD -E 1e-3 \
    --tblout ${PWD}/annotation_assembly/${SUBTELO}.sto.raw.tbl \
    ${PWD}/annotation_assembly/${SUBTELO}.sto ${PWD}/merged_assembly/${REF_ASSEMBLY_NAME}.tmp.fasta

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

echo "generate annotation bed"
sqlite3 $PWD/${REF_ASSEMBLY_NAME}.sqlite '.separator "\t"' "
SELECT
    am.ref_chr chr,
    aa.start,
    aa.end - 1,
    aa.name,
    0,
    aa.strand
FROM
    assembly_mapping am
    join assembly_annotation aa on aa.assembly_chr=am.ref_chr
ORDER BY am.ref_chr asc, aa.start asc, aa.end asc;
" > ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.annotation.bed

# chr: chromosome id
# chr_len: longueur total du chr
# XC_nbr: # total de XC pour ce chormosome
# XCR_nbr: # total de XCR pour ce chormosome
# Y_nbr: # total de Y pour ce chormosome
# chr_str_repr: une string representant l'organisation du chromosome. Ca ressemblerais a ceci: TEL5-XCR-ITS-XC-Y-Y-ITS-TEL3
# chr_loc_list: meme chose que chr_str_repr mais avec les coords (start:end): (1:230)-(1456:4565)-(4700:4750)-(15000:17000)-(30000:35000)-(415000:417000)-(420000:422000)-(500000:500450)

mkdir -p $PWD/report
perl ${E2EAssembler}/output_annotation_report.pl ${PWD}/annotation_assembly/${REF_ASSEMBLY_NAME}.annotation.bed > ${PWD}/report/${REF_ASSEMBLY_NAME}.annotation.report.tsv

echo "done"
