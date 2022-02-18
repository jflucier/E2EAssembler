#!/bin/bash

set -e

if ! command -v seqtk &> /dev/null
then
    echo "seqtk could not be found. Please install and put in PATH. export PATH=/path/to/seqtk:\$PATH"
    exit 1
fi

if ! command -v merge_wrapper.py &> /dev/null
then
    echo "quickmerge could not be found. Please install quickmerge (https://github.com/mahulchak/quickmerge) and put in PATH. export PATH=/path/to/quickmerge:\$PATH"
    exit 1
fi

if ! command -v stats.sh &> /dev/null
then
    echo "BBMap could not be found. Please install BBMap and put in PATH. export PATH=/path/to/bbmap:\$PATH"
    exit 1
fi

if ! command -v samtools &> /dev/null
then
    echo "Samtools could not be found. Please install Samtools and put in PATH. export PATH=/path/to/bbmap:\$PATH"
    exit 1
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
        # file .fastq exists
        NANOPORE_BASE=${file_no_gz%.fastq}
    fi
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
    if [ ! -f "$GENOME.chromsizes" ]; then
        echo "generating $GENOME.chromsizes using samtools"
        samtools faidx $GENOME
        cut -f1,2 $GENOME.fai > $GENOME.chromsizes
    fi
    CHROMSIZES=$GENOME.chromsizes
    echo "CHROMSIZES file is $CHROMSIZES"
fi

if [[ -z "${CANU_OUTPATH}" ]]; then
    echo "CANU_OUTPATH is not defined. Will set CANU_OUTPATH to $PWD/canu_assembly/"
    CANU_OUTPATH=$PWD/canu_assembly/
fi


echo "removing previous run results in merged_assembly"
rm -fr merged_assembly/*  2>/dev/null
mkdir -p merged_assembly

for d in ${CANU_OUTPATH}/*.0
do
    ref_name=$(basename $d)
    exp=${ref_name%.fasta.0}
    echo "cleaning previous runs temp files for $exp"
    rm $exp.tmp.incomplete_contigs.fasta $exp.complete_contigs.fasta $exp.incomplete_contigs.fasta

    echo "running assembly merge on $exp"
    for f in canu_assembly/${exp}.*
    do
        echo "****** running $f ******"
        name=$(basename $f)
        mkdir -p merged_assembly/$name
        mkdir -p telomotif/$exp/

        if [ "$f" = "$d" ]; then
            echo "Using $f as reference for 1st pass"
            merged_ass_fasta=$d/${name}.contigs.fasta

        else
            echo "Merging $name to previous incomplete contings"
            HYBRID_ASS=${exp}.incomplete_contigs.fasta
            SELF_ASS=$f/${name}.contigs.fasta

            N50=$(stats.sh in=$SELF_ASS format=6 | perl -ne '
            chomp($_);
            if($_ !~ /^\#/){
                my @t = split("\t",$_);
                print $t[6] . "\n";
            }
            ')

            echo "Assembly N50=$N50"
            cd merged_assembly/$name
            /ip29/jflucier/service/externe/wellinger/program/quickmerge/merge_wrapper.py -l $N50 --prefix $name ../../$HYBRID_ASS ../../$SELF_ASS
            cd ../../

            merged_ass_fasta=merged_assembly/${name}/merged_${name}.fasta
            #merge_wrapper.py hybrid_assembly.fasta self_assembly.fasta
            # nucmer -l 100 -prefix merged_assembly/$d $SELF_ASS $HYBRID_ASS
            # delta-filter -i 95 -r -q merged_assembly/$d.delta > merged_assembly/$d.rq.delta
            # quickmerge -d merged_assembly/$d.rq.delta -q complete_contigs.fasta -r $__REF -hco 5.0 -c 1.5 -l n -p merged_assembly/merged_assembly

        fi

        perl -e '
        use Bio::SeqIO;

        my $fa_in = Bio::SeqIO->new(
            -file => "<'$merged_ass_fasta'",
            -format => "fasta"
        );

        my $good_seq = Bio::SeqIO->new(
            -file   => ">>'${exp}'.complete_contigs.fasta",
            -format => "fasta"
        );

        my $bad_seq = Bio::SeqIO->new(
            -file   => ">'${exp}'.tmp.incomplete_contigs.fasta",
            -format => "fasta"
        );

        while (my $seq = $fa_in->next_seq) {
            my $seq_str = $seq->seq();
            if($seq_str =~ /CACACCCACACAC.+GTGTGTGGGTGTG/ ){
                $good_seq->write_seq($seq);
            }
            else{
                $bad_seq->write_seq($seq);
            }
        }

        '

        #cat tmp.incomplete_contigs.fasta $__UNNASS_REF > incomplete_contigs.fasta
        rm  ${exp}.incomplete_contigs.fasta
        mv  ${exp}.tmp.incomplete_contigs.fasta ${exp}.incomplete_contigs.fasta

        ## backup curretn assembly state
        cp ${exp}.complete_contigs.fasta merged_assembly/$name/
        cp ${exp}.incomplete_contigs.fasta merged_assembly/${name}/

        read good_contigs g_words g_chars <<< $(grep -e '>' ${exp}.complete_contigs.fasta | wc)
        read bad_contigs g_words g_chars <<< $(grep -e '>' ${exp}.incomplete_contigs.fasta | wc)
        echo "valid end to end scaffols with telomeres = $good_contigs"
        echo "incomplete scaffolds = $bad_contigs"

        echo "select telomeric sequences for $f"
        seqkit grep -s -R -200:-1 -r -p GTGTGTGGGTGTG canu_assembly/${name}/${name}.correctedReads.fasta.gz > telomotif/$exp/${name}.filtered.GTGTG.fastq
        seqkit grep -s -R 1:200 -r -p CACACCCACACAC canu_assembly/${name}/${name}.correctedReads.fasta.gz > telomotif/$exp/${name}.filtered.CACAC.fastq
        ##Combine the reads that are tailed
        cat telomotif/$exp/${name}.filtered.GTGTG.fastq telomotif/$exp/${name}.filtered.CACAC.fastq > telomotif/$exp/${name}.filtered.telomotif.fastq

    done

    read good_contigs g_words g_chars <<< $(grep -e '>' ${exp}.complete_contigs.fasta | wc)
    read bad_contigs g_words g_chars <<< $(grep -e '>' ${exp}.incomplete_contigs.fasta | wc)
    echo "*************** ${exp} *******************"
    echo "valid end to end scaffols with telomeres = $good_contigs"
    echo "incomplete scaffolds = $bad_contigs"
    echo "******************************************"


    echo "done merging assemblies. Moving results merged assembly results to wang/${exp}"
    mkdir -p wang/${exp}
    mv ${exp}.complete_contigs.fasta wang/${exp}/
    mv ${exp}.incomplete_contigs.fasta wang/${exp}/
    mv merged_assembly wang/${exp}/

    echo "processing telomeric sequences for $exp"
    cat telomotif/$exp/*.filtered.telomotif.fastq >  telomotif/$exp/all.filtered.telomotif.fastq
    minimap2 -t 32 -x ava-ont wang/${exp}/${exp}.incomplete_contigs.fasta telomotif/$exp/all.filtered.telomotif.fastq | gzip -1 > telomotif/$exp/$exp.filtered.telomotif.paf.gz

    /ip29/jflucier/service/externe/wellinger/program/miniasm/miniasm -1 -2 -c 1 \
    -f telomotif/$exp/all.filtered.telomotif.fastq telomotif/$exp/${exp}.filtered.telomotif.paf.gz > telomotif/$exp/${exp}.filtered.telomotif.gfa
    awk '/^S/{print ">"$2"\n"$3}' telomotif/$exp/${exp}.filtered.telomotif.gfa | fold > telomotif/$exp/${exp}.out.fa

    echo "generating tracks for hub for $exp"
    minimap2 -a -t 32 $GENOME wang/${exp}/${exp}.complete_contigs.fasta > wang/${exp}/${exp}_complete_contigs.fasta.sam
    samtools view --reference $GENOME -bS wang/${exp}/${exp}_complete_contigs.fasta.sam > wang/${exp}/${exp}_complete_contigs.fasta.bam
    bamToBed -i wang/${exp}/${exp}_complete_contigs.fasta.bam > wang/${exp}/${exp}_complete_contigs.fasta.bed
    bedSort wang/${exp}/${exp}_complete_contigs.fasta.bed wang/${exp}/${exp}_complete_contigs.fasta.sorted.bed
    bedToBigBed wang/${exp}/${exp}_complete_contigs.fasta.sorted.bed $CHROMSIZES wang/${exp}/${exp}_complete_contigs.fasta.sorted.bb


    echo "done assembly for $exp"
done


cp wang/*/*.bb hub/s288c_telo/saccer3_telo/bbi/
chmod a+r hub/s288c_telo/saccer3_telo/bbi/*
