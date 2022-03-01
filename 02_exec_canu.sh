
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

mkdir -p /home/jflucier/localhost/projet/E2EAssembler/canu_assembly
for f in /home/jflucier/Documents/service/externe/wellinger/20211116_subtelomere_analysis/data/21029_bar6_wt.*.fastq
do
    __fastq_file=$(basename $f)
    __fastq_group=${__fastq_file%.fastq}
    mkdir -p /home/jflucier/localhost/projet/E2EAssembler/canu_assembly/${__fastq_group}
    rm -fr  /home/jflucier/localhost/projet/E2EAssembler/canu_assembly/${__fastq_group}/*
    echo "running assembly on $f"
    echo "ouptutting resulting assembly in /home/jflucier/localhost/projet/E2EAssembler/canu_assembly/${__fastq_group}"
    /home/jflucier/app/canu/build/bin/canu/canu useGrid=false \
    -p ${__fastq_group} -d /home/jflucier/localhost/projet/E2EAssembler/canu_assembly/${__fastq_group} \
    genomeSize=12.21m maxThreads=4 maxMemory=8 \
    -nanopore $f
    echo "done assembly of $f"
done

echo "done running all assemblies"

