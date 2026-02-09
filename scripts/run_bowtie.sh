#! /bin/bash -l

#$ -cwd
#$ -q all.q
#$ -tc 16
#$ -pe smp 16 
#$ -l h_rt=72:00:00
#$ -l max_mem=32G

conda activate run_bowtie2

## JOB_FILE
JOB_FILE=$1

REF=$(awk "NR==$SGE_TASK_ID" $JOB_FILE | awk -F '\t' '{print $1}')
READ_1=$(awk "NR==$SGE_TASK_ID" $JOB_FILE | awk -F '\t' '{print $2}')
READ_2=$(awk "NR==$SGE_TASK_ID" $JOB_FILE | awk -F '\t' '{print $3}')
OUT_DIR=$(awk "NR==$SGE_TASK_ID" $JOB_FILE | awk -F '\t' '{print $4}')

## outdir
outdir=${OUT_DIR}/bowtie
mkdir -p ${outdir}

exec >"$outdir/stdout.txt" 2>"$outdir/stderr.txt"

## BASE_NAME (without extension)
BASE_NAME=$(basename "$REF" | rev | cut -d'.' -f2- | rev)

## PREFIX
PREFIX=${outdir}/${BASE_NAME}

## CHECK SAMPLE_IDs MATCH
B_1=$(basename $READ_1 | sed -r 's/_R.+//g')
B_2=$(basename $READ_2 | sed -r 's/_R.+//g')

if [[ ! $B_1 == $B_2 ]]; then 
    echo "Read files do not match."
    exit 1
fi

cat <<EOF

RUNNING: $0

    JOB_FILE: $JOB_FILE
    SGE_TASK_ID: $SGE_TASK_ID
    REF: $REF
    READ_1: $READ_1
    READ_2: $READ_2
    OUT_DIR: $OUT_DIR
    outdir: $outdir
    PREFIX: $PREFIX

EOF

## FASTQC
FASTQC_OUTDIR=${outdir}/fastqc
mkdir -p ${FASTQC_OUTDIR}

fastqc --threads $NSLOTS --outdir ${FASTQC_OUTDIR} $READ_1 $READ_2

## SEQKIT
SEQKIT_STATS=${PREFIX}_SEQKIT_STATS.txt

seqkit version

seqkit stats -a $READ_1 $READ_2 > $SEQKIT_STATS

## Replace spaces with tabs
sed -i 's/ \+/\t/g' $SEQKIT_STATS

## BOWTIE2
export PERL5LIB=$(which perl)

bowtie2 --version

## BUILD BOWTIE2 INDEX
BT2_INDEX=${outdir}/${BASE_NAME}_BOWTIE2_INDEX

bowtie2-build --threads $NSLOTS $REF ${BT2_INDEX}

## OUT_BAM
OUT_BAM=${PREFIX}.bam

bowtie2 --threads $NSLOTS --quiet -x ${BT2_INDEX} -1 $READ_1 -2 $READ_2 | \
    samtools view --threads $NSLOTS -bS - | samtools collate --threads $NSLOTS -O -u - | \
    samtools fixmate --threads $NSLOTS -m -u - - | samtools sort --threads $NSLOTS -u - | \
    samtools markdup --threads $NSLOTS - ${OUT_BAM}

if [[ ! -s ${OUT_BAM} ]]; then
	echo "$0: BAM file does not exist or is empty\n"
    echo -e "\tOUT_BAM: ${OUT_BAM}\n"
	exit 1
fi

## SAMTOOLS INDEX (.BAI)
samtools index $OUT_BAM

## SAMTOOLS FA_IDX
ST_FAI_IDX=${outdir}/$(basename ${REF}).fai

samtools faidx --fai-idx ${ST_FAI_IDX} $REF

## SAMTOOLS FLAGSTAT
ST_FLAGSTAT=${PREFIX}_SAMTOOLS_FLAGSTAT.txt

samtools flagstat -O tsv $OUT_BAM > $ST_FLAGSTAT

## SAMTOOLS IDXSTATS
ST_IDXSTATS=${PREFIX}_SAMTOOLS_IDXSTATS.txt

samtools idxstats $OUT_BAM > ${ST_IDXSTATS}

## SAMTOOLS COVERAGE
ST_COVERAGE=${PREFIX}_SAMTOOLS_COVERAGE.txt

samtools coverage --depth 0 $OUT_BAM > ${ST_COVERAGE}

sed -i 's/#//g' $ST_COVERAGE

## SAMTOOLS DEPTH
ST_DEPTH=${PREFIX}_SAMTOOLS_DEPTH.txt

samtools depth --threads $NSLOTS -a $OUT_BAM | awk 'BEGIN {OFS="\t"} {$1=$1; print}' > ${ST_DEPTH}

## SAMTOOLS MEAN DEPTH
MEAN_DEPTH=$(awk '{sum+=$3} END { print sum/NR}' $ST_DEPTH)

N_SEQ=$(awk '{print $1}' $ST_DEPTH | uniq | wc -l)

TOTAL_LEN=$(awk '{print $1}' $ST_DEPTH | wc -l)

TOTAL_MAPPED=$(samtools view -c -F 4 $OUT_BAM)

TOTAL_UNMAPPED=$(samtools view -c -f 4 $OUT_BAM)

ST_MEAN_DEPTH=${PREFIX}_SAMTOOLS_MEAN_DEPTH.txt

## HEADERS
echo -e "BASE_NAME\tN_SEQ\tTOTAL_LEN\tMEAN_DEPTH\tTOTAL_MAPPED\tTOTAL_UNMAPPED\tOUT_BAM\tREF\tST_FAI_IDX\tREAD_1\tREAD_2" > ${ST_MEAN_DEPTH}

paste -d '\t' <(echo $BASE_NAME) <(echo $N_SEQ) <(echo $TOTAL_LEN) \
    <(echo $MEAN_DEPTH) <(echo $TOTAL_MAPPED) <(echo $TOTAL_UNMAPPED) \
    <(echo $OUT_BAM) <(echo $REF) <(echo $ST_FAI_IDX) <(echo $READ_1) \
    <(echo $READ_2) >> ${ST_MEAN_DEPTH}


## RUN QUAST
QUAST_OUTDIR=${outdir}/quast

mkdir -p ${QUAST_OUTDIR}

quast --output-dir ${QUAST_OUTDIR} --min-contig 400 --threads ${NSLOTS} \
    --est-ref-size 18555 --contig-thresholds 0,500,1000,1500,2500,3000,3500 \
    --bam ${OUT_BAM} ${REF}

