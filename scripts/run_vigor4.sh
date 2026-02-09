#! /bin/bash -l

#$ -cwd
#$ -q all.q 
#$ -tc 16
#$ -pe smp 16
#$ -l h_rt=72:00:00
#$ -l max_mem=32G

conda activate run_vigor

## JOB_FILE
JOB_FILE=$1

## IN_FASTA
IN_FASTA=$(awk "NR==$SGE_TASK_ID" $JOB_FILE | awk -F '\t' '{print $1}')

## OUT_DIR
OUT_DIR=$(awk "NR==$SGE_TASK_ID" $JOB_FILE | awk -F '\t' '{print $4}')

## outdir
outdir=${OUT_DIR}/vigor
mkdir -p ${outdir}

exec >"$outdir/stdout.txt" 2>"$outdir/stderr.txt"

if [[ ! -s $JOB_FILE ]]; then
  echo "$0: JOB_FILE does not exist OR is an empty file."
  echo -e "\tUsage: $0 JOB_FILE\n\n"
  echo -e "JOB_FILE should be a tab-delimited file with the following columns:\n\n"
  echo -e "IN_FASTA\tREAD_1\tREAD_2\tOUT_DIR\n\n"
  exit 1
fi

if [[ ! -s $IN_FASTA ]]; then
  echo "$0: IN_FASTA file does not exist"
  echo -e "\tIN_FASTA: $IN_FASTA\n"
  exit 1
fi

## BASE_NAME (without extension)
BASE_NAME=$(basename "$IN_FASTA" | rev | cut -d'.' -f2- | rev)

## PREFIX
PREFIX=${outdir}/vigor4_${BASE_NAME}

## VIGOR_DB
VIGOR_DB=rtva_db

cat <<EOF

RUNNING: $0
  JOB_FILE: $JOB_FILE
  SGE_TASK_ID: $SGE_TASK_ID
  IN_FASTA: $IN_FASTA
  OUT_DIR: $OUT_DIR
  outdir: $outdir
  PREFIX: $PREFIX
  VIGOR_DB: $VIGOR_DB
  NSLOTS: $NSLOTS
EOF

## RUN VIGOR4
vigor4 --version

vigor4 --input-fasta $IN_FASTA --output-prefix ${PREFIX} --reference-database $VIGOR_DB

## PARSE VIGOR4 OUTPUT
./scripts/parse_vigor4_output.R --in_dir ${outdir}

conda deactivate
