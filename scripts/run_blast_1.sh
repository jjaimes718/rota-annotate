#! /bin/bash -l

#$ -cwd
#$ -q all.q 
#$ -tc 16
#$ -pe smp 16
#$ -l h_rt=72:00:00
#$ -l max_mem=32G

conda activate run_blast

JOB_FILE=$1

IN_FASTA=$(awk "NR==$SGE_TASK_ID" $JOB_FILE | awk -F '\t' '{print $1}')
OUT_DIR=$(awk "NR==$SGE_TASK_ID" $JOB_FILE | awk -F '\t' '{print $4}')

## outdir
outdir=${OUT_DIR}/blast_1
mkdir -p ${outdir}

exec >"$outdir/stdout.txt" 2>"$outdir/stderr.txt"

## BLAST_DB
BLAST_DB="/scicomp/reference-pure/ncbi-blast-databases/nt"

## check file exists & is not empty
if [[ ! -s $IN_FASTA ]]; then
  echo "$0: IN_FASTA file does not exist or is empty."
  echo -e "\tIN_FASTA: $IN_FASTA\n"
  exit 1
fi

## BASE_NAME (without extension)
BASE_NAME=$(basename "$IN_FASTA" | rev | cut -d'.' -f2- | rev)

## OUT_BLAST
OUT_BLAST=${outdir}/blast_${BASE_NAME}.txt

cat <<EOF

RUNNING: $0
  JOB_FILE: $JOB_FILE
  SGE_TASK_ID: $SGE_TASK_ID
  IN_FASTA: $IN_FASTA
  OUT_DIR: $OUT_DIR
  outdir: $outdir
  OUT_BLAST: $OUT_BLAST
  BLAST_DB: $BLAST_DB
  NSLOTS: $NSLOTS
EOF

## BLAST VERSION
blastn -version

## BLAST PARAMETERS
blastn -task megablast \
	-db $BLAST_DB \
	-query $IN_FASTA \
	-out $OUT_BLAST \
	-max_target_seqs 10 \
	-evalue .001 \
	-num_threads $NSLOTS \
	-outfmt "6 qseqid sskingdoms stitle sacc pident qlen length slen mismatch gaps qstart qend sstart send sframe evalue staxids sscinames scomnames sblastnames sseqid sallseqid sgi sallgi sallacc bitscore score nident positive gapopen ppos frames qframe salltitles sstrand qcovs qcovhsp qcovus"

## COUNT
COUNT=$(awk 'NR>1 {print}' $OUT_BLAST | wc -l)

if [ "$COUNT" -gt 0 ]; then
  # remove commas from blast output
  sed -i "s/,//g" $OUT_BLAST
  
HEADERS=$(echo "qseqid sskingdoms stitle sacc pident qlen length slen mismatch gaps qstart qend sstart send sframe evalue staxids sscinames scomnames sblastnames sseqid sallseqid sgi sallgi sallacc bitscore score nident positive gapopen ppos frames qframe salltitles sstrand qcovs qcovhsp qcovus"| sed 's/ /\t/g')

sed -i "1i$HEADERS" "$OUT_BLAST"

else
  echo "$0: No BLASTN hits for $IN_FASTA" > $OUT_BLAST
  exit 1
fi

