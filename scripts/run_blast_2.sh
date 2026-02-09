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
outdir=${OUT_DIR}/blast_2
mkdir -p ${outdir}

exec >"$outdir/stdout.txt" 2>"$outdir/stderr.txt"

## BLAST_DB
BLAST_DB="./db/rv_virus_variation/rv_virus_variation.fasta"

## TYPE (nucl or prot)
TYPE=nucl

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
  TYPE: $TYPE
  NSLOTS: $NSLOTS
EOF

if [[ $TYPE == nucl ]]; then
  blastn -version
  blastn -task megablast \
    -db $BLAST_DB \
    -query $IN_FASTA \
    -out ${outdir}/tmp.txt \
    -max_target_seqs 20 \
    -evalue .001 \
	  -num_threads $NSLOTS \
    -outfmt "6 qseqid pident length mismatch gaps qstart qend sstart send sframe evalue stitle"
elif [[ $TYPE == prot ]]; then
  blastp -version
  blastp \
    -db $BLAST_DB \
    -query $IN_FASTA \
    -out ${outdir}/tmp.txt \
    -max_target_seqs 20 \
    -evalue .001 \
    -outfmt "6 qseqid pident length mismatch gaps qstart qend sstart send sframe evalue stitle"
else
  echo -e "$0: invalid type argument\n"
  echo -e "\tTYPE: $TYPE"
  exit 1
fi

## COUNT
COUNT=$(awk 'NR>1 {print}' ${outdir}/tmp.txt | wc -l)

if [ "$COUNT" -gt 0 ]; then
  ## split blast hit definition column
  paste -d '\t' <(cut -f 1-11 ${outdir}/tmp.txt) \
    <(cut -f 12 ${outdir}/tmp.txt | sed 's/,//g; s/|/\t/g') > $OUT_BLAST
  
  rm ${outdir}/tmp.txt

  ## HEADERS
  HEADERS=$(echo "qseqid pident length mismatch gaps qstart qend sstart send \
    sframe evalue sacc host country date segment genotype def source \
    gi"|sed 's/ /\t/g')
  
  ## PREPEND HEADERS TO BLAST OUTPUT
  sed -i "1i $HEADERS" $OUT_BLAST

else
  echo "$0: No BLAST hits for $IN_FASTA"
  
  ## ADD HEADERS TO EMPTY BLAST OUTPUT
  HEADERS=$(echo "qseqid pident length mismatch gaps qstart qend sstart send \
    sframe evalue sacc host country date segment genotype def source \
    gi"|sed 's/ /\t/g')
  
  echo $HEADERS > $OUT_BLAST

  rm ${outdir}/tmp.txt
  exit 1
fi
