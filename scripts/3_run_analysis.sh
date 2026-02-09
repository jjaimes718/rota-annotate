#! /bin/bash -l

## USAGE
usage() {
    echo -e "\t\tUsage: $0 --in_dir <PROJECT_DIRECTORY>\n"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --in_dir)
            if [[ -n "${2:-}" && ! "$2" =~ ^-- ]]; then
                IN_DIR="$2"
                shift 2
            else
                echo "$0: --in_dir requires a directory path."
                usage
            fi
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

## Validate IN_DIR
if [[ -z "$IN_DIR" ]]; then
    echo "$0: --in_dir is required."
    usage
fi

if [[ ! -d "$IN_DIR" ]]; then
    echo "$0: '$IN_DIR' is not a valid directory."
    exit 1
fi

## PROJECT_NAME
PROJECT_NAME=$(basename $IN_DIR)

## JOB_FILE
JOB_FILE=$(find ${IN_DIR} -type f -name "job_file*.txt" | grep -v "bowtie")

if [[ ! -s $JOB_FILE ]]; then
  echo -e "$0: JOB_FILE does not exist or is empty\n"
  echo -e "\tJOB_FILE: $JOB_FILE\n"
  echo -e "JOB_FILE should be a tab-delimited file with the following columns:\n\n"
  echo -e "IN_FASTA\tREAD_1\tREAD_2\tOUT_DIR\n\n"
  exit 1
fi

## N_THREADS
N_THREADS=$(wc -l $JOB_FILE | awk '{print $1}')

## ANALYSIS_DIR
ANALYSIS_DIR=$(dirname $JOB_FILE)

LOG_DIR=${ANALYSIS_DIR}/logs
mkdir -p $LOG_DIR

cat <<EOF

RUNNING: $0
	IN_DIR: $IN_DIR
	PROJECT_NAME: $PROJECT_NAME
    JOB_FILE: $JOB_FILE
	N_THREADS: $N_THREADS
	ANALYSIS_DIR: $ANALYSIS_DIR
	LOG_DIR: $LOG_DIR
EOF

## VIGOR4
N_VIGOR=run_vigor_${PROJECT_NAME}

qsub -N ${N_VIGOR} \
	-t 1-$N_THREADS \
	-o ${LOG_DIR} \
	-e ${LOG_DIR} \
	./scripts/run_vigor4.sh $JOB_FILE

## BLAST_1
N_BLAST_1=run_blast_1_${PROJECT_NAME}

qsub -N ${N_BLAST_1} \
	-t 1-$N_THREADS \
	-o ${LOG_DIR} \
	-e ${LOG_DIR} \
	./scripts/run_blast_1.sh $JOB_FILE

## BLAST_2
N_BLAST_2=run_blast_2_${PROJECT_NAME}

qsub -N ${N_BLAST_2} \
	-t 1-$N_THREADS \
	-o ${LOG_DIR} \
	-e ${LOG_DIR} \
	./scripts/run_blast_2.sh $JOB_FILE

## BOWTIE
## Clear old job_file_bowtie 
B_NAME=$(basename $JOB_FILE | sed 's/job_file/job_file_bowtie/')

JOB_FILE_BOWTIE=${ANALYSIS_DIR}/${B_NAME}

if [[ -s $JOB_FILE_BOWTIE ]]; then
	echo -e "$0: Clearing old JOB_FILE_BOWTIE\n"
	echo -e "\tJOB_FILE_BOWTIE: $JOB_FILE_BOWTIE"
	rm $JOB_FILE_BOWTIE
fi

## Verify 'READ' files were listed in the original JOB_FILE
awk -F '\t' '$2 != "NA" && $3 != "NA" {print $0}' $JOB_FILE > ${JOB_FILE_BOWTIE}

## N_THREADS_BOWTIE
N_THREADS_BOWTIE=$(wc -l $JOB_FILE_BOWTIE | awk '{print $1}')

if [ "$N_THREADS_BOWTIE" -gt 0 ]; then
	## BOWTIE JOB NAME
	N_BOWTIE=run_bowtie_${PROJECT_NAME}

	qsub -N ${N_BOWTIE} \
		-t 1-$N_THREADS_BOWTIE \
		-o ${LOG_DIR} \
		-e ${LOG_DIR} \
		./scripts/run_bowtie.sh $JOB_FILE_BOWTIE
else
	echo -e "$0: No 'READ' files provided, therefore read coverage statistics \
			will not be reported."
	echo -e "\tN_THREADS_BOWTIE = ${N_THREADS_BOWTIE}\n"
fi
