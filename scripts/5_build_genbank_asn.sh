#! /bin/bash -l

conda activate run_table2asn

## USAGE
usage() {
    echo -e "\n\t\tUsage: $0 --in_dir <PROJECT_DIRECTORY>\n"
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
                echo -e "\n$0: --in_dir requires a directory path"
                usage
            fi
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo -e "\n$0: Unknown option: $1"
            usage
            ;;
    esac
done

## Validate IN_DIR
if [[ -z "$IN_DIR" ]]; then
    echo -e "\n$0: --in_dir (PROJECT_DIRECTORY) is required."
    usage
fi

## Validate that IN_DIR exists and is a directory
if [[ ! -d "$IN_DIR" ]]; then
    echo -e "\n$0: '$IN_DIR' is not a valid directory."
    exit 1
fi

## DIR_OUTPUT
DIR_OUTPUT=${IN_DIR}/5_table2asn_output

if [[ -d ${DIR_OUTPUT} ]]; then
  #echo -e "\n$0: Clearing previous table2asn output files..."
  #echo -e "\n\tDIR_OUTPUT: $DIR_OUTPUT\n"
  rm ${DIR_OUTPUT}/*
fi

mkdir -p ${DIR_OUTPUT}

exec >"$DIR_OUTPUT/STDOUT.txt" 2>"$DIR_OUTPUT/STDERR.txt"

## PROJECT_NAME
PROJECT_NAME=$(basename $IN_DIR)

## DIR_INPUT
DIR_INPUT=${IN_DIR}/4_table2asn_input

check_sbt_input=$(find ${DIR_INPUT} -type f -name "*.sbt" | wc -l)

if [[ $check_sbt_input -ne 0 ]]; then
  rm ${DIR_INPUT}/*.sbt
fi

## IN_SBT_1 (GenBank submission template)
IN_SBT_1=$(find ${IN_DIR}/GenBank_Submission_Template -type f -name "*.sbt")

if [[ ! -s $IN_SBT_1 ]]; then
  echo -e "\n$0: Cannot find input GenBank submission template file (IN_SBT_1)\n"
  exit 1
fi

## Rename IN_SBT_1 to IN_SBT_2 in OUT_DIR
## table2asn requires all input files to be in the same directory 
## with identical base names (except for file extensions)
IN_SBT_2=${DIR_INPUT}/${PROJECT_NAME}.sbt	

cat <<EOF

RUNNING: $0
	IN_DIR: $IN_DIR
	PROJECT_NAME: $PROJECT_NAME

    DIR_INPUT: $DIR_INPUT
    DIR_OUTPUT: $DIR_OUTPUT
    IN_SBT_1: $IN_SBT_1
    IN_SBT_2: $IN_SBT_2

EOF

cp -p ${IN_SBT_1} ${IN_SBT_2}

## LOG_FILE
LOG_FILE=${DIR_OUTPUT}/LOG.txt

## TABLE2ASN VERSION
echo -e "\nTABLE2ASN VERSION:\n"

table2asn -version-full

## TABLE2ASN
table2asn -t ${IN_SBT_2} \
  -indir ${DIR_INPUT} \
  -outdir ${DIR_OUTPUT} \
  -a s \
  -V vb \
  -j "[moltype=genomic] [molecule=rna]" \
  -T \
  -logfile ${LOG_FILE} \
  -H y

## options
# -H <String> : Hold Until Publish 
#     y : Hold for one year
#     mm/dd/yyyy

# -linkage-evidence-file

