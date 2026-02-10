# rota-annotate

Analysis and annotation of rotavirus genome assemblies.

## Usage

### Input

```text
./data/PROJECT_NAME/1_original_input/
└── INSTRUMENT_PLATFORM_ID
    └── NGS_RUN_ID
        └── NGS_METHOD_ID
            ├── assembly/
            └── reads/
```

Input directory structure:

- Create a project directory ( PROJECT_NAME ):
- Create a sub directory ( "1_original_input" ) containing original assembly and read files
  - assembly (FASTA)
  - reads (FASTQ)

At the moment, the code is executed in 4 dfferent steps.

### 1. Generate Sample Info Table

```bash

./scripts/1_generate_sample_info_table.R --in_dir ./data/PROJECT_NAME

```

- Reads input assemblies and read files
- Extracts NGS run/method/sample IDs
- Generates a "Sample_Info" table (.xlsx) which should be completed before continuing.
  - If the table is not completed, default values are set.

### 2. Setup Analysis Input/Output

```bash

./scripts/2_analysis_input.R --in_dir ./data/PROJECT_NAME

```

- Reads sample info table and merges with input assembly info
- Generates unique assembly and contig IDs
- Cleans filenames and sequence headers
- Cleans and filters input sequences (trims Ns, length >= 400)
- Generates input/output directories (for subsequent analysis step)
- Generates analysis input file list (tab-delimited "job_file"):

```text

IN_ASSEMBY  READ_1  READ_2  OUT_DIR
PROJECT_NAME/2_analysis_input/1/ASSEMBLY_1.fas  PROJECT_NAME/2_analysis_input/1/ASSEMBLY_1_R1.fastq  PROJECT_NAME/2_analysis_input/1/ASSEMBLY_1_R2.fastq	PROJECT_NAME/2_analysis_input/1
PROJECT_NAME/2_analysis_input/2/ASSEMBLY_2.fas  PROJECT_NAME/2_analysis_input/2/ASSEMBLY_2_R1.fastq  PROJECT_NAME/2_analysis_input/2/ASSEMBLY_2_R2.fastq	PROJECT_NAME/2_analysis_input/2
PROJECT_NAME/2_analysis_input/3/ASSEMBLY_3.fas  PROJECT_NAME/2_analysis_input/3/ASSEMBLY_3_R1.fastq  PROJECT_NAME/2_analysis_input/3/ASSEMBLY_3_R2.fastq	PROJECT_NAME/2_analysis_input/3

```

### 3. Run Sequence Analysis

```bash

./scripts/3_run_analysis.sh --in_dir ./data/PROJECT_NAME

```

- Submits 4 job arrays:
  - **run_vigor4.sh**: [VIGOR4](<https://github.com/JCVenterInstitute/VIGOR4>) is a tool used for rotavirus gene prediction
  - **run_blast_1.sh**: runs blastn against NCBI's "nt" database
  - **run_blast_2.sh**: runs blastn against a custom database (contains rotavirus sequences with known genotypes)
  - **run_bowtie.sh**: runs [bowtie2](<https://github.com/BenLangmead/bowtie2>) (read coverage statistics)

### 4. Generate Analysis Summary Report

```bash

./scripts/generate_analysis_summary.R --in_dir ./data/PROJECT_NAME

```

Once all jobs/tasks are completed. This step concatenates all results and generates a final Excel summary report.
