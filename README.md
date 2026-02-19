
# rota-annotate

Analysis and annotation of rotavirus genome assemblies.

---

## Table of Contents

1. [About](#about)
2. [General Workflow](#general-workflow)
3. [Input](#input)

---

## About

Since forming back in 2008, the [Rotavirus Classification Working Group] (RCWG) has developed guidelines for:

   1.  Genotype classification group A rotaviruses ([Matthijnssens et al., 2008])
   2.  Standardized rotavirus strain nomenclature ([Matthijnssens et al., 2011])

Additionally, these guidelines established a standardized format for documenting associated strain data. Although their adoption has enhanced uniformity in the description and reporting of rotavirus strains, implementation in formatting of sequence records remains inconsistent. This often complicates querying and retrieval of data from sequence databases.

Rota-annotate aims to automate the implementation of RCWG guidelines, offering tools for the analysis of genome assembly sequences and the generation of finalized sequence records suitable for publication (GenBank).

[Rotavirus Classification Working Group]: <https://rega.kuleuven.be/cev/viralmetagenomics/virus-classification/rcwg>

[Matthijnssens et al., 2008]: https://pubmed.ncbi.nlm.nih.gov/18604469/

[Matthijnssens et al., 2011]: <https://pmc.ncbi.nlm.nih.gov/articles/PMC3398998/#abstract1>


## General Workflow

![Rota-annotate workflow diagram](assets/images/Work_Flow.png)

___

## Input

### A. Input Files

   1. Assembly contig sequence files (FASTA)
   2. Optional sequencing read files (FASTQ)

### B. Input directory structure

```text
./data/PROJECT_NAME/
├── 1_original_input/
│   └── INSTRUMENT_PLATFORM_ID
│       └── NGS_RUN_ID
│           └── NGS_METHOD_ID
│               ├── assembly/
│               └── reads/
└── GenBank_Submission_Template/
```

- Create a project directory (**PROJECT_NAME**)
- Create a 2 sub directories named: 
  - **1_original_input** - containing original **assembly** files (FASTA) and **read** files (FASTQ)
  - **GenBank_Submission_Template** - containing completed submission template

## Usage

___

At the moment, the code is executed in 5 dfferent steps.

### 1. Generate Sample Info Table

```bash

./scripts/1_generate_sample_info_table.R --in_dir ./data/PROJECT_NAME

```

- Reads input assemblies and read files
- Extracts NGS run/method/sample IDs
- Generates a "**Sample_Info**" table (.xlsx) which should be completed before continuing.
  - Contains sample related information required for generating sequence annotation
  - If the table is not completed, default values are set

- Sample Information Table - Column Descriptions

| Column Name | Description |
| :----------- | :----------------------- |
| **collection.date** | The date the speciman was collected (equivalent to GenBank source modifier [collection_date](<https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html#modifiers>)) |
| **country** | Geographic location of the collected sample (equivalent to GenBank source modifier [geo_loc_name](<https://www.ncbi.nlm.nih.gov/genbank/collab/country/>)) |
| **iso_alpha3_code** | Three-letter country code ([ISO 3166-1 alpha-3](<https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3>)) |
| **host** | Host species using binomial nomenclature (*Genus + species*) <br><br> Ex: Homo sapiens |
| **host.common** | Refers to species of origin common name (first component) as described by [Matthijnssens et al. 2011](<https://pmc.ncbi.nlm.nih.gov/articles/PMC3398998/>) <br><br> {*RV group*} / {***species of origin***} / {*country of identification*} / {*common name*} / {*year of identification*} / {*G-type*}{*P-type*}<br><br> Ex: RVA/**Human**-wt/USA/OM46/1998/G9P[8] |
| **type** | Refers to species of origin sample type (second component) as described by [Matthijnssens et al. 2011](<https://pmc.ncbi.nlm.nih.gov/articles/PMC3398998/>) <br><br> {*RV group*} / {***species of origin***} / {*country of identification*} / {*common name*} / {*year of identification*} / {*G-type*}{*P-type*} <br><br> Ex: RVA/Human-**wt**/USA/OM46/1998/G9P[8] <br><br> Where *"**wt**"=wild type, "**tc**"=tissue-culture adapted, "**lab**"=lab-genereated or lab-engineered, "**X**"=unknown* |
| **common_name.prefix** | Optional prefix for "common name" component of proposed rotavirus strain nomenclature as described by [Matthijnssens et al. 2011](<https://pmc.ncbi.nlm.nih.gov/articles/PMC3398998/#S3>): <br><br> {*RV group*} / {*species of origin*} / {*country of identification*} / {***common name***} / {*year of identification*} / {*G-type*}{*P-type*} <br><br> Ex: RVA/Human-wt/HUN/**BP**1062/2004/G8P[14] <br><br> Where ***common_name.prefix***="BP" results in ***common name***="**BP**1062" |
| **isolation.source** | Describes the physical, environmental and/or local geographical source of the biological sample from which the sequence was derived (equivalent to GenBank source modifier [Isolation_source](<https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html>)) |


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
PROJECT_NAME/2_analysis_input/1/ASSEMBLY_1.fas  PROJECT_NAME/2_analysis_input/1/ASSEMBLY_1_R1.fastq  PROJECT_NAME/2_analysis_input/1/ASSEMBLY_1_R2.fastq PROJECT_NAME/2_analysis_input/1
PROJECT_NAME/2_analysis_input/2/ASSEMBLY_2.fas  PROJECT_NAME/2_analysis_input/2/ASSEMBLY_2_R1.fastq  PROJECT_NAME/2_analysis_input/2/ASSEMBLY_2_R2.fastq PROJECT_NAME/2_analysis_input/2
PROJECT_NAME/2_analysis_input/3/ASSEMBLY_3.fas  PROJECT_NAME/2_analysis_input/3/ASSEMBLY_3_R1.fastq  PROJECT_NAME/2_analysis_input/3/ASSEMBLY_3_R2.fastq PROJECT_NAME/2_analysis_input/3

```

### 3. Run Sequence Analysis

```bash

./scripts/3_run_analysis.sh --in_dir ./data/PROJECT_NAME

```

- Submits 4 job arrays:
  - **run_vigor4.sh**: [VIGOR4] is a tool used for rotavirus gene prediction
  - **run_blast_1.sh**: runs [blastn] against NCBI's [nt] database
  - **run_blast_2.sh**: runs [blastn] against a custom database (contains rotavirus sequences with known genotypes)
  - **run_bowtie.sh**: runs [bowtie2] (read coverage statistics)

[VIGOR4]: <https://github.com/JCVenterInstitute/VIGOR4>
[bowtie2]: <https://github.com/BenLangmead/bowtie2>
[blastn]: <https://www.ncbi.nlm.nih.gov/books/NBK569856/>
[nt]: <https://ftp.ncbi.nlm.nih.gov/blast/db/README>


### 4. Generate Analysis Summary Report

```bash

./scripts/4_generate_analysis_summary.R --in_dir ./data/PROJECT_NAME

```

Once all jobs/tasks are completed. This step concatenates all results and generates a final Excel summary report.

### 5. Generate GenBank Submission File

```bash

./scripts/5_build_genbank_asn.sh --in_dir ./data/PROJECT_NAME

```

- [table2asn]\: tool used to generate sequence records for submission to GenBank

[table2asn]: <https://www.ncbi.nlm.nih.gov/genbank/table2asn/?utm_source=ncbi_insights&utm_medium=referral&utm_campaign=table2asn-updated-20230706>


```text

./data/PROJECT_NAME/
├── 1_original_input
├── 2_analysis_input
├── 3_analysis_output
├── 4_table2asn_input/
│   ├── PROJECT_NAME.fsa
│   ├── PROJECT_NAME.src
│   └── PROJECT_NAME.tbl
├── 5_table2asn_output/
│   ├── PROJECT_NAME.gbf
│   └── PROJECT_NAME.sqn
└── GenBank_Submission_Template/
    └── GenBank_Submission_Template.sbt
    
```

- **Table2asn Input:**
  - [Nucleotide sequence file] (***.fsa**)
  - [Feature Table] (***.tbl**)
  - [Source Table] (***.src**)
  - [GenBank Submission Template] (***.sbt**)

[Nucleotide sequence file]: <https://www.ncbi.nlm.nih.gov/genbank/table2asn/?utm_source=ncbi_insights&utm_medium=referral&utm_campaign=table2asn-updated-20230706#fsa>

[Feature Table]: <https://www.ncbi.nlm.nih.gov/genbank/table2asn/?utm_source=ncbi_insights&utm_medium=referral&utm_campaign=table2asn-updated-20230706#tbl>

[Source Table]: <https://www.ncbi.nlm.nih.gov/genbank/table2asn/?utm_source=ncbi_insights&utm_medium=referral&utm_campaign=table2asn-updated-20230706#src>

[GenBank Submission Template]: <https://www.ncbi.nlm.nih.gov/genbank/table2asn/?utm_source=ncbi_insights&utm_medium=referral&utm_campaign=table2asn-updated-20230706#Template>


- **Table2asn Output:**
  - GenBank flatfile (***.gbf**)
  - Abstract Syntax Notation 1 (ASN.1) text file (***.sqn**)
