# Camelina Transcriptome Analysis Operation Manual

## 1. Manual Overview

This manual provides a detailed guide to the complete Camelina transcriptome analysis workflow based on the `fastp+hisat2+samtools+stringtie2` toolchain. It covers environment preparation, data preprocessing, sequence alignment, transcript assembly, and expression quantification. The workflow is designed to process RNA-seq data from developing Camelina seeds in the PRJNA771413 project, ultimately generating transcript/gene expression matrices to support subsequent functional analysis.

## 2. Environment and Dependency Preparation

### 2.1 Hardware Requirements

- **CPU**: 8 cores or more (the script is configured with 12 threads by default; adjust the `THREADS` parameter based on hardware specifications).
- **Memory**: At least 16GB (alignment and assembly steps require significant memory).
- **Storage**: Reserve 50GB or more (adjust based on raw data size to store raw data, intermediate files, and results).

### 2.2 Software Dependencies

Install the following software in advance and ensure they can be called normally (add to system environment variables or specify absolute paths in the script if needed):

| Software Tool    | Version Requirement         | Core Function                                         |
| ---------------- | --------------------------- | ----------------------------------------------------- |
| fastp            | v0.20.0+                    | Sequencing data quality control and adapter trimming  |
| Hisat2           | v2.2.1+                     | Genome sequence alignment                             |
| Samtools         | v1.10+                      | BAM file processing (sorting, indexing)               |
| StringTie2       | v3.0.0+                     | Transcript assembly and expression quantification     |
| Python           | v3.7+                       | Running scripts for expression extraction and merging |
| Python Libraries | pandas v1.0+, re (built-in) | Data processing and regular expression matching       |

### 2.3 Reference File Preparation

Place the following reference files in the `0ref` directory (download independently and ensure file integrity):

1. Camelina genome index files: `CS_genome` (built in advance using Hisat2; reference command: `hisat2-build Camelina_sativa_genome.fasta 0ref/CS_genome`).
2. Genome annotation file: `Camelina_sativa.Cs.57-plus1-7.gtf` (used to guide transcript assembly).

## 3. File Directory Structure

### 3.1 Initial Directory Creation

Run the following script to automatically create all directories required for the analysis (script file: `MD_TXT.sh`):

bash

```bash
# Grant execution permission to the script
chmod +x MD_TXT.sh
# Run the script to create directories
./MD_TXT.sh
```

Directories created and their functional descriptions:

| Directory Name | Functional Description                                       |
| -------------- | ------------------------------------------------------------ |
| `0ref`         | Stores reference genome index and annotation files           |
| `1raw_data`    | Stores raw FASTQ data files                                  |
| `2clean_data`  | Stores cleaned data after fastp quality control (intermediate directory; intermediate files are automatically deleted after quality control) |
| `3align_data`  | Stores sorted BAM files and indexes from Hisat2 alignment    |
| `5stringtie`   | Stores StringTie2 transcript assembly results (GTF files)    |
| `fastp`        | Stores fastp quality control reports (JSON and HTML formats) |

### 3.2 Placement of Raw Data and Sample List

1. **Raw Data**: Place paired-end FASTQ compressed files from the PRJNA771413 project (format: `${sample}_1.fastq.gz`, `${sample}_2.fastq.gz`) in the `1raw_data` directory. For example: `SRR16356205_1.fastq.gz`, `SRR16356205_2.fastq.gz`.
2. **Sample List**: Place the sample ID list file `PRJNA771413_developping_seed.txt` in the root analysis directory. The file contains the prefix of each sample (e.g., `SRR16356205`) for batch sample processing.

## 4. Core Analysis Workflow (Script: Camelina-PEceshi.sh)

### 4.1 Script Parameter Modification (Mandatory)

Open `Camelina-PEceshi.sh` with a text editor and adjust the following variables based on actual file paths (skip if reference file/tool paths match the default):

bash

```bash
# 1. Path to the sample list (default in the root directory; modify if the location changes)
SAMPLE_LIST="./PRJNA771413_developping_seed.txt"
# 2. Path to reference files (ensure index and GTF files exist)
REF_INDEX="0ref/CS_genome"  # Genome index prefix
GTF_ANNOTATION="0ref/Camelina_sativa.Cs.57-plus1-7.gtf"  # Annotation file path
# 3. Absolute path to StringTie2 (specify the full path if not added to environment variables)
# Example: /mnt/d/5.software_down/stringtie-3.0.0.Linux_x86_64/stringtie-3.0.0.Linux_x86_64/stringtie
# 4. Number of threads (adjust based on CPU core count; recommended not to exceed actual core count)
THREADS=12
```

### 4.2 Script Execution Steps

1. **Grant Execution Permission**:

   bash

   ```bash
   chmod +x Camelina-PEceshi.sh
   ```

2. **Start the Analysis Workflow**:

   bash

   ```bash
   # Run in the foreground (suitable for real-time log viewing)
   ./Camelina-PEceshi.sh
   # Or run in the background (suitable for long-term analysis; logs saved to log.txt)
   nohup ./Camelina-PEceshi.sh > log.txt 2>&1 &
   ```

### 4.3 Description of Core Workflow Steps

The script processes each sample in batches in the following order. Each step outputs logs; if a step fails, an error is prompted, and the script proceeds to the next sample:

1. **Step 1: fastp Quality Control**
   - Function: Removes low-quality bases (default Q20 threshold), trims sequencing adapters, and filters short sequences.
   - Input: `1raw_data/${sample}_1.fastq.gz`, `1raw_data/${sample}_2.fastq.gz`.
   - Output:
     - Cleaned data: `2clean_data/${sample}_good1.fastq.gz`, `2clean_data/${sample}_good2.fastq.gz` (automatically deleted after quality control to save space).
     - Quality control reports: `fastp/${sample}.json` (detailed data), `fastp/${sample}.html` (visual report; open with a browser to view).
2. **Step 2: Hisat2 Sequence Alignment and BAM Sorting**
   - Function: Aligns cleaned data to the Camelina genome and directly generates sorted BAM files (skips intermediate SAM files to improve efficiency).
   - Key Parameters:
     - `--dta`: Optimizes alignment results for downstream transcript assembly.
     - `--rna-strandness RF`: Specifies RNA strand specificity (RF mode is commonly used for paired-end sequencing; change to F/R for single-end sequencing).
     - `--sensitive`: Uses high-sensitivity alignment mode.
   - Output: `3align_data/${sample}.sorted.bam` (sorted BAM file).
3. **Step 3: Samtools BAM Index Construction**
   - Function: Creates an index for the sorted BAM file to enable fast reading by StringTie2.
   - Output: `3align_data/${sample}.sorted.bam.bai` (BAM index file).
4. **Step 4: StringTie2 Transcript Assembly and Expression Quantification**
   - Function: Assembles transcripts and calculates FPKM values for each transcript based on reference annotations and BAM alignment results.
   - Key Parameters:
     - `-G $GTF_ANNOTATION`: Guides assembly using the reference GTF.
     - `-m 200`: Filters transcripts shorter than 200bp.
   - Output: `5stringtie/${sample}.gtf` (GTF file containing transcript structure and FPKM information).

## 5. Expression Matrix Generation

### 5.1 Transcript FPKM Matrix (Script: Extra_FPKM_trans-from gtf.py)

#### 5.1.1 Script Parameter Modification

Open the script and confirm the GTF file directory path is correct (defaults to reading `./5stringtie/`; modify if the directory changes):

python行

```python
gtf_dir = "./5stringtie/"  # Directory containing GTF files (consistent with the core script output directory)
metric = "FPKM"           # Expression metric (default FPKM; no modification needed)
```

#### 5.1.2 Run the Script

bash

```bash
# Run directly with Python
python Extra_FPKM_trans-fromgtf.py
```

#### 5.1.3 Output Result

Generates the `transcript_FPKM_matrix.csv` file. Rows are transcript IDs (transcript_id), columns are sample names, and cell values are the FPKM values of the corresponding transcript in the sample (missing values are filled with 0).

### 5.2 Gene TPM Matrix (Script: merge_stringtie_tab.py)

#### 5.2.1 Prerequisite

Ensure StringTie2 outputs gene expression files in `.tab` format. If not generated by the core script, add the `-A ${STRINGTIE_DIR}/${sample}.tab` parameter to the StringTie2 command in the core script and re-run the core script.

#### 5.2.2 Script Parameter Modification

Confirm the `.tab` file directory path:

python运行

```python
input_dir = "./5stringtie"  # Directory containing .tab files output by StringTie
```

#### 5.2.3 Run the Script

bash

```bash
python merge_stringtie_tab.py
```

#### 5.2.4 Output Result

Generates the `TPM_expression_matrix.csv` file. Rows are gene IDs (Gene ID), columns are sample names, and cell values are the TPM values of the corresponding gene in the sample (duplicate genes are summed by TPM; missing values are filled with 0).