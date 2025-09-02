#!/bin/bash
# RNA-seq data analysis pipeline script
# Optimized version: Directly generate sorted BAM files from HiSAT2, avoiding intermediate SAM files

# Define variables
SAMPLE_LIST="./PRJNA771413_developping_seed.txt"
FASTQ_DIR="./1raw_data"
CLEAN_DIR="./2clean_data"
BAM_DIR="./3align_data"
FASTP_DIR="./fastp"
STRINGTIE_DIR="./5stringtie"
THREADS=12
REF_INDEX="0ref/CS_genome"
GTF_ANNOTATION="0ref/Camelina_sativa.Cs.57-plus1-7.gtf"

# Create necessary directories
mkdir -p $CLEAN_DIR $BAM_DIR $FASTP_DIR $STRINGTIE_DIR

# Check if input files exist
if [ ! -f "$SAMPLE_LIST" ]; then
    echo "Error: Sample list file $SAMPLE_LIST does not exist!"
    exit 1
fi

# Check if GTF annotation file exists
if [ ! -f "$GTF_ANNOTATION" ]; then
    echo "Error: GTF annotation file $GTF_ANNOTATION does not exist!"
    exit 1
fi

# Process each sample in a loop
while IFS= read -r sample; do
    echo "Starting processing sample: $sample"
    
    # Step 1: Quality control and trimming using fastp
    echo "  Performing quality control..."
    fastp -i ${FASTQ_DIR}/${sample}_1.fastq.gz -I ${FASTQ_DIR}/${sample}_2.fastq.gz \
          -o ${CLEAN_DIR}/${sample}_good1.fastq.gz -O ${CLEAN_DIR}/${sample}_good2.fastq.gz \
          -j ${FASTP_DIR}/${sample}.json -h ${FASTP_DIR}/${sample}.html || {
        echo "Error: fastp processing failed for sample $sample!"
        continue
    }
    
    # Step 2: Sequence alignment using hisat2 and generate sorted BAM file directly
    echo "  Performing sequence alignment and generating sorted BAM file..."
    hisat2 --dta --rna-strandness RF -p $THREADS -x $REF_INDEX --sensitive \
           -1 ${CLEAN_DIR}/${sample}_good1.fastq.gz -2 ${CLEAN_DIR}/${sample}_good2.fastq.gz | \
    samtools view -@$THREADS -b - | \
    samtools sort -@$THREADS -o ${BAM_DIR}/${sample}.sorted.bam || {
        echo "Error: hisat2 alignment or sorting failed for sample $sample!"
        continue
    }
    
    # Delete intermediate files
    rm -f ${CLEAN_DIR}/${sample}_good*.fastq.gz
    
    # Step 3: Prepare index for StringTie
    echo "  Preparing BAM index for StringTie..."
    samtools index -@$THREADS ${BAM_DIR}/${sample}.sorted.bam || {
        echo "Error: samtools indexing failed for sample $sample!"
        continue
    }
    
    # Step 4: Transcript assembly and expression quantification using StringTie
    echo "  Performing transcript assembly with StringTie..."
    /mnt/d/5.software_down/stringtie-3.0.0.Linux_x86_64/stringtie-3.0.0.Linux_x86_64/stringtie ${BAM_DIR}/${sample}.sorted.bam \
              -o ${STRINGTIE_DIR}/${sample}.gtf \
              -p $THREADS \
              -G $GTF_ANNOTATION \
              -m 200 \
              || {
        echo "Error: StringTie processing failed for sample $sample!"
        continue
    }
    
    echo "Sample $sample processing completed!"
    echo "------------------------"
    
done < "$SAMPLE_LIST"

echo "All samples processing completed!"
