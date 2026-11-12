#!/bin/bash

FASTQ_DIR=$1
FINAL_QC_OUT=$2
THREADS=24

mkdir -p ${FINAL_QC_OUT}/fastqc_out
mkdir -p ${FINAL_QC_OUT}/fastp_out
mkdir -p ${FINAL_QC_OUT}/multiqc_out

FASTQC_OUT="${FINAL_QC_OUT}/fastqc_out"
FASTP_OUT="${FINAL_QC_OUT}/fastp_out"
MULTIQC_OUT="${FINAL_QC_OUT}/multiqc_out"
FILENAMES="${FINAL_QC_OUT}/sample_list"

ls ${FASTQ_DIR}/*_1.fastq.gz | sed 's/_1.fastq.gz//; s|.*/||' | sort | uniq > ${FILENAMES}
cat ${FILENAMES}

module load fastqc/0.12.1
fastqc -t ${THREADS} ${FASTQ_DIR}/*.fastq.gz -o ${FASTQC_OUT}

module load multiqc/v1.31
multiqc ${FASTQC_OUT} -o ${MULTIQC_OUT}

module load fastp/v1.0.1

while read SAMPLE; do
  echo "Processing ${SAMPLE}"

  fastp -w ${THREADS} \
    --in1 ${FASTQ_DIR}/${SAMPLE}_1.fastq.gz \
    --in2 ${FASTQ_DIR}/${SAMPLE}_2.fastq.gz \
    --out1 ${FASTP_OUT}/${SAMPLE}_1_trim.fq.gz \
    --out2 ${FASTP_OUT}/${SAMPLE}_2_trim.fq.gz \
    --unpaired1 ${FASTP_OUT}/${SAMPLE}_1_se.fastq.gz \
    --unpaired2 ${FASTP_OUT}/${SAMPLE}_2_se.fastq.gz \
    --low_complexity_filter \
    -l 36 \
    --average_qual 20 \
    --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    > ${FASTP_OUT}/${SAMPLE}_fastp.log 2>&1
done < ${FILENAMES}

multiqc ${FINAL_QC_OUT} -o ${FINAL_QC_OUT}/multiqc_combined
