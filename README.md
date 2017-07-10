# cintia_snp_call_pipeline

## Requirements
Rake, FASTQC, trimmomatic, bowtie2, samtools

## Usage:
rake -f Rakefile R1=READ1.fastq R2=READ2.fastq reference=Refseq.fasta samplename=SAMPLENAME sampleid=SAMPLEID :default

The command will run FASTQC, quality trimming with trimmomatic, reads mapping with bowtie2 and SAM - BAM convertion with samtools.

Repeat the command above for every sample FASTQ files.

rake -f Rakefile reference=Refseq.fasta sampleid=SAMPLEID samtools:merge

This command will merge the BAM files generated for samples with same sampleid, sort the BAM file and index it.

