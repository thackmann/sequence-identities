#!/bin/bash

if [ x$usearch == x ] ; then
	echo Must set \$usearch >> /dev/stderr
	exit 1
fi

rm -rf ../out
mkdir -p ../out
cd ../out

#Merge paired reads (minmergelen and maxmergelen set to encompass length range of reference sequences +/-5 bp) 
$usearch -fastq_mergepairs ../data/Mock1_S1_L001_R1_001.fastq -reverse ../data/Mock1_S1_L001_R2_001.fastq -fastq_minmergelen 400 -fastq_maxmergelen 436  -fastq_maxdiffs 500 -fastqout ../out/Mock1_S1_L001_Merged_001_V34.fq #For V34
$usearch -fastq_mergepairs ../data/Mock1_S1_L001_R1_001.fastq -reverse ../data/Mock1_S1_L001_R2_001.fastq -fastq_minmergelen 248 -fastq_maxmergelen 259  -fastq_maxdiffs 500 -fastqout ../out/Mock1_S1_L001_Merged_001_V4.fq #For V4
$usearch -fastq_mergepairs ../data/Mock1_S1_L001_R1_001.fastq -reverse ../data/Mock1_S1_L001_R2_001.fastq -fastq_minmergelen 368 -fastq_maxmergelen 383  -fastq_maxdiffs 500 -fastqout ../out/Mock1_S1_L001_Merged_001_V45.fq #For V34

#Annotate reads with UPARSE-REF
$usearch -uparse_ref ../out/Mock1_S1_L001_Merged_001_V34.fq -db ../data/reference_sequences_V34_no_primers.fasta -strand plus -uparseout Mock1_S1_L001_Merged_001_V34.up #For V34
$usearch -uparse_ref ../out/Mock1_S1_L001_Merged_001_V4.fq -db ../data/reference_sequences_V4_no_primers.fasta -strand plus -uparseout Mock1_S1_L001_Merged_001_V4.up #For V4
$usearch -uparse_ref ../out/Mock1_S1_L001_Merged_001_V45.fq -db ../data/reference_sequences_V45_no_primers.fasta -strand plus -uparseout Mock1_S1_L001_Merged_001_V45.up #For V45
