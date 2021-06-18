#!/bin/bash

# modified from https://umi-tools.readthedocs.io/en/latest/QUICK_START.html#
# dependencies
# - STAR aligner
# - samtools
# - fsubread package
# - UMI tools
# - genome and annotation files
# - cutadapt
# - bbtools (optionally)

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# create genome index for STAR, replace paths with actual paths to fasta and gtf files
# change the overhang to read length -1 
STAR \
--runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /shared/genome/ncbi_grcm38/star \
--genomeFastaFiles /shared/genome/ncbi_grcm38/fasta/ncbi_grcm38.fa \
--sjdbGTFfile /shared/genome/ncbi_grcm38/annotation/ncbi_grcm38.gtf \
--sjdbOverhang 49;

# this part is done if there are too many fastq files (I had 960). You pool all the fastq files into 
# one so it makes the processing easier, it introduces a specific artificial barcode into the read name
# rest of the pipeline assumes 1 merged fastq file, but it can be updated to read individual ones
printf "%s\t%s\n" "file" "barcode" >> barcode_file.txt

# initial value gives 5 bases
A=256
for f in *R2.fastq.gz
do
file=${f%.R2.fastq.gz}
f2=${f%R2.fastq.gz}R1.fastq.gz
barcode=$(echo "obase=4;$A" | bc | awk '{gsub(/0/, "G"); gsub(/1/, "A"); gsub(/2/, "T"); gsub(/3/, "C"); print}')

awk -v bar=$barcode '{ if (NR%4==1) {print $1"#"bar } else {print} }' <(gzip -cdk $f) >> R2_reads_mod.fastq
awk -v bar=$barcode '{ if (NR%4==1) {print $1"#"bar } else {print} }' <(gzip -cdk $f2) >> R1_reads_mod.fastq
printf "%s\t%s\n" $file $barcode >> barcode_file.txt
let "A=$A+1" 
done

# compress to save space
pigz --best -p 16 R1_reads_mod.fastq;
pigz --best -p 16 R2_reads_mod.fastq;    

# filtering all the barcode files that don't have specified length (this case 8 for UMI + 10 for RT barcode)
cutadapt -j 16 -m 18: -o R1_reads_mod_fil.fastq.gz -p R2_reads_mod_fil.fastq.gz R1_reads_mod.fastq.gz R2_reads_mod.fastq.gz;

# create a RT barcode whitelist from R1 reads
# - change the pattern for the actual one
# - change the cell number to actual amount of RT barcodes
umi_tools whitelist --stdin R1_reads_mod_fil.fastq.gz \
                    --bc-pattern=NNNNNNNNCCCCCCCCCC \
                    --set-cell-number=384 \
                    --subset-reads=300000000 \
                    --error-correct-threshold 3 \
					--log2stderr > whitelist.txt;
         
            
# PUT UMI and barcodes into read names based on generated whitelist
# alternatively use whitelist_true.txt which contains the actual barcodes
umi_tools extract --bc-pattern=NNNNNNNNCCCCCCCCCC \
                  --stdin R1_reads_mod_fil.fastq.gz \
                  --stdout R1_reads_extracted.fastq.gz \
                  --read2-in R2_reads_mod_fil.fastq.gz \
                  --read2-out=R2_reads_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist=whitelist_true.txt;

# this changes the names of fastq reads to concatenate the artificial into the RT read barcode
# that will be 'demultiplexed' later in R during analysis
awk '{ if (NR%4==1) {gsub(/#/, ":"); gsub(/_/, ":"); print $0 } else {print}}' <(gzip -cdk R2_reads_extracted.fastq.gz) > R2_reads_extracted_1.fastq
awk -F: '{ if (NR%4==1) {print $1":"$2":"$3":"$4":"$5":"$6":"$7"_"$8$9"_"$10 } else {print} }' R2_reads_extracted_1.fastq > R2_reads_extracted_2.fastq

# remove temporarily generated files that will no longer be needed
rm -f R1_reads_extracted.fastq.gz;
rm -f R2_reads_extracted.fastq.gz;
rm -f R1_reads_mod_fil.fastq.gz;
rm -f R2_reads_mod_fil.fastq.gz;
rm -f R2_reads_extracted_1.fastq;

# split mouse and human reads to see the proportions
# in the end I aligned all the reads, since there were too few mouse specific ones
bbsplit.sh -Xmx60g minid=0.76 maxindel=16k minhits=1 ambig2=split in=R2_reads_extracted_2.fastq ref=/shared/genome/GRCm38/fasta/GRCm38.fa,/shared/genome/GRCh38/fasta/GRCh38.fa basename=out_%.fastq
             
# compress to save space
pigz --best -p 16  R2_reads_extracted_2.fastq;

# align reads to genome
# parameters are modified according to:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5792058/bin/NIHMS902588-supplement-Supplemental.pdf    
STAR --runThreadN 16 \
     --genomeDir /shared/genome/mm10/star \
     --readFilesIn R2_reads_extracted_2.fastq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFilterMismatchNmax 33 \
     --alignEndsType Local \
     --alignSJoverhangMin 8 \
     --outFilterScoreMinOverLread 0 \
     --outFilterMatchNminOverLread 0 \
     --outFilterMatchNmin 0 \
     --seedSearchStartLmax 12 \
     --outFilterType BySJout \
     --sjdbOverhang 50 \
     --quantMode TranscriptomeSAM \
     --limitBAMsortRAM 20000000000;
 
 
samtools sort -@ 16 -m 2G Aligned.sortedByCoord.out.bam -o Aligned.sortedByCoord.out2.bam;

# count features with subread package
# output as BAM file for the UMI tools deduplication
# in this case it would be interesting to substitute for mmquant
# due to relaxed aligning parameters, many reads are multimapped
featureCounts -a /shared/genome/mm10/annotation/mm10.gtf \
			  -o gene_assigned.txt \
			  -R BAM Aligned.sortedByCoord.out.bam \
			  -M \
			  -O \
		 	  -s 0 \
			  -T 16;

              
# extract just gene counts for exploratory analysis
cut -d"   " -f1,7 gene_assigned.txt > gene_assigned_counts.txt;

# sort and index the counted reads
samtools sort -@ 16 -m 2G Aligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam;
samtools index assigned_sorted.bam;

# count and deduplicate per gene from bam file
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell -I assigned_sorted.bam -L count.log -S counts_long.tsv;
