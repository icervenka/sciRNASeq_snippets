# Processing of sciRNASeq data

This a is basic shell based pipeline to process single-cell combinatorial indexing RNASeq data (sciRNASeq) as described in [Comprehensive single-cell transcriptional profiling of a multicellular organism](https://science.sciencemag.org/content/357/6352/661.full) (Cao et. al 2017) that identifies single cell transcriptomes without single cell isolation.

For the overview of the library generation and sequencing, please refer to the [following schematics](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq.html), unlike regular scRNASeq, this one requires 2-step barcode demultiplexing. Example data is available [in GEO archive](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2599700).

## Usage
`sciRNAseq.sh -i . -o ./results -x ~/star_index -g ~/annotation.gtf`

## Commandline arguments
- `-i` input directory with fastq files, requires *R1.fastq.gz and *R2.fastq.gz files where UMI barcodes are stored in R1 file
- `-o` output directory for results
- `-x` location of genome index for STAR aligner, example code for index generation is included sciRNAseq.sh
- `-g` location of corresponding gtf file

## Known issues
Currently umi tools estimates the barcode whitelist from the data and user supplied barcodes are not currently supported.

## Dependencies

`cutadapt, STAR, samtools, featureCounts, umi_tools, pigz`
