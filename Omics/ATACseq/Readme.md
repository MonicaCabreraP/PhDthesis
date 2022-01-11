# ATAC-seq pipeline
--------------------

**Preprocessing of the p53 ATAC-seq datasets**
*The ATAC-seq data were processed (trimmed, aligned, filtered, and quality controlled) using the ATAC-seq pipeline from the Kundaje lab31,32 (Table 2). The model-based analysis of ChIP-seq (MACS2)33 version 2.1.2 was used to identify the peak regions with options -B, -q 0.01–nomodel, -f BAM, -g mm. The Irreproducible Discovery Rate (IDR) method34 was used to identify reproducible peaks between two technical replicates (Fig. 1b). Only peaks reproducible between the two technical replicates were retained for downstream analyses. Peaks for all tissues were then merged together into a standard peak list. The number of raw reads mapped to each standard peak were counted using the intersect function of BedTools35 version 2.26.0. The raw count matrix32 was normalized by Reads Per Million mapped reads (RPM). Pearson correlation coefficients between technical or biological replicates across tissues were calculated based on the Log10 RPM matrix.*

## Introduction
This pipeline is designed for automated end-to-end quality control and processing of ATAC-seq data.  The pipeline can be run on compute clusters using parallelized computing with job submission engines as well as on stand alone machines modifying the cluster job part. The pipeline can be run end-to-end, starting from raw FASTQ files all the way to peak calling and signal track generation using a single command. One can also start the pipeline from intermediate stages (for example, using alignment files as input). The pipeline supports paired-end data as well as replicated or non-replicated datasets. 

There are 4 important steps:
- [Pre-analysis](#Pre-analysis)

The results produced by the pipeline include: 
1) formatted HTML reports that include quality control measures specifically designed for ATAC-seq data
2) analysis of reproducibility
3) stringent and relaxed thresholding of peaks
4) fold-enrichment and pvalue signal tracks. 


(1) pre-analysis (quality control (QC) and alignment)
(2) core analysis (peak calling)
(3) advanced analysis at the level of peaks, motifs, nucleosomes, and TF footprints
(4) integration with multiomics data to reconstruct regulatory networks. 


## Pre-analysis

### Pre-alignment quality control
------------------------------------

The pre-alignment QC and read alignment steps are standard for most high-throughput sequencing technologies. For example, FastQC [39] can be used to visualize base quality scores, GC content, sequence length distribution, sequence duplication levels, k-mer overrepresentation and contamination of primers and adapters in the sequencing data. An overall high base quality score with a slight drop towards the 3′ end of sequencing reads is acceptable. No obvious deviation from expected GC content and sequence read length should be observed. Moreover, the metrics should be homogeneous among all samples from the same experimental batch and sequencing run.

Currently, due to the ubiquitous use of Illumina’s Nextera library for ATAC-seq, overrepresentation of Nextera sequencing adapters is often observed and should be removed for accurate read alignment. Most adapter removal tools employ different variations of dynamic programming, such as cutadapt [40], AdapterRemoval v2 [41], Skewer [42], and trimmomatic [43] all requiring input of known adapter sequences. For example, using trimmomatic with built-in adapter sequences for Nextera and TruSeq library would be a straightforward step. Low-quality bases can also be eliminated using these tools. From our experience, read trimming tools are generally comparable in performance of efficient removal of low-quality and contaminating adapter sequences.
Alignment

After read trimming, FastQC can be performed again to check the successful removal of adapter and low-quality bases. Trimmed reads are then mapped to a reference genome. BWA-MEM [44] and Bowtie2 [45] aligners are memory-efficient and fast for short paired-end reads. The soft-clip strategy from both aligners allows the overhang of bases on both ends of reads which can further increase unique mapping rates [46]. We suggest that a unique mapping rate over 80% is typical for a successful ATAC-seq experiment. For mammalian species, the suggested minimum number of mapped reads is 50 million for open chromatin detection and differential analysis, and 200 million for TF footprinting based on empirical and computational estimations [10, 12, 47,48,49].

**Post-alignment processing and quality control**

After sequence alignment, as in most DNA sequencing data, basic metrics of the aligned BAM file, such as unique mapping reads/rates, duplicated read percentages, and fragment size distribution can be collected using Picard [50] and SAMtools [51]. Additionally, reads should be removed if they are improperly paired or of low mapping quality. The mitochondrial genome, which is more accessible due to the lack of chromatin packaging [52], and the ENCODE blacklisted regions [53, 54] often have extremely high read coverage, and should also be discarded [33]. Duplicated reads, which are likely to have arisen as PCR artifacts, should also be removed to significantly improve biological reproducibility [48]. These steps will together improve the power of open chromatin detection and produce fewer false positives.

There are additional ATAC-seq-specific quality metrics that need to be evaluated. Typically, a successful ATAC-seq experiment should generate a fragment size distribution plot with decreasing and periodical peaks corresponding to the nucleosome-free regions (NFR) (< 100 bp) and mono-, di-, and tri-nucleosomes (~ 200, 400, 600 bp, respectively) (Fig. 1b) [9, 55]. Fragments from the NFR are expected to be enriched around the transcription start site (TSS) of genes, while fragments from nucleosome-bound regions are expected to be depleted at TSS with a slight enrichment of flanking regions around TSS (Fig. 1c) [55]. These can be evaluated with the tool ATACseqQC [55]. Lastly, reads should be shifted + 4 bp and − 5 bp for positive and negative strand respectively, to account for the 9-bp duplication created by DNA repair of the nick by Tn5 transposase and achieve base-pair resolution of TF footprint and motif-related analyses [9, 33, 56]. Most aforementioned QC and analysis reports can be integrated using MultiQC [57] for an aggregated, user-friendly, and interactive presentation.

A major consideration for the appropriate tools to choose here is often time to result. Read trimming and alignment can be time consuming, and there is always a trade-off between speed and accuracy. In our experience, the following pipeline performs reasonably well: FastQC➔ trimmomatic➔BWA-MEM➔ATACseqQC, and we would suggest this as a good starting point for processing of ATAC-seq data.
