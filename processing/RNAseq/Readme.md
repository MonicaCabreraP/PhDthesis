# RNA-seq processing

Detailed description of how RNA-seq data was processed in this thesis:

## Dependencies

* [Trim Galore](https://github.com/FelixKrueger/TrimGalore)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [STAR](https://github.com/alexdobin/STAR/tree/master)
* [sambamba](https://github.com/biod/sambamba)
* [featureCounts](https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html)
* [DESeq2 (R package)](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

## Summary of the workflow

1. **Quality Control**: FastQC
2. **Trimming**: Trim Galore
3. **Alignment**: STAR
4. **Read Counting**: featureCounts
5. **Differential Expression Analysis**: DESeq2

---

## 1. Quality Control (FastQC)
FastQC performs quality assessment on raw sequencing data.

```bash
fastqc sample1.fq.gz -o /path/quality/
```
**Expected Results**: All green flags; for explanation of each quality check you can read the FastQC manual. NOTE: Red flags in Sequence Duplication Levels are typical for RNA-seq due to varying transcript expression levels.

## 2. Trimming (Trim Galore)
Trim Galore removes adapter sequences and poor-quality ends from reads, enhancing data quality for downstream analysis.

```bash
trim_galore $FASTQ1 $FASTQ2 --paired --output_dir $FASTQ_DIR --basename $NAME -a $ADAPTERS -A $ADAPTERS --fastqc_args '--outdir $OUTDIR/quality' --cores $THREADS
```
## 3. Alignment (STAR)

### Preparing Genome Index
Download reference genome [GRCh37](https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/), and  the [GTF](http://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/) file of this version.

Generates genome index using STAR, specifying genome FASTA files and GTF annotations.
Note: Adjust --sjdbOverhang based on read length (read length - 1).

```bash
STAR --runThreadN 20 --runMode genomeGenerate \
     --genomeDir $(pwd) --genomeFastaFiles chr1.fa chr2.fa ... \
     --sjdbGTFfile /path/GTF --sjdbOverhang 99 \
     --sjdbGTFtagExonParentTranscript transcript_id \
     --sjdbGTFtagExonParentGene gene_id
```

###Running STAR Alignment
Align trimmed reads to the indexed reference genome using STAR, following ENCODE-recommended parameters and converts SAM to BAM format using samtools.

```bash
STAR --runThreadN 15 --genomeDir /path/to/index/ \
     --readFilesIn sample1.fq.gz sample2.fq.gz --outFilterType BySJout \
     --readFilesCommand zcat --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
     --outFileNamePrefix /path/alignment/sample1 --outStd SAM | \
     samtools view -bS - > /path/alignment/sample1.bam
```

## 3. Read Counting (featureCounts)

### Sorting BAM File
Sort BAM file for efficient read counting.

```bash
sambamba sort -t 15 -o /gpfs/projects/bsc08/bsc08471/p53/results/HCT116/Omics/RNA/HCT116_WTNutlin1h/results/RNA/alignment/HCT116_WTNutlin1h_RNA_1.sort.bam /gpfs/projects/bsc08/bsc08471/p53/results/HCT116/Omics/RNA/HCT116_WTNutlin1h/results/RNA/alignment/HCT116_WTNutlin1h_RNA_1.bam
```
### Counting Reads
Use featureCounts to generate a count matrix for downstream analysis, specifying GTF and BAM files.

```bash
featureCounts -a $GTF -p -B -C -T 15 -o $OUTFILE.csv $BAM1 $BAM2
```

## 4. Differential Expression Analysis (DESeq2)

### Creating Count Matrix
Merge individual sample count files into a unified count matrix, preparing for DESeq2 analysis.

```r
# Merge count files into a single count matrix
files <- list.files(path = "/path/to/featureCounts/out", pattern = "*.csv$", recursive = TRUE, full.names = TRUE)
countdata <- matrix(data.table::fread(files[1])$Geneid)
colnames(countdata) <- "Geneid"

for (i in files) {
  a <- data.table::fread(i)
  a <- a[, -(2:6)]
  colnames(a)[-1] <- gsub("\\.[sb]am$", "", basename(colnames(a)[-1]))
  a <- as.matrix(a)
  countdata <- merge(countdata, a, by = "Geneid")
}

countdata[, -1] <- apply(countdata[, -1], 2, as.numeric)

# Write count data to file
data.table::fwrite(countdata, file = "count_table.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
```

### Differential Expression Analysis
Normalize data for library size factors using estimateSizeFactors.
Perform differential expression analysis using DESeq2.
Summarize and identify significant genes based on fold change and adjusted p-values.

```r
# Load DESeq2 library
library(DESeq2)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = countdata[, -1], colData = coldata, design = ~condition)

# Normalize for library size factors
dds <- estimateSizeFactors(dds)

# Perform differential expression analysis
dds <- DESeq(dds)

# Extract significant results
res <- results(dds, contrast = c("condition", "cond1", "cond2"))

# Summarize differential expression results
summary(res)

# Identify significantly upregulated and downregulated genes
up <- rownames(res)[which(res$log2FoldChange > 2 & res$padj < 0.05)]
down <- rownames(res)[which(res$log2FoldChange < -2 & res$padj < 0.05)]

```
