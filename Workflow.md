**Step 1: Preparing .gtf files (optional)**

Convert .gff files provided by microbesNG to .gtf files so they can be used with STAR

>Note: This .gtf file does not carry over uniprotKB annotations but they will be merged in later steps in R



```
gffread genome1.gff -T -o genome1.gtf
```
**Step 2: Use STAR to map bacterial transcripts to bacterial genomes**

Build a genome index

_Repeat for all bacterial strains_

```
STAR \
  --runMode genomeGenerate \
  --genomeDir APLE2372_index \
  --genomeFastaFiles genome1.fasta \
  --sjdbGTFfile genome1.gtf \
  --sjdbOverhang 6328143 \
  --runThreadN 40 \
  --genomeSAindexNbases 10
```

Map each set of transcripts to the genome 

_Repeat for all transcript sets_
```
STAR --runThreadN 40 \
  --genomeDir genome1_index \
  --readFilesIn transcript1_forward.fq.gz transcript1_reverse.fq.gz \
  --readFilesCommand gunzip \
  -c \
  --alignIntronMax 1 \
  --outFileNamePrefix ./genome1/genome1T1 \
  --quantMode GeneCounts \
```

**Step 3: Use FeatureCounts to count the number of reads that map to each genomic feature**

_repeat for all bacterial strains_

```
cd genome1
featureCounts -a genome1.gtf -o genome1_counts.txt -T 4 -s 2 -p -B -C -g gene_id -O -f -F GTF -t transcript genome1T1Aligned.out.sam 
cd ..
```

**Step 4: EdgeR gene expression analysis**
Load the appropriate packages 
```
library(Rsubread)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(edgeR)
```

Read the counts file and remove the first 5 columns of gene annotation

```
counts_genome1 <- read.table("genome1_counts.txt", header = TRUE, row.names = 1, sep = "\t")
head(counts_genome1)
counts_data_genome1 <- counts_genome1[, -c(1:5)]
```
Make a DGElist data class depending on treatment for each sample in the counts file 
```
group_1 <- c("NT", "NT", "NT", "T", "T", "T")
dge_1 <- DGEList(counts<-counts_data_genome1, group = group1)
```
Filter out lowly expressed genes

```
keep_1 <- filterByExpr(dge_1)
dge_1 <- dge_1[keep, , keep.lib.sizes=FALSE]
```
Normalise library sizes to minimise log-fold changes between samples of most genes 
_Below uses edgeR's default method of trimmed mean of M-values (TMM) between each pair of samples and is recommended if least half of the genes in your RNAseq data are likely to not be differentially expressed_

```
dge_1 <- normLibSizes(dge_1)
```

Create an MDS plot of each sample to check for outliers 
```
label_colors <- c(rep("red", 3), rep("blue", 3))
plotMDS(dge_1, labels=c("G1_NT1", "G1_NT3", "G1_NT4", "G1_T1", "G1_T3", "G1_T4"), cex = 1.25, col = label_colors)
```
Estimate dispersion 
_This is needed to be able for the following steps_
```
design_1 <- model.matrix(~group_1)
dge_1 <- estimateDisp(dge_1, design_1)
```
_Optional: plot the dispersion_
```
plotBCV(dge_1)
```

