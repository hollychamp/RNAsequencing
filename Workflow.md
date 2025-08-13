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

_Repeat for all bacterial strains_

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
_This is needed to be able fit the quasi-likelihood negative binomal generalised linear model_
```
design_1 <- model.matrix(~group_1)
dge_1 <- estimateDisp(dge_1, design_1)
```
_Optional: plot the dispersion_
```
plotBCV(dge_1)
```
Use trended NB dispersion and fit the QL GLM
```
fit_1 <- glmQLFit(dge_1, design_1)
```
Use the quasi-likelihood f-test to determine differential expression 
_'coef' specifies which group you would like to test for differential gene expression, below I used coef = 2 to test whether the treated samples (second line in my design matrix) is significantly different from the baseline (group 1)_
```
DE_1 <- glmQLFTest(fit_1, coef = 2)
```
Use the Benjamini-Hochberg procedure to adjust p-values for all genes - these adjusted p-values will be the False Discovery Rate
```
DE_1$table$FDR <- p.adjust(DE_1$table$PValue, method = 'BH')
```
**Step 5: Analysing data**

Create a volcano plot highlighting genes that have an FDR < 0.05 and an absolute log fold change of >1
```
ggplot(DE_1$table, aes(
  x = logFC,
  y = -log10(FDR),
  color = case_when(
    FDR < 0.05 & logFC > 1 ~ "Upregulated",
    FDR < 0.05 & logFC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  )
)) +
  geom_point(size = 2, alpha = 0.8, shape = 16) +
  scale_color_manual(
    values = c("Upregulated" = "#D32F2F", "Downregulated" = "#1976D2", "Not Significant" = "gray80")
  ) +
  labs(
    title = "Strain 1 ",
    x = "logFC",
    y = expression(-log[10](FDR)),
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )
```
Calculate the number of significantly upregulated and downregulated genes (FDR < 0.05 and absolute logFC > 1)
```
sum(DE_1$table$FDR < 0.05 & (DE_1$table$logFC) > 1)
sum(DE_1$table$FDR < 0.05 & (DE_1$table$logFC) < -1)
```
Filter DE_1 to create a data frame that only includes genes with an FDR <0.05. The order this list by logFC
```
sig_DE_1 <- DE_1$table %>% filter(FDR < 0.05)
sig_DE_1_ordered <- sig_DE_1[order(sig_DE_1$logFC, decreasing=TRUE),]
```
Add a gene ID column to the sig_DE_1_ordered so we can merge two data frames by gene ID in later steps 
```
sig_DE_1_ordered$gene_id <- rownames(sig_DE_1_ordered)
```
Import the .gff file from microbesNG as a data frame to be able to get product names of each gene ID
```
gff_file_1 <- "genome1.gff"
gff_data_1<- import.gff(gff_file_1, format="gff")
gff_data_1 <- as.data.frame(gff_data_1)
```
Rename the ID column to gene_ID (to align with the format of the sig_DE_1_ordered data frame)
```
colnames(gff_data_1)[colnames(gff_data_1) == "ID"] ="gene_id"
```
Make a data frame with gene ID and product name only 
```
gene_annots_1 <- gff_data_1 %>%
  dplyr::select(gene_id, product) %>%
  dplyr::distinct()
```
Merge gene_annots_1 with sig_DE_1_ordered
```
sig_DE_1_annotated <- sig_DE_1_ordered  %>%
  left_join(gene_annots_1, by = "gene_id")
```
Create a data frame that lists the log counts per million (logCPM) of all genes per sample 
```
log_cpm_1 <- cpm(dge_1, log = TRUE)
```
**Step 5.5: Create a heatmap of the top 30 significantly enriched genes**

Extract gene_ids and product names of top 30 significantly enriched genes (FDR<0.05)
```
top_DE_genes_1 <- head(sig_DE_1_annotated, 30)
top_30_gene_id_1 <- top_DE_genes_1$gene_id
```
Subset logCPM data to only include the top 30 genes 
```
top30_log_cpm_1 <- log_cpm_1[rownames(log_cpm_1) %in% top_30_gene_id_1, ]
```
Create a vector which maps gene IDs to product name 
```
id_to_product_1 <- setNames(top_DE_genes_1$product, top_DE_genes_1$gene_ID)
```
Replace rownames in top30_log_cpm_1 to product names 
```
rownames(top30_log_cpm_1) <- paste(id_to_product_1)
```
Create a heatmap of the top 30 enriched genes 
```
Library(pheatmap)
heatmap_1_top30 <- pheatmap(top30_log_cpm_1, fontsize_row = 18, cluster_cols = FALSE, cluster_rows = TRUE)
```



