# Peak Annotation and Functional Analysis

After performing differential accessibility analysis with DESeq2, we will annotate peaks based on their genomic location using ChIPseeker and perform functional enrichment analysis.

## Overview

In this section, we will:
1. Load the necessary libraries for peak annotation
2. Prepare our differential accessibility results for annotation
3. Visualize peak distribution across the genome
4. Annotate peaks with genomic features
5. Perform functional enrichment analysis

## 1. Load Required Libraries

First, let's load the necessary R packages for peak annotation and functional analysis:

```{r}
# Load libraries
library(ChIPseeker)
require("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(clusterProfiler)
library("org.Mm.eg.db")

# Set up annotation databases
TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
AnnoDb <- 'org.Mm.eg.db'
```


## 2. Prepare Differential Accessibility Results

Now let's prepare our data by combining the differential accessibility statistics with the count data:

```{r}
# Load again the counts table
counts_file <- "results/05_counts_reads/feature_Counts.txt"
cts_table <- read.table(counts_file, header = T)

# Load the results from DESeq2 and extract relevant columns from DESeq2 results
res <- read.table("results/06_DA_analysis/DA_results.txt")
DA_stats <- as.data.frame(res)[c(2,3,5,6)]

# Clean column names from count table
colnames(cts_table) <- gsub(colnames(cts_table), pattern = "results.01_filtered_bams.", replacement = "") %>% 
  gsub(., pattern=".qc_bl_filt.sorted.bam", replacement="")

# Merge differential accessibility stats with count data
peaks_df <- merge(DA_stats, cts_table, by.x="row.names", by.y="Geneid")

# Add 'chr' prefix to chromosome names for consistency
peaks_df$Chr <- paste0("chr", peaks_df$Chr)

# Create GRanges object for downstream analysis
gr <- makeGRangesFromDataFrame(peaks_df, keep.extra.columns=TRUE)
```

## 3. Visualize Peak Distribution

Let's visualize how our peaks are distributed across the genome and their fold changes:

```{r}
# Create coverage plot showing peak positions and their fold changes
covplot(gr, weightCol = "log2FoldChange")
```

This plot shows the distribution of peaks across chromosomes, with y-axis representing the log2 fold change values.

## 4. Prepare Data for Annotation

```{r}
# Sort the GRanges object by genomic coordinates
gr <- sort(gr, by = ~ seqnames + start + end)

# Split peaks into upregulated and downregulated based on significance and fold change
gr_list <- list(
  up = gr[gr$padj < 0.01 & gr$log2FoldChange > 2,], 
  down = gr[gr$padj < 0.01 & gr$log2FoldChange < -2,]
)




# Check the structure of our peak lists
head(gr_list)
```

**Task 1** 

- Check how many significantly upregulated and downregulated peaks we have:

```{r}
# Count peaks in each category
cat("Number of upregulated peaks:", length(gr_list$up), "\n")
cat("Number of downregulated peaks:", length(gr_list$down), "\n")
```

## 5. Annotate Genomic Overlap with ChIPseeker

Now we'll annotate our peaks to determine their genomic context (promoters, exons, introns, etc.):

```{r}
# Annotate upregulated peaks
peakAnno_up <- annotatePeak(gr_list$up, 
                           tssRegion = c(-1000, 1000), 
                           TxDb = TxDb, 
                           annoDb = AnnoDb, 
                           overlap = "TSS")

# Annotate downregulated peaks
peakAnno_down <- annotatePeak(gr_list$down, 
                             tssRegion = c(-1000, 1000), 
                             TxDb = TxDb, 
                             annoDb = AnnoDb, 
                             overlap = "TSS")

                             
# Save the objects
# Save the objects
saveRDS(peakAnno_up, "results/06_DA_analysis/Annotated_peaks_up.rds")
saveRDS(peakAnno_down, "results/06_DA_analysis/Annotated_peaks_down.rds")
saveRDS(gr, "results/06_DA_analysis/gr_peaks.rds")
```

**Task 2**  

- Examine the annotation results


```{r}
# Look at the upregulated peaks annotation summary
peakAnno_up

# Count how many downregulated peaks are in promoter regions
sum(peakAnno_down@detailGenomicAnnotation$Promoter)
```

### Visualize Peak Annotations

Create pie charts to visualize the genomic distribution of our peaks:

```{r}
# Plot pie charts for peak categories
plotAnnoPie(peakAnno_up, main = "\n\nDA up") 
plotAnnoPie(peakAnno_down, main = "\n\nDA down")
```

These plots show the percentage of peaks falling into different genomic categories (promoters, exons, introns, intergenic regions, etc.).

## 6. Functional Enrichment Analysis

Next, we will focus in peaks overlaping promoter regions (TSS) of genes, and we will perform functional enrichment analysis on those genes to understand which biological processes or metabolic pathways may be afected by the changes in chromatin accessibility.

### Extract Genes Associated with Promoter Peaks

```{r}
# Extract gene symbols for peaks in promoter regions
genes_tss_up <- peakAnno_up@anno[peakAnno_up@detailGenomicAnnotation$Promoter, "SYMBOL"]
genes_tss_down <- peakAnno_down@anno[peakAnno_down@detailGenomicAnnotation$Promoter, "SYMBOL"]

# Create universe of genes with all promoter-overlapping peaks for background
peakAnno <- annotatePeak(gr, 
                        tssRegion = c(-1000, 1000), 
                        TxDb = TxDb, 
                        annoDb = AnnoDb, 
                        overlap = "TSS")
universe <- peakAnno@anno[peakAnno@detailGenomicAnnotation$Promoter,]$SYMBOL
```

### Gene Ontology (GO) Enrichment Analysis

Perform GO enrichment analysis for upregulated peaks:

```{r}
# GO enrichment for upregulated genes
ego_up <- enrichGO(gene          = genes_tss_up$SYMBOL,
                   universe      = unique(universe),
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "SYMBOL",
                   ont           = "BP",  # Molecular Function
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# Plot results
ego_up
```

**Task 3**  

- Did you find significantly enriched GO terms?  

- Try to run the same analysis considering all genes as a universe, why do you think there is a difference?  

- Visualise the GO terms using a function from [ClusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html) such as: `barplot()`


??? success "Solution"
    ```{r}
    # GO enrichment for upregulated genes
    ego_up <- enrichGO(gene          = genes_tss_up$SYMBOL,
                    universe      = unique(universe),
                    OrgDb         = org.Mm.eg.db,
                    keyType       = "SYMBOL",
                    ont           = "BP",  # Molecular Function
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)

    # Plot results
    print(barplot(ego_down, showCategory = 20))
    ```

**Task 4**

- Perform GO enrichment analysis for downregulated peaks using all genes as universe:

??? success "Solution"
    ```{r}
    # GO enrichment for downregulated genes
    ego_down <- enrichGO(gene          = genes_tss_down$SYMBOL,
                        OrgDb         = org.Mm.eg.db,
                        keyType       = "SYMBOL",
                        ont           = "BP",  # Molecular Function
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

    # Plot results
    print(barplot(ego_down, showCategory = 20))
    ```



!!! tip "Next Steps"
    You can change the `ont` parameter to explore different aspects:
    - "MF": Molecular Function
    - "BP": Biological Process  
    - "CC": Cellular Component

## 7. Additional Visualizations

### Peak Heatmap

Create a heatmap showing peak signal around gene bodies:

```{r}
# Generate peak heatmap around gene bodies
peakHeatmap(peak = gr,
            TxDb = TxDb,
            upstream = rel(0.2),
            downstream = rel(0.2),
            by = "gene",
            type = "body",
            nbin = 800)
```

### Load Additional Reference Data

```{r}
# Load TSS reference data for additional analysis
tss_bed <- read.table("data/references/ENCODE_mm10_M21_TSS_reference.bed")
```

## Summary

In this tutorial, we have:

1. ✅ Loaded necessary libraries for peak annotation
2. ✅ Prepared differential accessibility results for analysis
3. ✅ Visualized peak distribution across the genome
4. ✅ Annotated peaks with genomic features using ChIPseeker
5. ✅ Performed functional enrichment analysis on genes associated with promoter peaks
6. ✅ Created visualizations to interpret the biological significance of our results

The results from this analysis help us understand:
- Where our differentially accessible peaks are located in the genome
- Which genes might be affected by changes in chromatin accessibility
- What biological processes or molecular functions are enriched in our peak sets

!!! tip "Next Steps"
    Consider exploring different GO ontology categories (BP, CC) and other enrichment databases like KEGG pathways for a more comprehensive functional analysis.
