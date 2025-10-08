In the next step, we will start importing the read counts on consensus peaks (output from featureCounts) into R and perform differential accessibility analysis using DESeq2.

## 1. Import read counts and metadata

We start by loading featureCounts output as a matrix of counts

Load necessary packages for the analysis
    ```r
    library(DESeq2)
    library(ggplot2)
    library(dplyr)
    ```

**Task 1**

- Read the output of featurecounts that contains counts for all peaks in all samples: `results/05_counts_reads/feature_Counts.txt`
- Keep the header of the table
- Process the table in order to keep: read-counts columns only, peak IDs as row.names, simple column names (sample names)
- Convert the table into matrix format

??? success "Solution"

    ```r
    # Read table of counts (output from FeatureCounts)
    counts_file <- "results/05_counts_reads/feature_Counts.txt"
    cts_table <- read.table(counts_file, header = T)
    
    # convert into matrix, with interval IDs as rownames
    cts <- cts_table[,7:10]
    row.names(cts) <- cts_table[,1]
    colnames(cts) <- gsub(colnames(cts), pattern = "results.01_filtered_bams.", replacement = "") %>% 
                     gsub(., pattern=".qc_bl_filt.sorted.bam", replacement="")
    ```

**Task 2**

Create metadata table
We need to create a data.frame containing metadata information about the samples.
Here we have two conditions, PC and RACM, with two replicates each. We need to create a data.frame containing the conditions information for each sample. 
If we would have a second factor, e.g. treatment (treated vs untreated), we would need to add this information as a second column in the data.frame.

!!! Warning
  
    Important: The order of the rows in the data.frame must match the order of the columns in the counts matrix
               The metadata data.frame must have factor variables (not character variables)


```r
condition <- factor( c(rep("Cerebrum",2), rep("Kidney",2)) )
colData <- data.frame(condition, row.names = colnames(cts))
        
```

## 2. DEseq object and analysis
Bring together counts and metadata to create a DESeq object

```r
dds <- DESeqDataSetFromMatrix(
  countData = cts, colData = colData, 
  design = ~ condition)
dim(dds)
```

Remove lowly exp peaks

```r
idx <- rowSums(counts(dds, normalized=FALSE) >= 50) >= 2
dds.f <- dds[idx, ]
dim(dds.f)
```

We perform the estimation of dispersions

```r
dds <- DESeq(dds)
```

And plot PCA of the samples
```r
vsd <- varianceStabilizingTransformation(dds, blind=TRUE )
pcaData <- plotPCA(vsd, intgroup=c("condition"), ntop="all")
pcaData + geom_label(aes(x=PC1,y=PC2,label=name))
plotPCA(vsd, intgroup=c("sizeFactor"), ntop="all")
```
    

## 3. Save differential Accessibility results
After the dispersion estimates have been calculated, we can proceed to test for differential accessibility between the two conditions (PC vs RACM).  

The `results()` function extracts a results table with the log2 fold changes, p-values and adjusted p-values for each peak.
```r
res <- results(dds)
summary( res )
```

We will save the results table in a new directory
```r
dir.create("results/06_DA_analysis")
DA_results <- as.data.frame(res)
write.table(DA_results, file="results/06_DA_analysis/DA_results.txt", quote=FALSE)
```

!!! Note

    By default, the results() function will extract the results for the last variable in the design formula (here: condition) and will perform a comparison of the second level of the factor over the first level (here: RACM over PC).
    If you want to extract results for a different comparison, you can specify the contrast argument. Have a look at [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) documentation for more information.

We can visualise the results of DA peaks with a Volcano plot, analogous to RNAseq DE genes results

```r

```


### Create GRanges object
We will convert the peak coordinates into a GRanges object (from [GenomicRanges package](https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html)). 
Granges class represents a collection of genomic ranges and associated data to them. In this case, we will use it to represent the peak coordinates, and we will add as metadata the log2 fold changes and adjusted p-values from the differential accessibility analysis.


```r
    # Prepare peak coordinates object 
    peaks_coord <- cts_table[,1:4]

    # Select information from DESeq2 we want to keep
    DA_stats <- as.data.frame(res)[c(2,3,5,6)]

    # Merge both objects and convert them into GRanges class
    peaks_df <- merge(DA_stats, peaks_coord, by.x="row.names", by.y="Geneid")
    peaks_gr = makeGRangesFromDataFrame(peaks_df, keep.extra.columns=T)

    # Have a look into the new object
    head(peaks_gr)
```
To this object `peaks_gr`, we will add other metadata along the downstream analysis. 

```r
rm(list = ls())
```