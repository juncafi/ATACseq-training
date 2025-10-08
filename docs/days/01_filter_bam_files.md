## 0. Dataset

We are going to work with a subset of the publicly available dataset from [Liu et al. 2019](https://www.nature.com/articles/s41597-019-0071-0).


We are going to analyse and compare ATAC-seq data from 2 different adult mouse tissues:      
**Kidney**: *Rep1*, *Rep2*   
**Cerebrum**: *Rep1*, *Rep2* 

> Raw data can be found in [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/SRP167062)


  
In order to save time, raw reads have already been trimmed using [Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) and mapped to the reference genome (mm10) using [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) in end-to-end mode. Therefore, we will start the analysis directly from the alignment output (**.bam** files).  

To avoid large waiting times during high-demanding computational steps, .bam files have been subset to keep only reads aligning to chromosome 6.  

!!! note

    You can find the .bam files in `/data/Liu_alignments_chr6/` folder


## 1. Filtering reads from alignments

### General NGS filtering steps 

In order to keep good quality data, read filtering criteria can be applied. These creteria are similar to the ones we would apply when analysing any other NGS related dataset (ie. remove duplicated reads, read mapping quality thresholds, etc.). They will depend on the quality of the data and the downstream analysis that will be applied.

We can use tools like [picard](https://broadinstitute.github.io/picard/) to indetify and mark read duplicates, which originate most likely from PCR amplifications. This step has already been performed for you, and read duplicates have been marked in the provided .bam files. If you want to know more about this tool and how it works you can have a look [here](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates).

**Task 1:** 

- Create new directory for filtered bams: `results/01_filtered_bams`

- Using `samtools view`, select/filter reads from .bam files in order to apply the following criteria: 

    * Keep only paired reads 
    * Remove unmapped reads
    * Remove reads who's mate is unmapped
    * Remove read duplicates
    * Remove not primary alignment reads
    * Keep mapping quality > 10

  Save the results in `results/01_filtered_bams` in BAM format with the following output name: `<sample name>.qc_filt.bam`

!!! Samtools view useful parameters

    Use the samtools view manual [here](https://www.htslib.org/doc/samtools-view.html) to understand which parameters you need, and this useful [page](https://broadinstitute.github.io/picard/explain-flags.html) from the broad institue to find what SAM Flag value would correspond to the combination of criteria we want to apply. 

<details>
<summary>Hint</summary>
You can use:  
 `-f` and `-F` flags to filter reads based on a combination of SAM Flags  
 `-q` for MapQ  
 `-h` to keep the .bam header  
 `-b` to get output in .bam format  
</details>

??? success "Solution"

    ```{bash}
    # create new directory for filtered bams

    mkdir -p results/01_filtered_bams

    # filter files using a for loop

    for bam in /data/Liu_alignments_chr6/*.bam; do
        echo "Processing file: $bam"
        sample_name=$(basename "$bam" .bam) # extract sample name without path and extension
        samtools view -h -b -q 10 -f 1 -F 1292 -o results/01_filtered_bams/$sample_name.qc_filt.bam $bam
    done
    ```

    -h: keep header; -b: output in bam format; -q 10: minimum mapping quality 10; -f 1: keep paired; -F 1292: exclude reads with any of the following flags: read unmapped, mate unmapped, not primary alignment, read is duplicate


After filtering we will sort and index the bam files for the next step


**Task 2:** 

- Use `samtools sort` to sort the filtered bam files (from previous task: `<sample name>.qc_filt.bam`). Save them in the same folder as: `<sample name>.qc_filt.sorted.bam`
- Use `samtools index` to index the sorted bam files 


??? success "Solution"

    ```{bash}
    for bam in results/01_filtered_bams/*qc_filt.bam; do
        echo "Processing file: $bam"
        samtools sort -o ${bam%.bam}.sorted.bam $bam
        samtools index ${bam%.bam}.sorted.bam
    done
    ``` 


To avoid confusion and large size files, we will keep only the sorted and indexed bams:
```{bash}
rm results/01_filtered_bams/*qc_filt.bam
```



### ATACseq related filtering steps

Mitochondria DNA is nucleosome free, therefore it is more accessible for Tn5 and several reads may have originated from mitochondrial DNA. After having assessed mitochondrial % on the QC, we can discard reads coming from chmMT to avoid biases in downstream analsis.

Since we are working only with chm 13 we don't need to do this step, but here is the command you could use for that:
```{bash}
samtools view -h input.bam | awk  '($3 != "MT")' | samtools view -hb - > output.bam
```
!!! Note
    The mitochondrial chromosome name may differ depending on the reference genome (e.g., "MT", "chrM", "chrMT").

Next, we will remove reads overlapping problematic regions of the genome. ENCODE consortium has created comprehensive lists of such regions (anomalous, unstructured or high signal in NGS experiments) for different genome species (including mouse mm10). These lists are called ENCODE Blacklists, and you can find them [here](https://github.com/Boyle-Lab/Blacklist/). 

!!! Note
    The regions for mm10 have been dowloaded as .bed file, and you will find it here: `/data/references/mm10-blacklist.v2.nochr.bed`

**Task 3**


- Using `bedtools intersect` filter out reads overlapping regions in the Blacklist .bed file
- Use previously filtered .bam files as input: `results/01_filtered_bams/*qc_filt.sorted.bam`, don't forget to specify your input is in BAM format
- Use the `data/references/mm10-blacklist.v2.nochr.bed` as regions to filter out reads from (it is already sorted)
- Save the results in `results/01_filtered_bams` in BAM format with the following output name: `<sample name>.qc_bl_filt.bam` 
- Sort and index the resulting .bam with `samtools` and save them in `results/01_filtered_bams` with name  `<sample name>.qc_bl_filt.sorted.bam` 

!!! Bedtools intersect useful information
    Bedtools intersect allows one to screen for overlaps between two sets of genomic features/regions, and then decide on which kind of information do you want to report.
    Here, we will intersect: a) aligned reads to the genome (information in your .bam files)l b) Problematic regions listed in Blacklist .bed file
    We do not want to keep the reads that overlap Blacklist regions.
    You can find documentation on which parameters to use [here](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)



??? success "Solution"
    ```{bash}
    blacklist="/data/references/mm10-blacklist.v2.nochr.bed"

    for bam in results/01_filtered_bams/*qc_filt.sorted.bam; do
        echo "Processing file: $bam"
        sample_name=$(basename "$bam" .qc_filt.sorted.bam) 
        bedtools intersect -v -abam $bam -b $blacklist > results/01_filtered_bams/$sample_name.qc_bl_filt.bam
        samtools sort -o results/01_filtered_bams/$sample_name.qc_bl_filt.sorted.bam results/01_filtered_bams/$sample_name.qc_bl_filt.bam
        samtools index results/01_filtered_bams/$sample_name.qc_bl_filt.sorted.bam
    done
    ```

     Parameter explanation:  
     -v: only report those entries in A that have no overlap with B  
     -abam: input is in bam format 


After sorting the .bam files, we don't need the unsorted .bam. To free some space we will remove the unsorted bam files (intermediate files).

**Task 4**

- Remove intermediate files that we don't need anymore. 

!!! Warning
    Run this script only if you named the files in the same way the tutorial mentioned, and only if you have finished the sorting and indexing steps from Task 2 and 3

```{bash}
    rm results/01_filtered_bams/*qc_filt.bam
    rm results/01_filtered_bams/*qc_bl_filt.bam
```
