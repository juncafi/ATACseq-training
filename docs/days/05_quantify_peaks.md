## 1. Count reads on peaks

After defining a consensus peak set, we will count the number of reads that overlap each peak in each sample. This step can be seen as counting reads in annotated genes for RNA-seq data, where instead of genes, we are counting reads in peaks.

There are different tools that can be used for this purpose, such as featureCounts from the Subread package or bedtools coverage. In this example, we will use featureCounts, which is a widely used tool for counting reads in genomic features. 

**Task 1**

We first need to convert the consensus peaks to SAF format, which is the format required by featureCounts
The SAF format is a tab-delimited file with the following columns: GeneID, Chr, Start, End, Strand


- Create a new folder named: `results/05_counts_reads`
- Use awk or an alternative way to convert the bed file to SAF format. Don't forget to add a header "GeneID    Chr Start   End Strand". 

??? success "Solution"

    ```{bash}
    mkdir results/05_counts_reads
    echo "GeneID    Chr Start   End Strand" > results/04_consensus_peaks/consensus_peaks.saf
    awk '{OFS = "\t"} {print "Interval_"NR,$1,$2,$3,"."}' results/04_consensus_peaks/consensus_peaks.bed >> results/04_consensus_peaks/consensus_peaks.saf
    ```

Now we can use featureCounts to count the reads in each peak for each sample

**Task 2**

- Run featureCounts to count filtered reads (from "*qc_bl_filt.sorted.bam") into consensus peak annotation

!!! FeatureCounts
    -F SAF: specify that the annotation file is in SAF format  
    -p: specify that the input files are paired-end  
    -a: specify the annotation file  
    -o: specify the output file  
    
    These are only some parameteres that we will apply. You can also adjust the parameters of featureCounts based on your specific requirements, such as setting a minimum mapping quality or handling multi-mapping reads.


??? success "Solution"

    ```{bash}
    path_bams="results/01_filtered_bams"
    featureCounts -F SAF -T 2 -p -a results/04_consensus_peaks/consensus_peaks.saf -o results/05_counts_reads/feature_Counts.txt $path_bams/*qc_bl_filt.sorted.bam 2> results/05_counts_reads/featureCounts.log
    ```

The output file will contain the counts of reads in each peak for each sample, which can be used for downstream analysis such as differential accessibility analysis.

## 2. Assess results and overall quality
After filtering, peak calling and counting reads on peaks, we can run MultiQC to aggregate the QC metrics from different steps and generate a comprehensive report.  
This will help us to assess further he overall quality of the ATAC-seq data and identify any potential issues that may need to be addressed before proceeding with downstream analysis.  

**Task 3**  

- Run multiqc on the entire `results` directory


??? success "Solution"

    ```{bash}
    multiqc --outdir results/multiQC_report --title multiQC_post_read_counting results

    ```
    
