## 1. Visualisation of peaks

We will start by having a look into the previous results. For that, we can use IGV

## 2. Buid a consensus annotation

We need to build a consensus peak set for downstream analysis (e.g. differential accessibility analysis), which will become like our reference annotation for counting reads in peaks.
We will do this by intersecting the peaks called in each replicate for each condition, and then merging the peaks from both conditions into a final consensus peak set.

Run the following commands to create intersected peak files for each condition (PC and RACM), requiring at least 25% overlap (-f 0.25) and reciprocal overlap (-r)

!!! Code
    Create a folder for consensus peaks. Use bedtools intersect to keep peaks present in both replicates

    ```{bash}
 
        mkdir results/04_consensus_peaks
        path_peaks="results/03_peak_calling/all_peaks"

        bedtools intersect -wa -a $path_peaks/all_peaks_Kidney_rep1/NA_peaks.narrowPeak -b $path_peaks/all_peaks_Kidney_rep2/NA_peaks.narrowPeak -f 0.25 -r > results/04_consensus_peaks/Kidney_intersect.bed
        bedtools intersect -wa -b $path_peaks/all_peaks_Kidney_rep1/NA_peaks.narrowPeak -a $path_peaks/all_peaks_Kidney_rep2/NA_peaks.narrowPeak -f 0.25 -r >> results/04_consensus_peaks/Kidney_intersect.bed
        bedtools intersect -wa -a $path_peaks/all_peaks_Cerebrum_rep1/NA_peaks.narrowPeak -b $path_peaks/all_peaks_Cerebrum_rep2/NA_peaks.narrowPeak -f 0.25 -r > results/04_consensus_peaks/Cerebrum_intersect.bed
        bedtools intersect -wa -b $path_peaks/all_peaks_Cerebrum_rep1/NA_peaks.narrowPeak -a $path_peaks/all_peaks_Cerebrum_rep2/NA_peaks.narrowPeak -f 0.25 -r >> results/04_consensus_peaks/Cerebrum_intersect.bed

    ```
    Parameters explanation:  
    -a: first input file (bed or narrowPeak format)  
    -b: second input file (bed or narrowPeak format)  
    -f 0.25: minimum overlap required as a fraction of A  
    -r: require reciprocal overlap  
    -wa: write the original entry in A for each overlap  

    To get all the columns from both files, we need to run bedtools intersect twice, swapping -a and -b


!!! Code 
    Merge the intersected peaks from both conditions into a final consensus peak set, allowing a maximum distance of 10bp between peaks to be merged (-d 10)

    ```{bash}
    cat results/04_consensus_peaks/Kidney_intersect.bed results/04_consensus_peaks/Cerebrum_intersect.bed | sort -k1,1 -k2,2n | bedtools merge -d 10 -i - > results/04_consensus_peaks/consensus_peaks.bed    ```


Add the `results/04_consensus_peaks/consensus_peaks.bed` track to IGV and have a look to the resulting peak annotation. 