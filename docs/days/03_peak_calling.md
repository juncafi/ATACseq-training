## 1. Subset .bams to keep only NF fragments

Before doing peak calling with MACS3, we will split fragments into nucleosome-free (NF) and nucleosome-associated fragments. 
For that, we will filter the BAM files based on fragment length (insert size).
We will then call the peaks based on NF reads only, and compare with peaks called on all reads (?)
Finally, we will also try peak calling with HMMRATAC, which is a peak caller specifically designed for ATAC-seq data and uses a Hidden Markov Model to identify accessible chromatin regions. (?)

**Task 1**

- Create an output folder named: `results/01_filtered_bams/NF_bams`
- Using `samtool view` keep only fragments with length: 1 - 100 bp. Keep output format as .bam file, and save them in `results/01_filtered_bams/NF_bams`. Follow the naming: ${sample_name}_NF.bam
- While doing that, keep the header of .bam file 

!!! note 
    In SAM format, column 9 is described as "TLEN: signed observed Template LENgth", which corresponds to the insert size 

<details>
<summary>Hint</summary>
You can first use samtools view, and then use awk commands to filter the reads for those which column 9 is: ($9>= 1 && $9<=100) || ($9<=-1 && $9>=-100)
</details>


??? success "Solution"

    ```{bash}
    # create new directory for peaks
    mkdir -p results/01_filtered_bams/NF_bams

    # filter for NF fragments only

    for bam in results/01_filtered_bams/*qc_bl_filt.sorted.bam; do
        sample_name=$(basename "$bam" .qc_bl_filt.sorted.bam)
        echo "Processing sample: $sample_name"
        samtools view -h $bam | awk 'substr($0,1,1)=="@" || ($9>= 1 && $9<=100) || ($9<=-1 && $9>=-100)' | \
        samtools view -b > results/01_filtered_bams/NF_bams/${sample_name}_NF.bam 
    done 

    ```
    Parameter explanation:  
     -h: keep header; -b: output in bam format; $9: insert size (TLEN field in SAM format); keep fragments with insert size between 1 and 100bp


**Task 2**


- Sort and index the bams files for next step
- Remove the unsorted file

??? success "Solution"

    ```{bash}
    for bam in results/01_filtered_bams/NF_bams/*_NF.bam; do
        sample_name=$(basename "$bam" _NF.bam)
        samtools sort -o results/01_filtered_bams/NF_bams/${sample_name}_NF.sorted.bam results/01_filtered_bams/NF_bams/${sample_name}_NF.bam
        rm results/01_filtered_bams/NF_bams/${sample_name}_NF.bam
        samtools index results/01_filtered_bams/NF_bams/${sample_name}_NF.sorted.bam
    done
    ```

## 2. Peak calling

### MACS3 based on NF fragments

Next we will call peaks using MACS3 on the NF reads only. This way, we focus on the fragments that are most likely to represent open chromatin regions only (nucleosome-free regions).

MACS3 is a widely used peak calling tool that can handle both narrow and broad peaks. We will start using the "callpeak" function of MACS3 to identify peaks in the NF reads.

**Task 3**

- Create a new folder named: `results/03_peak_calling`
- Create a subfolder inside called: `NF_peaks`
- Call peaks on NF reads using MACS3 callpeak function. Save the results inside r`esults/03_peak_calling/NF_peaks/` and name the files as: `NF_peaks_${sample_name}`


!!! MACS3 parameters
    You can have a look at the MACS3 documentation for more details on the parameters used [here](https://macs3-project.github.io/MACS/docs/callpeak.html)

<details>
<summary>Hint</summary>
    For ATAC-seq data, we will use the BAMPE format, which is suitable for paired-end data and we will set the genome size to "mm" for mouse. We will also set a q-value cutoff of 0.01 to control the false discovery rate.
</details>

??? success "Solution"

    ```{bash}
    mkdir -p results/03_peak_calling/NF_peaks

    samples=(PC_rep1 PC_rep2 RACM_rep1 RACM_rep2)
    samples=(Kidney_rep1 Kidney_rep2 Cerebrum_rep1 Cerebrum_rep2)

    # process NF reads
    bams_path="results/01_filtered_bams/NF_bams"

    for sample_name in "${samples[@]}"; do
        echo "Processing sample: $sample"
        macs3 callpeak -f BAMPE -t $bams_path/${sample_name}_NF.sorted.bam -g mm -q 0.01 --outdir results/03_peak_calling/NF_peaks/NF_peaks_${sample_name}/
    done    
    ```
    parameters explanation:  
    -f BAMPE: input file format is BAM paired-end  
    -t: input file (BAM file)  
    -g mm: It’s the mappable genome size or effective genome size (some are pre-computed, like mouse, and you can specify "mm"
    -q 0.01: q-value cutoff for peak detection  
    --outdir: output directory for peak files  

### MACS3 based on all fragments

We will do the same, but using all fragments in filtered .bam files 

**Task 4**

- Create a subfolder inside called: `all_peaks`
- Call peaks on all filtered reads using MACS3 callpeak function. Save the results inside r`esults/03_peak_calling/all_peaks/` and name the files as: `all_peaks_${sample_name}`

??? success "Solution"

    ```{bash}
    bams_path="results/01_filtered_bams/"
    mkdir -p results/03_peak_calling/all_peaks

    for sample in "${samples[@]}"; do
        echo "Processing sample: $sample"
        macs3 callpeak -f BAMPE -t $bams_path/${sample}.qc_bl_filt.sorted.bam -g mm -q 0.01 --outdir results/03_peak_calling/all_peaks/all_peaks_${sample}/
    done
    ```

