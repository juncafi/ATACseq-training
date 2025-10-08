
In order to assess the quality of the ATACseq experiment, there are several metrics that can be calculated from the aligned reads.

Important metrics include:
 - Fragment size distribution: The distribution of fragment sizes can indicate the quality of the library preparation. A good ATAC-seq library should have a characteristic pattern with peaks corresponding to nucleosome-free regions (NFRs) and mono-, di-, and tri-nucleosomes.
 - TSS enrichment: The enrichment of reads at transcription start sites (TSS) is a key indicator of data quality. High-quality ATAC-seq data should show a strong enrichment of reads at TSSs.

 `ATAQV` is a tool that calculates a variety of QC metrics for ATAC-seq data. It provides a comprehensive report that includes fragment size distribution, TSS enrichment, and other important metrics. 


**Task 1**

- Create a new folder for QC results called: `results/02_QC_post_aligment`
- Run ATAQV on each filtered .bam file 
- Using ATAQV tool, create a report from all ATAQV sample results


!!! ATAQV
    You can find further information about ATAQC tool and commands in the Usage section [here](https://github.com/ParkerLab/ataqv)
    Ataqv needs a Transcription Start Site (TSS) reference file to compute TSS enrichment score. You will find this file in: `data/references/ENCODE_mm10_M21_TSS_reference.bed`
    We have downloaded the TSS reference gile from ENCODE project (https://www.encodeproject.org/files/ENCFF498BEJ/)


<details>
<summary>Hint</summary>
   First run ataqv for each "*qc_bl_filt.sorted.bam" file, important parameters are: --tss-file, mouse, --metrics-file
   Second run mkarv with all "*json" files
</details>

??? success "Solution"
    ```{bash}
    mkdir -p results/02_QC_post_aligment
    TSS_bed="/data/references/ENCODE_mm10_M21_TSS_reference.bed"

    for bam in results/01_filtered_bams/*qc_bl_filt.sorted.bam; do
        echo "Processing file: $bam"
        sample_name=$(basename "$bam" .qc_bl_filt.sorted.bam) # extract sample name without path and extension
        echo $sample_name
        ataqv --name $sample_name --metrics-file results/02_QC_post_aligment/$sample_name.ataqv.json --tss-file $TSS_bed mouse $bam > results/02_QC_post_aligment/$sample_name.ataqv.out
    done

    cd results/02_QC_post_aligment
    mkarv summary_ataqv Kidney_rep1.ataqv.json Kidney_rep2.ataqv.json Cerebrum_rep1.ataqv.json Cerebrum_rep2.ataqv.json 
    cd ../../

    ```
    parametr explanation:  
     --name: sample name for output files  
     --metrics-file: output file for metrics in json format  
     --tss-file: bed file with TSS locations  
     mouse: genome (can be mouse or human)  
     $bam: input bam file  
    mkarv is a tool to create a summary report from multiple ataqv json files


**Task 2**

Open the QC report `results/02_QC_post_aligment/summary_ataqv/index.html` and have a look at the QC metrics. Do you think the experiment worked well? 
Is there some metrics that would concern you?

