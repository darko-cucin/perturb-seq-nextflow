# perturb-seq pipeline

## Introduction

**Perturb-seq** is an experimental approach to map the transcriptional effects of genetic perturbations by combining CRISPR-based genetic screening with information-rich, single-cell RNA-sequencing phenotypes.

**Perturb-seq** accurately identifies individual gene targets, gene signatures, and cell states affected by individual perturbations and their genetic interactions.

Applied at the genome-scale, **perturb-seq** enables systematic prediction of gene function and in-depth study of complex cellular phenotypes such as aneuploidy, differentiation, RNA processing, and stress-specific regulation of the mitochondrial genome.

## Getting started - installing dependencies

To run **perturb-seq** pipeline, you should install:

* 1. **Nextflow** - as this pipeline is written in **Nextflow** language it has to be run with the same one. Requirements for installing and running **Nextflow**:

  * **Bash 3.2** (or later); 

  * **Java 17** (or later, up to 25);
 
    * If you find that your **Java** version is below v17 after installation, you may encounter difficulties running **Nextflow** since a lower **Java** version is set as the default. To resolve this, execute the following command:

      * ```export JAVA_HOME=$(/usr/libexec/java_home -v 17) ```
 
Installing steps are as follows:
  * Download the executable package by copying and pasting one of the following commands: 

    * ```wget -qO- https://get.nextflow.io | bash  ```
    * ```curl -s https://get.nextflow.io | bash  ```

  * Make the binary executable:
    
    * ```chmod +x nextflow```

    An executable script for **Nextflow** can be installed in the directory which is not accessible by your ```$PATH```. When **Nextflow** is installed, the following message is printed:

``` the executable file Nextflow has been created in the folder: /path/to/folder. ```

To make **Nextflow** executable from the terminal either move the **Nextflow** to a location that is in your ```$PATH```, or add the file location to the ```$PATH``` using the following command:

```PATH=$PATH:"path/to/folder"```


* 2. **Docker** - **Nextflow** language can be run with different container runtimes. We decided to use **Docker** to run all modules/tools that are contained within **perturb-seq** pipeline. Consequently, there is no need to install tools locally and installing local dependencies. To install **Docker**, you can refer to the following [page](https://docs.docker.com/engine/install/).
 

## Running a pipeline

If you want to run **Perturb-seq** pipeline, the following input files are required:

* 1. **Samplesheet file** (```--input```) - you will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row.

Keep in mind that **perturb-seq** analysis requires read files from single-cell and CRISPR assays/experiments. This type of analysis, where single-cell gene expression data is coupled with other data, which gives an additional layer of information, is called feature barcode analysis. For the **perturb-seq** data feature barcode analysis is performed by coupling results of cell pertubrations affected by CRISPR experiments with single-cell RNA-seq. This information will be important for creating **Samplesheet file**.

* 2. **Index for GEX libraries**(```--index_reference_gex```) or **GEX Reference genome fasta** (```--fasta_reference_gex```) and **GEX Annotation file** (```--gtf_reference_gex```) - path of folder containing 10x-compatible genome reference  for **STAR Solo** or Path to FASTA file containing your genome reference for **GEX Reference genome fasta** and path to input genes GTF file for **GEX Annotation file** input.

 * 3. **Index for CRISPR libraries**(```--index_reference_crispr```) or **CRISPR Reference genome fasta** (```--fasta_reference_crispr```) and **CRISPR Annotation file** (```--gtf_reference_crispr```) - path of folder containing genome reference for CRISPR reference in **STAR Solo** or Path to FASTA file containing your CRISPR genome reference for **GEX Reference genome fasta** and path to input CRISPR genes GTF file for **GEX Annotation file** input.

* **Important note** References for different libraries (GEX and CRISPR) have to be separated. For GEX library type, you provide a standard reference genome of the species. For the CRISPR library you need to make your own reference FASTA and GTF files. CRISPR reference fasta should contains sequence of guide RNA / spacer and the scaffold sequence, which is a constant. Every sequence is separated contig. GTF file should contain information about these guide RNAs. Examples of these files can be found in assets/test_data directory.

* 4. **Barcode whitelist file**(```--barcode_whitelist```) or **Version of 10x chemistry** (```--tenx_version```) and **Annotation file** (```--gtf```) - File of valid cell barcodes for 10x experiments. If this file is not specified **Version of 10x chemistry** (```--tenx_version```) can be specified to get whitelist (from **10XV1** to **10XV4**).
### Multiple runs of the same sample

In some cases, you will run the same sample on multiple flow cells. In this case the samplesheet file has to have the following structure:

```samplesheet.csv
sample,fastq_1,fastq_2,lib_type
SC3_v3,Directory/SC3_v3_gex/SC3_v3_gex_S5_L001_R1_001.fastq.gz,Directory/SC3_v3_gex/SC3_v3_gex_S5_L001_R2_001.fastq.gz,GEX
SC3_v3,Directory/SC3_v3_gex/SC3_v3_gex_S5_L001_I1_001.fastq.gz,Directory/SC3_v3_gex/SC3_v3_gex_S5_L001_I2_001.fastq.gz,GEX
SC3_v3,Directory/SC3_v3_gex/SC3_v3_gex_S5_L002_R1_001.fastq.gz,Directory/SC3_v3_gex/SC3_v3_gex_S5_L002_R2_001.fastq.gz,GEX
SC3_v3,Directory/SC3_v3_gex/SC3_v3_gex_S5_L002_I1_001.fastq.gz,Directory/SC3_v3_gex/SC3_v3_gex_S5_L002_I2_001.fastq.gz,GEX
SC3_v3,Directory/SC3_v3_crispr/SC3_v3_crispr_S4_L001_R1_001.fastq.gz,Directory/SC3_v3_crispr/SC3_v3_crispr_S4_L001_R2_001.fastq.gz,CRISPR
SC3_v3,Directory/SC3_v3_crispr/SC3_v3_crispr_S4_L001_I1_001.fastq.gz,Directory/SC3_v3_crispr/SC3_v3_crispr_S4_L001_I2_001.fastq.gz,CRISPR
SC3_v3,Directory/SC3_v3_crispr/SC3_v3_crispr_S4_L002_R1_001.fastq.gz,Directory/SC3_v3_crispr/SC3_v3_crispr_S4_L002_R2_001.fastq.gz,CRISPR
SC3_v3,Directory/SC3_v3_crispr/SC3_v3_crispr_S4_L002_I1_001.fastq.gz,Directory/SC3_v3_crispr/SC3_v3_crispr_S4_L002_I2_001.fastq.gz,CRISPR
```

| Column         | Description                                                                                                                                                                            |
|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample`       | Sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. It is important to mention that sample name has to be the same for both single-cell read files and CRISPR reads files which will be subsequently separated in the pipeline. The sample name has to have the following name convention: `file_name - ${gex/crispr}_S{i}_L00{i}_R/L{i}_001.fastq.gz` |
| `fastq_1`      | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                              |
| `fastq_2`      | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `lib_type`      | Library type for a given pair of fastq files. (GEX or CRISPR)                                                             |


###  Running multiple samples

In some cases, more single-cell/CRISPR samples will be run, which will require parallelization of execution. In this use case, the samplesheet will have the following format:

```samplesheet.csv
sample,fastq_1,fastq_2,lib_type
SC3_v3,Directory/SC3_v3_gex/SC3_v3_gex_S5_L001_R1_001.fastq.gz,Directory/SC3_v3_gex/SC3_v3_gex_S5_L001_R2_001.fastq.gz,GEX
SC3_v3,Directory/SC3_v3_crispr/SC3_v3_crispr_S4_L001_R1_001.fastq.gz,Directory/SC3_v3_crispr/SC3_v3_crispr_S4_L001_R2_001.fastq.gz,GEX
SC3_v4,Directory/SC3_v4_gex/SC3_v4_gex_S5_L001_I1_001.fastq.gz,Directory/SC3_v3_gex/SC3_v4_gex_S5_L001_I2_001.fastq.gz,CRISPR
SC3_v3,Directory/SC3_v4_crispr/SC3_v4_crispr_S4_L001_R1_001.fastq.gz,Directory/SC3_v3_crispr/SC3_v4_crispr_S4_L001_R2_001.fastq.gz,CRISPR
```
  
Example of CRISPR Fasta reference file 

 ```
>Non-Targeting-5
ACTCGAAATCACCTATGGTAGTTTAAGAGCTAAGCTGGAA
>Non-Targeting-7
TTATGTGAGCACGCCATTACGTTTAAGAGCTAAGCTGGAA
>Non-Targeting-8
CGACGGTAATGCACCTACTAGTTTAAGAGCTAAGCTGGAA
>APH1A-1
GGCAACGCGACCCCACGAGGTTTAAGAGCTAAGCTGGAA
>APH1A-2
ATGTCACCCCCAGACCCCGGTTTAAGAGCTAAGCTGGAA
>CDKN3-1
TGCAGCGCCGGCGACTCACGTTTAAGAGCTAAGCTGGAA
>CDKN3-2
CGGGGCACCGGTGAGTCGCGTTTAAGAGCTAAGCTGGAA
>EZR-1
CACTCGGCGGACGCAAGGGGTTTAAGAGCTAAGCTGGAA
>EZR-2
GCGCACTCGGCGGACGCAAGTTTAAGAGCTAAGCTGGAA
>GRB2-1
TGCTGCTTCGGCGACCGGGGTTTAAGAGCTAAGCTGGAA
```
Example of CRISPR GTF file 

```
Non-Targeting-5	CRISPR	exon	1	40	.	+	.	gene_id "Non-Targeting-5"; gene_name "Non-Targeting"; gene_type "gRNA"
Non-Targeting-7	CRISPR	exon	1	40	.	+	.	gene_id "Non-Targeting-7"; gene_name "Non-Targeting"; gene_type "gRNA"
Non-Targeting-8	CRISPR	exon	1	40	.	+	.	gene_id "Non-Targeting-8"; gene_name "Non-Targeting"; gene_type "gRNA"
APH1A-1	CRISPR	exon	1	39	.	+	.	gene_id "APH1A-1"; gene_name "APH1A"; gene_type "gRNA"
APH1A-2	CRISPR	exon	1	39	.	+	.	gene_id "APH1A-2"; gene_name "APH1A"; gene_type "gRNA"
CDKN3-1	CRISPR	exon	1	39	.	+	.	gene_id "CDKN3-1"; gene_name "CDKN3"; gene_type "gRNA"
CDKN3-2	CRISPR	exon	1	39	.	+	.	gene_id "CDKN3-2"; gene_name "CDKN3"; gene_type "gRNA"
EZR-1	CRISPR	exon	1	39	.	+	.	gene_id "EZR-1"; gene_name "EZR"; gene_type "gRNA"
EZR-2	CRISPR	exon	1	39	.	+	.	gene_id "EZR-2"; gene_name "EZR"; gene_type "gRNA"
GRB2-1	CRISPR	exon	1	39	.	+	.	gene_id "GRB2-1"; gene_name "GRB2"; gene_type "gRNA"
GRB2-2	CRISPR	exon	1	39	.	+	.	gene_id "GRB2-2"; gene_name "GRB2"; gene_type "gRNA"
GSK3A-1	CRISPR	exon	1	39	.	+	.	gene_id "GSK3A-1"; gene_name "GSK3A"; gene_type "gRNA"
GSK3A-2	CRISPR	exon	1	39	.	+	.	gene_id "GSK3A-2"; gene_name "GSK3A"; gene_type "gRNA"
HRAS-1	CRISPR	exon	1	39	.	+	.	gene_id "HRAS-1"; gene_name "HRAS"; gene_type "gRNA"
 ```

 
Other parameters that can be specified are as follows (we will list the most important ones):

* **Start of Unique Molecular Identifier (UMI)** (```--umi_start```) - Start base of UMI.
* **Length of UMI** (```umi_len```) - UMI Length
* **Start of Cell Barcode** (```--cb_start```) - Cell barcode start base.
* **Length of Cell Barcode** (```--cb_len```) - Cell barcode length.
* **Number of matched bases CRISPR** (```--number_of_bases_to_match_crispr```) - alignment will be output only if the number of matched bases is higher
than or equal to this value. Important parameter for CRISPR library types, as one CRISPR read will never match the reference with full length (because reference contigs are shorter than a read)
* **Percentage of matched bases CRISPR** (```--percent_of_bases_to_match_crispr```) - alignment will be output only if the number of matched bases is higher
than or equal to this value. Important parameter for CRISPR library types, as one CRISPR read will never match the reference with full length (because reference contigs are shorter than a read).

The example command line for running a pipeline is as follows:

```nextflow run main.nf -profile docker --input assets/samplesheet.csv --index_reference_gex assets/test_data/chr20 --fasta_reference_crispr assets/test_data/crispr_reference.fa.gz --gtf_reference_crispr assets/test_data/crispr_reference.gtf.gz --outdir testing```

You can also run with a test profile ```nextflow run main.nf -profile,test```


## Pipeline outputs

This section describes the output produced by the pipeline.

### FastQC outputs

**FastQC** gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

### STAR Solo outputs

* Most of **STAR Solo** outputs are stored within directory with **.Solo.out** suffix. We will list a couple of files:
* **features.tsv.gz** file - The **features.tsv.gz** file serves as the annotation index for your single-cell dataset, mapping the row indices of the count matrix to specific genomic identifiers. It typically contains three columns (**Gene ID**, **Gene Name**, and **Feature Type**). It ensures that downstream analysis tools can correctly identify and label the expression data.
* **barcodes.tsv.gz** file - The **barcodes.tsv.gz** file is a compressed list containing all the unique cellular identifiers (sequences) that were successfully identified and quantified in the sample. Each line in this file corresponds to a single cell, acting as the "column names" for the gene expression matrix.
* **matrix.tsv.gz** file - The **matrix.mtx.gz** file is a compressed sparse matrix that stores the actual quantification values (UMIs or counts) for the experiment. It uses a coordinate system to map expression levels to specific intersections of the features (rows) and barcodes (columns), omitting zeros to significantly reduce file size.
* **Sorted Coordinated BAM files** - This file contains the raw sequence alignments for the experiment, sorted by genomic position to allow for rapid indexing and visualization.
* Other files are mostly log and metrics files that can be used to check the overall performance of STAR solo. 

