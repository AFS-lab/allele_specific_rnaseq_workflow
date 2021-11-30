# Pipeline for Allele-Specific Transcript Quantification of RNA-seq

## Summary

This workflow will process raw Illumina paired-end reads for allele-specific transcript quantification to produce matrices of counts summarised in four ways: 

- at the gene level 
- at the isoform level 
- at the gene level with specific counts for each allele
- at the isoform level with specific counts for each allele

The output of the pipeline can then be used in downstream analysis, either for regular differential analysis (using _DESeq2_) or for allele-specific bias analysis (for example using ISoLDE). 
These downstream steps need to be customised depending on the dataset being used.
For an example of downstream analysis see [this repository](https://github.com/tavareshugo/mammary_cellsort_zfp57_rnaseq).

## Details

The pipeline is implemented using _snakemake_ and, briefly, performs the following steps: 

- quality-trimming and adapter removal from reads using _cutadapt_
- alignment to the reference genome using _STAR_ (this is only used for quality assessment)
- generate a "diploid" transcriptome based on the reference genome (FASTA), gene annotation (GTF) and variants (VCF) - see `workflow/scripts/diploidomes.sh` for details of how this is done
- transcript quantification using _salmon_ against the "diploid" transcriptome
- create four `DESeqDataSet` objects, with count matrix information summarised as detailed above
- several quality control metrics (e.g. _FastQC_ and _rseqc_) compiled to a _multiqc_ report

## Setup & Instructions 

- Install _snakemake_ using Conda: `conda create --name snakemake snakemake`
- Clone the repository: `git clone `
- Edit the `read_info.csv` file with information relevant to your data. Further details about this file given below.
- Download reference genome and transcript annotation to a directory named `data/external/reference`. An example script that downloads the reference genome for the mouse GRCm38 assembly is given in `scripts/download_reference.sh`.
- Prepare or download a VCF file with SNPs/Indels for the hybrid individual used in the experiment. We provide a script that downloads mouse SNPs for the CAST strain and creates a synthetic hybrid with the reference B6 strain (`scripts/download_variants.sh`).
- Edit the `workflow/config.yaml` file with any options relevant to your analysis.
- Run the workflow locally using: `snakemake --use-conda`. If you are running this on a HPC SLURM environment, see the example given in [this repository](https://github.com/tavareshugo/mammary_cellsort_zfp57_rnaseq). 


### Read files' information

The file `read_info.csv` should have at least four columns named as follows: 

- `id` - unique identifier for a sequencing sample
- `sample` - identifier for a biological sample; samples with the same name will be merged together
- `library` - identifier for the library; samples with the same library ID are de-duplicated together
- `fq1` - path to read 1 FastQ file - can be compressed
- `fq2` - path to read 2 FastQ file - can be compressed

Other columns can be included, for example information about genotypes, treatments, timepoints, etc. 

If the same library was re-sequenced, the files will be merged together. 
If the same RNA sample was used to make a new library and then resequenced, this will also be merged together (but PCR duplicates will be correctly assessed at the library level). 

An example is included in the `read_info.csv` file given in this repository. 
In this example we have 3 samples: sample 1 has one library that was sequenced twice; sample 2 was only sequenced once; sample 3 was also sequenced twice, but from two different libraries.


**Note about duplicates |**
Note that, despite there being a step assessing the level of duplication (using `picard MarkDuplicates`), the duplicates are in fact included in the transcript quantification using Salmon. 
The de-duplication metrics are mostly for quality control (and will be included in the MultiQC report).


## Limitations of this pipeline

At the moment the pipeline does not include the downstream steps to run _ISoLDE_ and _DESeq2_. 
This is mainly because these steps are somewhat specific to each experiment (depending on what conditons/groups of samples there are). 
But the output files can be readily used with those tools. 

An example of downstream analysis is given in the [accompanying repository on mammary gland data](https://github.com/tavareshugo/mammary_cellsort_zfp57_rnaseq). 