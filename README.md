<!--
---
title: Evolution of tissue-specific expression of ancestral genes across vertebrates and insects
output:
  html_document: default
  pdf_document: default
---
-->

# Evolution of tissue-specific expression of ancestral genes across vertebrates and insects

<figure>
  <img align="middle" src="https://github.com/fedemantica/bilaterian_GE/tree/main/docs/Fig1.png" width=600 height=400 />
  <figcaption>Dataset used for the analysis</figcaption>
<figure>

Content overview
-------

This is a github repository containing the code used to generate the results and images associated with the final version of the following publication: [Mantica et al, biorxive, 2022](). 
Something else about the general content.  

Table of contents
-------
* [Analysis](#analysis)  
  + [Preliminary annotation fixes](#preliminary-annotation-fixes)  
  + [Gene orthologies and sets](#gene-orthologies-and-sets)  
  + [Gene expression preprocessing](#gene-expression-preprocessing)  
  + [Gene expression quantification](#gene-expression-quantification)    
  + [GO transfers](#go-transfers)  
  + [Bilaterian ancestral modules](#bilaterian-ancestral-modules)  
  + [Bilaterian ancestral motif analyses](#bilaterian-ancestral-motif-analysis)
  + [Tissue-specificity call](#tissue-specificity-call)  
  + [Tissue-specificity inferences](#tissue-specificity-inferences)  
  + [Tissue-specificity inferences validation](#tissue-specificity-inferences-validation)  
  + [Supplementary analysis](#supplementary-analysis)
* [Figures](#figures)  
  + [Main Figures](#main-figures)  
    - [Figure 1](#figure-1)  
    - [Figure 2](#figure-2)  
    - [Figure 3](#figure-3)  
    - [Figure 4](#figure-4)  
    - [Figure 5](#figure-5)  
  + [Extended Data Figures](#extended-data-figures)  
    - [Extended Data Figure 1](#extended-data-figure-1)  
    - [Extended Data Figure 2](#extended-data-figure-2)  
    - [Extended Data Figure 3](#extended-data-figure-3)  
    - [Extended Data Figure 4](#extended-data-figure-4)  
    - [Extended Data Figure 5](#extended-data-figure-5)  
    - [Extended Data Figure 6](#extended-data-figure-6)  
    - [Extended Data Figure 7](#extended-data-figure-7)  
    - [Extended Data Figure 8](#extended-data-figure-8)  
    - [Extended Data Figure 9](#extended-data-figure-9)  
    - [Extended Data Figure 10](#extended-data-figure-10)  
  + [Supplementary figures](#supplementary-figures)
    - [Supplementary Figure 1](#supplementary-figure-1)  
    - [Supplementary Figure 2](#supplementary-figure-2)  
    - [Supplementary Figure 3](#supplementary-figure-3)  
    - [Supplementary Figure 4](#supplementary-figure-4)  
  

Analysis
------------

### Preliminary annotation fixes  
This folder contains a snakemake pipeline (i.e., **Snakefile**) used to fix annotation errors, together with the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the **bin** subfolder. The Snakefile includes rules to:  

* Generate protein fasta files for all species starting from the original annotation (best protein isoform per gene only).    
* Perform a preliminary gene orthology run with [Broccoli](), which will be used to identify problematic gene annotations.   
* Identify and fix broken genes, as described in XXX.   
* Identify and fix chimeric genes, as described in XXX.   
* Generate new gene annotations for all species, which include the corrections performed on broken and chimeric genes.  
* Run [BUSCO]() to evaluate the quality of the new gene annotations.  
* Correct the existing gene orthologies including the changes introduced for the broken/chimeric genes.  
* Merge the gene orthogroups based on known human onhologs.  

### Gene orthologies and sets  
This folder contains a snakemake pipeline (i.e., **Snakefile**) used to generate the different sets of orthogroups used throughout the paper. The folder also contains the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the **bin** subfolder. The Snakefile includes rules to:  

* Filter the original orthogroups based on the total number of genes and proportions per species.  
* Filter the original orthogroups based on the minimum expression.  
* Isolate bilaterian, vertebrate and insect conserved genes. Only bilaterian-conserved genes have been used for the analysis included in this paper.  
* Compute pairwise sequence similarities and pearson expression correlations between all possible pairs of genes (from different species) within an orthogroup. These values will be used to select the best-ancestral orthologs.  
* Select the best-ancestral orthologs for each species within the bilaterian-conserved orthogroups.  

### Gene expression preprocessing  
This folder contains a snakemake pipeline (i.e., **Snakefile**) used to generate gene expression quantification at the sample level starting from the fastq files (both for the public available and in-house generated data). The folder also contains the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the **bin** subfolder. The Snakefile includes rules to:  

* Download the public available fastq files.  
* Run quality control on all fastqs (public available and in-house).  
* Get the transcriptome index fasta needed for quantification of gene expression by [Kallisto]().    
* Quantify counts at the transcript level through [Kallisto]() for all the samples.  
* Merge counts at the gene level. 
* Normalize raw counts with DESeq2. 
* Generate gene expression quantification in logged (log2) and unlogged TPMs or RPKMs, using as input the DESeq2 normalized counts. Only the logged TPM quantification has been considered for the final results presented in the paper.  

### Gene expression quantification  
This folder contains a snakemake pipeline (i.e., **Snakefile**) used to generate gene expression quantification at the sample, meta-sample and tissue level starting from normalized counts across species. The folder also contains the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the **bin** subfolder. The Snakefile includes rules to:  

* Compute gene expression quantification at the meta-sample level, which is represented by the median expression across all samples included in a given meta-sample.  
* Compute gene expression quantification at the tissue level, which is represented by the average expression across all meta-samples of a given tissue.
* Generate gene expression matrices for the best-ancestral orthogroups (see [Gene orthologies and sets](#gene-orthologies-and-sets)) with orthogroup IDs on row, metasamples/tissues from all species on columns, and expression of the relative best-ancestral orthogroup as value.  
* Quantile-normalize the best-ancestral orthogroup gene expression matrices.  
* Perform an SVA normalization by species (from [Fukushima et al, NatCom, 2021]()), which I used for one of the controls.  

### GO transfers  
This folder contains a snakemake pipeline (i.e., **Snakefile**) used to generate GO annotations at the orthogroup level starting from the annotation of selected reference species (e.g., human and drosophila). The folder also contains the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the **bin** subfolder. The Snakefile includes rules to:  

* Generate annotations at the orthogroup level, both in text and gmt format. These will be used for GO enrichments throughout the paper.  
* Generate backgrounds to be used for GO enrichments.  

### Bilaterian ancestral modules  
This folder contains a snakemake pipeline (i.e., **Snakefile**) used to derive ancestral bilaterian modules of conserved tissue-specific genes. The folder also contains the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the **bin** subfolder. The Snakefile includes rules to:  

* Run the sparse least square discriminant analysis (spls-DA) on all tissues.  
* Select loadings (i.e., genes belonging to ancestral tissue-specific modules) and perform GO enrichments.  
* Evaluate association with Neural- and Testis-related phenotypes for orthogroups belonging to the Neural- and Testis-specific ancestral modules.  

### Bilaterian ancestral motif analyses  
This folder contains a snakemake pipeline (i.e., **Snakefile**) used to generate motif annotations and subsequent motif enrichment analysis for ancestral bilaterian modules. The folder also contains the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the **bin** subfolder. The Snakefile includes rules to:  

* <span style="color: red;">This is Luisito's code. I need to add it.</span>  

### Tissue-specificity call  
This folder contains a snakemake pipeline (i.e., **Snakefile**) used to perform tissue-specificity call across species. The folder also contains the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the **bin** subfolder. The Snakefile includes rules to:  

* Compute tau for all protein-coding genes in each species.  
* Compute relative expression across tissues for each protein-coding gene. This information will be used to identify the associated tissue and top tissue(s) for each protein-coding gene.  
* Identify the associated tissue and top tissue(s) for each protein-coding gene.  

### Tissue-specificity inference  
This folder contains a snakemake pipeline (i.e., **Snakefile**) used to perform inferences of tissue-specificity gains and losses for all tissues and throughout the phylogeny. The folder also contains the relative configuration file (**config.yml**) and a file with the parameters needed for cluster job submission (**cluster.json**). All relevant scripts can be found in the **bin** subfolder. The Snakefile includes rules to:  

* Select the best-TS orthogroups for each tissue.  
* Perform inferences of tissue-specificity gains for each of the tissues.  
* Perform inferences of tissue-specificity losses for each of the tissues.  

### Tissue-specificity inferences validation  
This folder contains an Rmarkdown file used to validate tissue-specificity gains with an orthogonal method.  


### Supplementary analysis  
These are all the analysis that were required during the revision phase of the paper. These analysis include:  

* Do something  
* Do something else  

Figures
------------

### Main Figures  

#### Figure 1  

<figure>
  <img align="middle" src="https://github.com/fedemantica/bilaterian_GE/tree/main/docs/Fig1.png" width=60 height=40 />
  <figcaption>Dataset used for the analysis</figcaption>
<figure>

#### Figure 2  

<figure>
  <img align="middle" src="https://github.com/fedemantica/bilaterian_GE/tree/main/docs/Fig2.png" width=600 height=400 />
  <figcaption>Dataset used for the analysis</figcaption>
<figure>

#### Figure 3  

<figure>
  <img align="middle" src="https://github.com/fedemantica/bilaterian_GE/tree/main/docs/Fig3.png" width=600 height=400 />
  <figcaption>Dataset used for the analysis</figcaption>
<figure>

#### Figure 4  

<figure>
  <img align="middle" src="https://github.com/fedemantica/bilaterian_GE/tree/main/docs/Fig4.png" width=600 height=400 />
  <figcaption>Dataset used for the analysis</figcaption>
<figure>

#### Figure 5  

<figure>
  <img align="middle" src="https://github.com/fedemantica/bilaterian_GE/tree/main/docs/Fig5.png" width=600 height=400 />
  <figcaption>Dataset used for the analysis</figcaption>
<figure>


### Extended Data Figures  

#### Extended Data Figure 1  
#### Extended Data Figure 2  
#### Extended Data Figure 3   
#### Extended Data Figure 4   
#### Extended Data Figure 5   
#### Extended Data Figure 6   
#### Extended Data Figure 7   
#### Extended Data Figure 8   
#### Extended Data Figure 9   
#### Extended Data Figure 10   


### Supplementary Figures    

#### Supplementary Figure 1  
#### Supplementary Figure 2  
#### Supplementary Figure 3  
#### Supplementary Figure 4  