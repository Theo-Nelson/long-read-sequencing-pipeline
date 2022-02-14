# A Long Read Sequencing Pipeline in the Cloud  

A Cloud Pipeline to Analyze Long Read Sequencing Data from Oxford Nanopore and PacBio Sequencers. This pipeline was originally implemented and has been extensively tested in Google Colaboratory.

---

### Content

  - [Installation](#installation)
  - [Features](#features)
  - [General Usage](#general-usage)
  - [Advanced options](#advanced-options)
  - [Details on the output](#details-on-the-output)
  - [Complementary functions](#complementary-functions)
  - [Citation](#citation)
  - [Contributors](#contributors)

---

### Installation

It is recommended that users install the pipeline in their Google Drive by following the instructions in this short Colab notebook (requires a Google Account, which is available for free): https://colab.research.google.com/drive/1CeGSw-tFIPaiXbELoTEvcraIfoeS646x?usp=sharing. If you wish to install the pipeline on your local machine you can run the following command (requires Git: https://github.com/git-guides/install-git). You will also need to procure your own reference genome. 
```bash
github clone https://github.com/Theo-Nelson/long-read-sequencing-pipeline
```
Once this is complete you can open the file ```long_read_rna_seq_analysis_prebuilt_indices.ipynb``` either in Google Colaboratory or on your local machine and begin working with the pipeline.

---

### Features

Links to documentation for each software package can be found in the notebook file. The features are colored by importance: red programs are critical for successful generation of aligned files, blue are optional but additive, while green denotes software suitable for advanced users. 

#### Parameter Input and User Instructions ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+)
This section allows you to setup and automate the remainder of the analysis. Most parameters are necessary for specific programs, but four specific variables: ```PIPELINE_FILE_PATH```, ```ACC```, ```INDEX_FILE_PATH``` and ```ANNOTATION_FILE_PATH``` are critical. Please note that ```PARTITIONED_INDEX_FILE_PATH``` takes the place of ```INDEX_FILE_PATH``` if you are using minimap2 featherweight alignment (described below). The pipeline file path describes where in your Google Drive file system or local file system the pipeline is set up. The ACC variable specifies either the FASTQ file on your machine containing long-read sequencing data or the run accession number (e.g. [SRR12389274](https://www.ebi.ac.uk/ena/browser/view/SRR12389274)) for the publically available file which you wish to analyze. You can search for samples relevant to your field of interest with the European Nucleotide Archive's Advanced Search Feature: https://www.ebi.ac.uk/ena/browser/advanced-search (tutorial: https://www.youtube.com/watch?v=ugLaYRgh1pE). The index file path describes where the reference genome exists within your file system; the annotation file path describes where the reference annotation exists in the same file system. 

#### Mounting your Google Drive ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+)
This will connect your current hardware instance to your Google Drive and allow you to permeantly store your analysis. Most academic users enjoy unlimited Google Drive Storage space, while basic users are given 30 GB of free storage space.

#### Managing Software via BioConda ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+)
BioConda is a centralized package manager necessary to install the remaining programs. 

#### Kingfisher: fast and flexible program for procurement of sequence files ![#c5f015](https://via.placeholder.com/15/c5f015/000000?text=+)
This program will download sequence files from the European Nucleotide Archive. 

#### FastQC: A quality control tool for high throughput sequence data ![#c5f015](https://via.placeholder.com/15/c5f015/000000?text=+)
This tool will generate basic high-level statistics regarding read length and sequencing quality. Long-read sequences generally have poor basepair quality. 

#### minimap2: A versatile pairwise aligner for genomic and spliced nucleotide sequences (Google Colab Pro required) ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+)
Alignment will occur either with minimap2 (available in Google Colab Pro due to memory constraints) or minimap2 featherweight alignment (available in the free version of Google Colab). minimap2 is a versatile aligner which maps reads onto a reference genome. 

#### minimap2 featherweight alignment ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+)
Alignment will occur either with minimap2 (available in Google Colab Pro due to memory constraints) or minimap2 featherweight alignment (available in the free version of Google Colab). minimap2 is a versatile aligner which maps reads onto a reference genome. The featherweight version saves memory by splitting the reference genome into four parts, aligning to each of these in turn and thereafter merging the results. 

#### samtools: Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format ![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+)
This tool stores your aligned reads in a compressed binary format, sorting and indexing them along the way for quick access by chromosomal location. These files could be downloaded and viewed in a program such as the Integrated Genomics Viewer available from the Broad Institute (https://software.broadinstitute.org/software/igv/). 

#### TranscriptClean: correct mismatches, microindels, and noncanonical splice junctions ![#c5f015](https://via.placeholder.com/15/c5f015/000000?text=+)
This feature will try to polish your reads to adhere to known biological principles and decrease variability among reads aligning to the same region.  

#### featureCounts: an efficient general purpose program for assigning sequence reads to genomic features ![#c5f015](https://via.placeholder.com/15/c5f015/000000?text=+)
This tool will produce a count matrix assigning reads to features within your reference annotation. 

#### TAMA: Transcriptome Annotation by Modular Algorithms (Google Colab Pro required) ![#c5f015](https://via.placeholder.com/15/c5f015/000000?text=+)
This memory-intensive program will produce a reference transcriptome (GTF file) based on your long-read sequencing data. 

#### svist4get: a simple visualization tool for genomic tracks from sequencing experiments ![#c5f015](https://via.placeholder.com/15/c5f015/000000?text=+)
This tool will generate a graph of read coverage (i.e. how many reads align to a given region) for a specific chromosomal region (the names for these must match the names of chromosome within your reference genome/annotation). 

#### AlignQC: Long read alignment analysis ![#c5f015](https://via.placeholder.com/15/c5f015/000000?text=+)
This time-intensive tool will generate statistics regarding the aligned long-read sequences, including the number of novel transcripts within your sample.  

#### MakeHub: Fully automated generation of UCSC assembly hubs ![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+)
This feature will create a track hub which you can host on a public file service allowing for byte-access. This can then be connected to the UCSC genome browser and viewed as a track (e.g.: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A116533435%2D116536513&hgsid=1282820889_LUJAMUR9DzBvxiKtV6M3Qhk7c0Iv). Please see the MakeHub documentation for more details: https://github.com/Gaius-Augustus/MakeHub#how-to-use-makehub-output-with-ucsc-genome-browser

#### MultiQC: Aggregate results from bioinformatics analyses across many samples into a single report ![#c5f015](https://via.placeholder.com/15/c5f015/000000?text=+)
This tool will collect results from other quality control tools into a viewable report.

---

### General Usage 

The default mode to run ***bambu** is using a set of aligned reads (bam files), reference genome annotations (gtf file, TxDb object, or bambuAnnotation object), and reference genome sequence (fasta file or BSgenome). ***bambu*** will return a summarizedExperiment object with the genomic coordinates for annotated and new transcripts and transcript expression estimates.

We highly recommend to use the same annotations that were used for genome alignment. If you have a gtf file and fasta file you can run ***bambu*** with the following options:

```rscript
test.bam <- system.file("extdata", "SGNex_A549_directRNA_replicate5_run1_chr9_1_1000000.bam", package = "bambu")
  
fa.file <- system.file("extdata", "Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa", package = "bambu")

gtf.file <- system.file("extdata", "Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf", package = "bambu")

bambuAnnotations <- prepareAnnotations(gtf.file)

se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file)

```
**Transcript discovery only (no quantification)**

```rscript
bambu(reads = test.bam, annotations = txdb, genome = fa.file, quant = FALSE)
```

**Quantification of annotated transcripts and genes only (no transcript/gene discovery)**

```rscript
bambu(reads = test.bam, annotations = txdb, genome = fa.file, discovery = FALSE)
```

**Large sample number/ limited memory**     
For larger sample numbers we recommend to write the processed data to a file:

```rscript
bambu(reads = test.bam, rcOutDir = "./bambu/", annotations = bambuAnnotations, genome = fa.file)
```

For very large samples (>100 million reads) where memory is limiting we recommend running Bambu in lowMemory mode:

```rscript
bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file, lowMemory = TRUE)
```

---


### Use precalculated annotation objects

You can also use precalculated annotations.

If you plan to run ***bambu*** more frequently, we recommend to save the bambuAnnotations object.

The bambuAnnotation object can be calculated from a *.gtf* file:

```rscript
annotations <- prepareAnnotation(gtf.file)
```

From *TxDb* object

```rscript
annotations <- prepareAnnotations(txdb)
```

---

### Advanced Options

**More stringent filtering thresholds imposed on potential novel transcripts**    
 
- Keep novel transcripts with min 5 read count in at least 1 sample: 

```rscript
bambu(reads, annotations, genome, opt.discovery = list(min.readCount = 5))
```

- Keep novel transcripts with min 5 samples having at least 2 counts:

```rscript
bambu(reads, annotations, genome, opt.discovery = list(min.sampleNumber = 5))
```

- Filter out transcripts with relative abundance within gene lower than 10%: 

```rscript
bambu(reads, annotations, genome, opt.discovery = list(min.readFractionByGene = 0.1))
```

- Set novel transcript discovery rate to 50% of the detected transcripts (lower is more): 

```rscript
bambu(reads, annotations, genome, NDR = 0.5)
```

**Quantification without bias correction**     

 The default estimation automatically does bias correction for expression estimates. However, you can choose to perform the quantification without bias correction.

```rscript
bambu(reads, annotations, genome, opt.em = list(bias = FALSE))
```

**Parallel computation**      
 ***bambu***  allows parallel computation.  

```rscript
bambu(reads, annotations, genome, ncore = 8)
```

See [our page](https://goekelab.github.io/bambu/) for a complete step-by-step workflow and manual on how to customize other condictions.

---

### Details on the output 

***bambu*** will output different results depending on whether *quant* mode is on. 

By default, *quant* is set to TRUE, 
so ***bambu*** will generate a *SummarizedExperiment* object that contains the transcript expression estimates.  

* access transcript expression estimates by ***counts()***, including a list of variables: counts, CPM, fullLengthCount, partialLengthCounts, and uniqueCounts, and theta
    + counts: expression estimates
    + CPM: sequencing depth normalized estimates
    + fullLengthCounts: estimates of read counts mapped as full length reads for each transcript
    + partialLengthCounts: estimates of read counts mapped as partial length reads for each transcript
    + uniqueCounts: counts of reads that are uniquely mapped to each transcript
    + theta: raw estimates
* access annotations that are matched to the transcript expression estimates by ***rowRanges()***
* access transcript to gene id map by ***rowData()***, *eqClass* that defines the equivalent class transcripts is also reported

In the case when *quant* is set to FALSE, i.e., only transcript discovery is performed, 
***bambu*** will report the *grangeslist* of the extended annotations

### Complementary functions

**Transcript expression to gene expression**

```rscript
transcriptToGeneExpression(se)
```

**Visualization**

 You can visualize the novel genes/transcripts using ***plotBambu*** function 

```rscript
plotBambu(se, type = "annotation", gene_id)

plotBambu(se, type = "annotation", transcript_id)
```

- ***plotBambu*** can also be used to visualize the clustering of input samples on gene/transcript expressions

```rscript
plotBambu(se, type = "heatmap") # heatmap 

plotBambu(se, type = "pca") # PCA visualization
```

- ***plotBambu*** can also be used to visualize the clustering of input samples on gene/transcript expressions with grouping variable

```rscript
plotBambu(se, type = "heatmap", group.var) # heatmap 

plotBambu(se, type = "pca", group.var) # PCA visualization
```

**Write bambu outputs to files**

- ***writeBambuOutput*** will generate three files, including a *.gtf* file for the extended annotations, and two *.txt* files for the expression counts at transcript and gene levels.

```rscript
writeBambuOutput(se, path = "./bambu/")
```

---

### Citation

to be added.

---

### Contributors

This pipeline was developed by [Theodore Nelson](https://github.com/Theo-Nelson) at the Columbia University Irving Medical Center. 
