# L-RAP: A Cloud-Computing Pipeline for the Analysis of Long-Read RNA Sequencing Data    

A Cloud Pipeline to Analyze Long Read Sequencing Data from Oxford Nanopore and PacBio Sequencers, originally implemented and tested in Google Colaboratory.

---

### Content

  - [Installation](#installation)
  - [Features](#features)
  - [General Usage](#general-usage)
  - [Details on the output](#details-on-the-output)
  - [Troubleshooting Guide](#troubleshooting-guide)
  - [Citation](#citation)
  - [Contributors](#contributors)

---

### Installation

#### Google Colaboratory 

Open the following notebook (requires a Google Account): https://colab.research.google.com/drive/1OORV8X1Qe3RP7u-bIaJh0PCsEF7Y_imJ (that's it!)

#### Reference Genomes and Annotations

Within Google Colaboratory, reference genomes and annotations can be downloaded for hg38 (human) and mm39 (mouse) from UCSC. For additional species, users should input their own download links. Users are expected to provide their own reference genomes if running the pipeline on their local machine. Ensembl reference annotations for a large number of species can be accessed at the following address: http://ftp.ensembl.org/pub/

#### For Users with Unlimited Google Drive Storage

It is recommended that users install the pipeline and reference genome directly into their Google Drive by following the instructions in this short Colab notebook: https://colab.research.google.com/drive/1CeGSw-tFIPaiXbELoTEvcraIfoeS646x?usp=sharing. Once installation is complete, users should open ```long_read_rna_seq_analysis.ipynb``` to begin working with the pipeline. For first-time users Colaboratory can be installed as follows in Google Drive:  ```NEW``` => ```MORE``` => ```+ Connect more apps``` => Search ```Colaboratory```. Please make sure that the ```$PIPELINE_FILE_PATH``` variable matches between your installation and pipeline script. By default, this is set to ```/content/drive/MyDrive```.

#### Local Installations

The pipeline can be installed on a local machine with the following command (requires Git: https://github.com/git-guides/install-git).
```bash
github clone https://github.com/Theo-Nelson/long-read-sequencing-pipeline
```

Please note that the installation scripts within the pipeline have not been tested outside of the Google Colaboratory environment, which is generally set up in the following manner, according to ```cat /etc/os-release```. 
```
NAME="Ubuntu"
VERSION="18.04.6 LTS (Bionic Beaver)"
ID=ubuntu
ID_LIKE=debian
PRETTY_NAME="Ubuntu 18.04.6 LTS"
VERSION_ID="18.04"
HOME_URL="https://www.ubuntu.com/"
SUPPORT_URL="https://help.ubuntu.com/"
BUG_REPORT_URL="https://bugs.launchpad.net/ubuntu/"
PRIVACY_POLICY_URL="https://www.ubuntu.com/legal/terms-and-policies/privacy-policy"
VERSION_CODENAME=bionic
UBUNTU_CODENAME=bionic
```
#### Optional Upgrades

No pipeline program requires these optional upgrades, which allow for speedier runtimes and increased cloud storage:
- [Google Colab Pro](https://colab.research.google.com/signup), a service which costs $10/month, provides up to 25 GB of RAM. Once purchased, this feature can be activated by selecting ```Runtime``` from the top menu => ```Change Runtime Type``` => ```Runtime shape``` => ```High RAM```.
- [Google One](https://one.google.com/u/2/storage), a service which costs $2.99/month for 200 GB of storage. Once purchased, this feature will be automatically applied. Please note that users with unlimited storage through academic or workspace accounts do not require this service. Once purchased, please see the section above entitled, ```For Users with Unlimited Google Drive Storage```.

---

### Features

Links to the official documentation for each software package can be found in ```long_read_rna_seq_analysis.ipynb```. Here, we describe the general function of the programs and color them by importance: red programs (![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/red_square.png)) are critical for successful generation of aligned reads, pink programs (![#c5f015](https://github.com/Theo-Nelson/squares/blob/main/pink_square.png)) are optional in this process but additive in contextualizing your long-read sample, and blue programs (![#1589F0](https://github.com/Theo-Nelson/squares/blob/main/blue_square.png)) are recommended for users with external web-hosting resources. Additionally, programs are rated based on speed, with options being "Super Quick", "Quick" and "Grab a Coffee ![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/SMirC-coffeebreak.svg.png)". The approximate timing should be less than 1 min, less than 10 minutes and greater than 10 minutes, respectively. As you increase the number of reads you analyze, the time will increase proportionally for many of the programs listed. A key exception is the installation script, which should consistently take betweeen 12 and 15 minutes.  

#### Parameter Input and User Instructions ![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/red_square.png) (Super Quick)
This section allows you to setup and automate the remainder of the analysis by providing key pieces of information. Most parameters are necessary for specific programs, but four variables in particular are more widely applicable: ```PIPELINE_FILE_PATH```, ```ACC```, ```INDEX_FILE_PATH``` and ```ANNOTATION_FILE_PATH```. Please note that ```PARTITIONED_INDEX_FILE_PATH``` is necessary if you are using minimap2 featherweight alignment which is described below. The pipeline file path describes where in your Google Drive file system or local file system the pipeline is set up. The ACC variable specifies either the FASTQ file on your machine containing long-read sequencing data or the run accession number (e.g. [SRR12389274](https://www.ebi.ac.uk/ena/browser/view/SRR12389274)) for the publically available file which you wish to analyze. You can search for samples relevant to your field of interest with the European Nucleotide Archive's Advanced Search Feature: https://www.ebi.ac.uk/ena/browser/advanced-search (tutorial: https://www.youtube.com/watch?v=ugLaYRgh1pE). The index file path describes where the reference genome exists within your file system; the annotation file path describes where the reference annotation exists in the same file system. 

#### Mounting your Google Drive ![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/red_square.png) (Super Quick)
This will connect your current hardware instance to your Google Drive and allow you to permeantly store your analysis. Most academic users enjoy unlimited Google Drive Storage space, while basic users are given 30 GB of free storage space.

#### Managing Software via BioConda ![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/red_square.png) (Grab a Coffee ![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/SMirC-coffeebreak.svg.png))
BioConda is a centralized package manager necessary to install the remaining programs. 

#### Kingfisher: fast and flexible program for procurement of sequence files ![#c5f015](https://github.com/Theo-Nelson/squares/blob/main/pink_square.png) (Quick)
This program will download sequence files from the European Nucleotide Archive. 

#### FastQC: A quality control tool for high throughput sequence data ![#c5f015](https://github.com/Theo-Nelson/squares/blob/main/pink_square.png) (Quick)
This tool will generate basic high-level statistics regarding read length and sequence quality. Long-read sequences generally have poor base-pair level quality. 

#### minimap2: A versatile pairwise aligner for genomic and spliced nucleotide sequences ![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/red_square.png) (Grab a Coffee ![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/SMirC-coffeebreak.svg.png))
Alignment will occur either with minimap2 (available in Google Colab Pro due to memory constraints) or minimap2 featherweight alignment (available in the free version of Google Colab). minimap2 is a versatile aligner which maps reads onto a reference genome. This pipeline utilizes program options suitable for alignment of reads from cDNA libraries or RNA. Please consult the minimap2 documentation (https://github.com/lh3/minimap2) to find options relating to your specific application.  

#### samtools: Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format ![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/red_square.png) (Super Quick)
This tool stores your aligned reads in a compressed binary format, sorting and indexing them along the way for quick access by chromosomal location. These files can be downloaded and viewed in a program such as the Integrated Genomics Viewer available from the Broad Institute (https://software.broadinstitute.org/software/igv/). 

#### TranscriptClean: correct mismatches, microindels, and noncanonical splice junctions ![#c5f015](https://github.com/Theo-Nelson/squares/blob/main/pink_square.png) (Grab a Coffee ![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/SMirC-coffeebreak.svg.png))
This feature will try to polish your reads to adhere to known biological principles and decrease variability among reads aligning to the same region.  

#### featureCounts: an efficient general purpose program for assigning sequence reads to genomic features ![#c5f015](https://github.com/Theo-Nelson/squares/blob/main/pink_square.png) (Super Quick)
This tool will produce a count matrix assigning reads to features within your reference annotation. 

#### svist4get: a simple visualization tool for genomic tracks from sequencing experiments ![#c5f015](https://github.com/Theo-Nelson/squares/blob/main/pink_square.png) (Quick)
This tool will generate a graph of read coverage (i.e. how many reads align to a given region) for a specific chromosomal region (the names for these must match the names of chromosome within your reference genome/annotation). 

#### Pistis: Quality control plotting for long reads ![#c5f015](https://github.com/Theo-Nelson/squares/blob/main/pink_square.png) (Grab a Coffee ![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/SMirC-coffeebreak.svg.png))
This quality control package will provide information regarding read quality, gc content and read alignment. 

#### MakeHub: Fully automated generation of UCSC assembly hubs ![#1589F0](https://github.com/Theo-Nelson/squares/blob/main/blue_square.png) (Grab a Coffee ![#f03c15](https://github.com/Theo-Nelson/squares/blob/main/SMirC-coffeebreak.svg.png))
This feature will create a track hub which you can host on a public file service. This can then be connected to the UCSC genome browser and viewed as a track (e.g.: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr12%3A116533435%2D116536513&hgsid=1282820889_LUJAMUR9DzBvxiKtV6M3Qhk7c0Iv). Please see the MakeHub documentation for more details: https://github.com/Gaius-Augustus/MakeHub#how-to-use-makehub-output-with-ucsc-genome-browser.

#### MultiQC: Aggregate results from bioinformatics analyses across many samples into a single report ![#c5f015](https://github.com/Theo-Nelson/squares/blob/main/pink_square.png) (Quick)
This tool will summarize results from other quality control tools into a one-page report.

---

### General Usage 

After setting the parameters as needed, you should hit "connect" in order to be provided a hardware allocation within Google Colaboratory. Save the set parameters into the runtime by hitting the arrow on the left side of each code bar. Thereafter, minimize the remaining sections as shown below. You can then run the relevant analysis automatically by utilizing the arrow referencing the entire section. 

![screenshot of the runtime view](https://u.cubeupload.com/MakeTheBrainHappy/ScreenShot20220214at.png)

---

### Details on the output 

The pipeline will save all outputs to the relevant folders. 

---

### Troubleshooting Guide

to be added.

---

### Citation

to be added.

---

### Contributors

* [Theodore Nelson](https://github.com/Theo-Nelson), Columbia University Irving Medical Center 
