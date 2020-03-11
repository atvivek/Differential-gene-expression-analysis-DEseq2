# Differential-gene-expression-analysis-DEseq2
Welcome to the RNA-seq Tutorial. Use this page to navigate your way through all exercises. Each page has a link at the bottom to bring you back to this table of contents.

Arabidopsis thaliana is a small flowering plant that is used as a model system in research of plant biology. It helped researchers to build basic understanding around molecular, biochemical and genetics processes in the plants. A wealth of knowledge and information is available around Arabidopsis genomics (genome sequence, transcriptome, genetic markers etc) and hence could be used as ideal system to develop fundamental understanding plant RNA-seq analysis before venturing in the transcriptomics world of non model species. 

In this tutorial we will be using RNAseq dataset of A. thaliana and the study was published in "APS publications" (http://apsjournals.apsnet.org/doi/full/10.1094/MPMI-07-15-0156-R). 

	Shanks CM et al., "The Role of Cytokinin During Infection of Arabidopsis thaliana by the Cyst Nematode Heterodera schachtii.", Mol Plant Microbe Interact, 2016 Jan;29(1):57-68

Total 6 RNAseq datasets representing 3 biological replicates each for Control and infected samples were used in this study The RNA profiles are archived in the SRA, and meta Information on each may be viewed through the SRP ID SRP063017 (https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP063017).

The Sequence Read Archive, or SRA, is a publicly available database containing read sequences from a variety of experiments. Scientists who would like their read sequences present on the SRA submit a report containing the read sequences, experimental details, and any other accessory meta-data.

Our data, SRR2221833,SRR2221837 , SRR2221841, SRR2221834, SRR2221838, SRR2221842 come from CON1, CON2, CON33, INF1, INF2, and INF3 respectively. Our objective is to identify genes which are differentially expressed between WT and EE sample conditions. We have 3 replicates for each condition CON1, CON2, CON3 and INF1, INF2, INF3 for wildtype control and 4 dpi of Beet Cyst Nematode (Heterodera schachtii) respectively. All the data sets are single-end sequences.

### Download the data using sra-toolkit
We know that the SRA contain the read sequences and accessory meta Information from experiments. Rather than downloading experimental data through a browser, we may use the sratoolkit's "fastq-dump" function to directly dump raw read data into the current terminal directory. Let's have a look at this function (it is expected that you have read the Xanadu tutorial, and are familiar with loading modules):

To check the options you can simply type fastq-dump once you load the module in the terminal window.

    fastq-dump
Which will show you the following options it has:

      fastq-dump [options] <path> [<path>...]
      fastq-dump [options] <accession>

        Use option --help for more  Information

        fastq-dump : 2.8.2 
        
    fastq-dump --split-files SRR2221833 SRR2221837 SRR2221841 SRR2221834 SRR2221838 SRR2221842
    mv SRR2221833.fastq con_Rep1.fastq
    mv SRR2221834.fastq inf_Rep1.fastq
    mv SRR2221837.fastq con_Rep2.fastq
    mv SRR2221838.fastq inf_Rep2.fastq
    mv SRR2221841.fastq con_Rep3.fastq
    mv SRR2221842.fastq inf_Rep3.fastq
    
    mkdir con_Rep1 con_Rep2 con_Rep3 inf_Rep1 inf_Rep2 inf_Rep3 
    
    mv con_Rep1.fastq con_Rep1
    mv inf_Rep1.fastq inf_Rep1
    mv con_Rep2.fastq con_Rep2
    mv inf_Rep2.fastq inf_Rep2
    mv con_Rep3.fastq con_Rep3
    mv inf_Rep3.fastq inf_Rep3
    
### Quality control using fastp

    fastqc con_Rep1/con_Rep1.fastq con_Rep2/con_Rep2.fastq  con_Rep3/con_Rep3.fastq  inf_Rep1/inf_Rep1.fastq inf_Rep2/inf_Rep2.fastq inf_Rep3/inf_Rep3.fastq

    fastp -i con_Rep1/con_Rep1.fastq  -o con_Rep1/con_Rep1_trimmed.fastq 
    fastp -i inf_Rep1/inf_Rep1.fastq  -o inf_Rep1/inf_Rep1_trimmed.fastq 
    fastp -i con_Rep2/con_Rep2.fastq  -o con_Rep1/con_Rep2_trimmed.fastq 
    fastp -i inf_Rep2/inf_Rep2.fastq  -o inf_Rep2/inf_Rep2_trimmed.fastq 
    fastp -i con_Rep3/con_Rep3.fastq  -o con_Rep3/con_Rep3_trimmed.fastq 
    fastp -i inf_Rep3/inf_Rep3.fastq  -o inf_Rep3/inf_Rep3_trimmed.fastq 
     
  
    Detecting adapter sequence for read1...
    No adapter detected for read1

    Read1 before filtering:
    total reads: 11876535
    total bases: 593826750
    Q20 bases: 588695215(99.1359%)
    Q30 bases: 578204450(97.3692%)

    Read1 after filtering:
    total reads: 11863788
    total bases: 593189400
    Q20 bases: 588382234(99.1896%)
    Q30 bases: 577931353(97.4278%)

    Filtering result:
    reads passed filter: 11863788
    reads failed due to low quality: 11139
    reads failed due to too many N: 1608
    reads failed due to too short: 0
    reads with adapter trimmed: 0
    bases trimmed due to adapters: 0

    Duplication rate (may be overestimated since this is SE data): 49.058%

    JSON report: fastp.json
    HTML report: fastp.html

    fastp -i con_Rep1.fastq -o con_Rep1_trimmed.fastq 
    fastp v0.20.0, time used: 53 seconds
    
    

    fastqc con_Rep1/con_Rep1_trimmed.fastq con_Rep2/con_Rep2_trimmed.fastq  con_Rep3/con_Rep3_trimmed.fastq  inf_Rep1/inf_Rep1_trimmed.fastq inf_Rep2/inf_Rep2_trimmed.fastq inf_Rep3/inf_Rep3_trimmed.fastq
    
    
https://plants.ensembl.org/Arabidopsis_thaliana/Info/Index
    
    wget  ftp://ftp.ensemblgenomes.org/pub/plants/release-46/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
    gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz 
    mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ath10.fa
    
    wget ftp://ftp.ensemblgenomes.org/pub/plants/release-46/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gff3.gz
    gunzip Arabidopsis_thaliana.TAIR10.46.gff3.gz 
    mv Arabidopsis_thaliana.TAIR10.46.gff3 ath10.gff3
    
    
    hisat2
    hisat2-build ath10.fa ath
    
    hisat2 -p 8 --dta -x ath10 -U con_Rep1/con_Rep1_trimmed.fastq -S con_Rep1/con_Rep1.sam
    hisat2 -p 8 --dta -x ath10 -U inf_Rep1/inf_Rep1_trimmed.fastq -S inf_Rep1/inf_Rep1.sam
    hisat2 -p 8 --dta -x ath10 -U con_Rep2/con_Rep2_trimmed.fastq -S con_Rep2/con_Rep2.sam
    hisat2 -p 8 --dta -x ath10 -U inf_Rep2/inf_Rep2_trimmed.fastq -S inf_Rep2/inf_Rep2.sam
    hisat2 -p 8 --dta -x ath10 -U con_Rep3/con_Rep3_trimmed.fastq -S con_Rep3/con_Rep3.sam
    hisat2 -p 8 --dta -x ath10 -U inf_Rep3/inf_Rep3_trimmed.fastq -S inf_Rep3/inf_Rep3.sam

command 
-p : number of processors been used
--dta: report alignments tailored for transcript assemblers
-x: path to index generated from previous step
-q: query input files in fastq format
-S: output SAM file

    Usage:   samtools  [options]

    Commands:
      -- Indexing
         dict           create a sequence dictionary file
         faidx          index/extract FASTA
         index          index alignment

      -- Editing
         calmd          recalculate MD/NM tags and '=' bases
         fixmate        fix mate  Information
         reheader       replace BAM header
         targetcut      cut fosmid regions (for fosmid pool only)
         addreplacerg   adds or replaces RG tags
         markdup        mark duplicates

      -- File operations
         collate        shuffle and group alignments by name
         cat            concatenate BAMs
         merge          merge sorted alignments
         mpileup        multi-way pileup
         sort           sort alignment file
         split          splits a file by read group
         quickcheck     quickly check if SAM/BAM/CRAM file appears intact
         fastq          converts a BAM to a FASTQ
         fasta          converts a BAM to a FASTA

      -- Statistics
         bedcov         read depth per BED region
         depth          compute the depth
         flagstat       simple stats
         idxstats       BAM index stats
         phase          phase heterozygotes
         stats          generate stats (former bamcheck)

      -- Viewing
         flags          explain BAM flags
         tview          text alignment viewer
         view           SAM<->BAM<->CRAM conversion
         depad          convert padded BAM to unpadded BAM
         
         
     Usage: samtools sort [options...] [in.bam]
 
    Options:
      -l INT     Set compression level, from 0 (uncompressed) to 9 (best)
      -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]
      -n         Sort by read name
      -t TAG     Sort by value of TAG. Uses position as secondary index (or read name if -n is set)
      -o FILE    Write final output to FILE rather than standard output
      -T PREFIX  Write temporary files to PREFIX.nnnn.bam
          --input-fmt-option OPT[=VAL]
                   Specify a single input file format option in the form
                   of OPTION or OPTION=VALUE
      -O, --output-fmt FORMAT[,OPT[=VAL]]...
                   Specify output format (SAM, BAM, CRAM)
          --output-fmt-option OPT[=VAL]
                   Specify a single output file format option in the form
                   of OPTION or OPTION=VALUE
          --reference FILE
                   Reference sequence FASTA FILE [null]
      -@, --threads INT
                   Number of additional threads to use [0]
                   
      
      
      samtools view -@ 8 -bhS con_Rep1/con_Rep1.sam -o con_Rep1/con_Rep1.bam
      samtools sort -@ 8 con_Rep1/con_Rep1.bam -o con_Rep1/con_Rep1_sort.bam
      
      samtools view -@ 8 -bhS inf_Rep1/inf_Rep1.sam -o inf_Rep1/inf_Rep1.bam
      samtools sort -@ 8 inf_Rep1/inf_Rep1.bam -o inf_Rep1/inf_Rep1_sort.bam
      
      samtools view -@ 8 -bhS con_Rep2/con_Rep2.sam -o con_Rep2/con_Rep2.bam
      samtools sort -@ 8 con_Rep2/con_Rep2.bam -o con_Rep2/con_Rep2_sort.bam
      
      samtools view -@ 8 -bhS inf_Rep2/inf_Rep2.sam -o inf_Rep2/inf_Rep2.bam
      samtools sort -@ 8 inf_Rep2/inf_Rep2.bam -o inf_Rep2/inf_Rep2_sort.bam
      
      samtools view -@ 8 -bhS con_Rep3/con_Rep3.sam -o con_Rep3/con_Rep3.bam
      samtools sort -@ 8 con_Rep3/con_Rep3.bam -o con_Rep3/con_Rep3_sort.bam
      
      samtools view -@ 8 -bhS inf_Rep3/inf_Rep3.sam -o inf_Rep3/inf_Rep3.bam
      samtools sort -@ 8 inf_Rep3/inf_Rep3.bam -o inf_Rep3/inf_Rep3_sort.bam
     
     
   ### Reference Guided Transcript Assembly
   
   Stringtie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. It can be executed in 3 different modes
   
   
  1.Exclusively reference guided : In this mode stringtie quantify the expression of known transcripts only.
  
  2.Reference guided transcript discovery mode : Quantify known transcripts and detect novel ones.
  
  3.De-novo mode : Detect and assemble transcripts.

We will be running stringtie using the option 2 .In the first step of this process stringtie along with sample bam file and reference gtf file generate a gtf file corresponding to the sample. This gtf file has information on expression levels of transcripts, exons and other features along with any novel transcripts. The syntax of the command is


    stringtie -p 4 -l label -G Reference.gtf -o sample.gtf sample.bam
    
In this command

    -p specifies the number of threads to use.
    -l label used in gtf file
    -G Reference GTF available from public domain databases
    -o output gtf corresponding to expression levels of features of the sample
    
 Once we have run this command through all our six samples (con1, con2, con3, inf1,inf2 and inf3) we will have 6 gtf files each corresponding to one of the sample containing feature expression values. Having 6 different gtf files is of no advantage as each may contain the same novel transcript but labelled differently. Ideally we would like to merge these 6 gtf files along with the reference GTF to achieve following goals-
 
 
-- Redundant transcripts across the samples should be represented once
-- Known transcripts should hold their stable gene ID's (assigned in Ensembl)

The command we will use to achieve this is stringtie --merge and the syntax is




    
  


### *Installation Requirements*

1.Download the most recent versions of R and RStudio for your laptop:
fastqc
    R
    RStudio

2.Install the following packages using the instructions provided below.

    NOTE:  When installing the following packages, if you are asked to select (a/s/n) or (y/n), please select “a” or "y" as applicable but know that it can take awhile.

(a) Install the below packages on your laptop from CRAN. You DO NOT have to go to the CRAN webpage; you can use the following function to install them one by one:

install.packages("insert_first_package_name_in_quotations")
install.packages("insert__second_package_name_in_quotations")
& so on ...

Packages to install from CRAN (note that these package names are case sensitive!):

    BiocManager
    RColorBrewer
    pheatmap
    ggrepel
    devtools
    tidyverse

(b) Install the below packages from Bioconductor, using BiocManager::install() function 7 times for the 7 packages:

    BiocManager::install("insert_first_package_name_in_quotations")
    BiocManager::install("insert_second_package_name_in_quotations") 

Packages to install from Bioconductor (note that these package names are case sensitive!):

    DESeq2
    clusterProfiler
    DOSE
    org.Hs.eg.db
    pathview
    DEGreport
    EnsDb.Hsapiens.v86
    AnnotationHub
    ensembldb

   Finally, please check that all the packages were installed successfully by loading them one at a time using the library() function.

    library(DESeq2)
    library(ggplot2)
    library(RColorBrewer)
    library(pheatmap)
    library(ggrepel)
    library(clusterProfiler)
    library(DEGreport)
    library(org.Hs.eg.db)
    library(DOSE)
    library(pathview)
    library(tidyverse)
    library(EnsDb.Hsapiens.v86)
    library(AnnotationHub)
    library(ensembldb)

Once all packages have been loaded, run sessionInfo().

    sessionInfo()
