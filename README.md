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
    
    wget ftp://ftp.ensemblgenomes.org/pub/release-46/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.46.gtf.gz
    gunzip Arabidopsis_thaliana.TAIR10.46.gtf.gz 
    mv Arabidopsis_thaliana.TAIR10.46.gtf ath10.gtf
    
    
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

	stringtie --merge -p 4 -o stringtie_merged.gtf -G Reference.gtf listOfSampleGTFs.txt

options used:

	-p specifies the number of threads to use
	-G Reference GTF available from public domain databases
	-o output merged gtf file

listOfSampleGTFs.txt : This is a text file with list of gtfs generated freom the samples in previous step.

    ls -1 ath*/*.gtf >> sample_assembly_gtf_list.txt
    
The command above is to generate listOfSampleGTFs.txt that will be used in stringtie --merge command. The merged GTF can be compared with Reference GTF to get some stats on the stringtie_merged.gtf. The above set of commands can be put together as shown below,

	stringtie -p 8 -l con_Rep1 -G ath10.gtf -o con_Rep1/transcripts.gtf con_Rep1_sort.bam
        stringtie -p 8 -l inf_Rep1 -G ath10.gtf -o inf_Rep1/transcripts.gtf inf_Rep1_sort.bam
        stringtie -p 8 -l con_Rep2 -G ath10.gtf -o con_Rep2/transcripts.gtf con_Rep2_sort.bam
	stringtie -p 8 -l inf_Rep2 -G ath10.gtf -o inf_Rep2/transcripts.gtf inf_Rep2_sort.bam
	stringtie -p 8 -l con_Rep3 -G ath10.gtf -o con_Rep3/transcripts.gtf con_Rep3_sort.bam
	stringtie -p 8 -l inf_Rep3 -G ath10.gtf -o inf_Rep3/transcripts.gtf inf_Rep2_sort.bam
	
	ls -1 *.gtf | grep transcripts >> sample_assembly_gtf_list.txt
	
	stringtie --merge -p 8 -o stringtie_merged.gtf -G ath10.gtf sample_assembly_gtf_list.txt
	
	gffcompare -r ath10.gtf -o gffcompare stringtie_merged.gtf

Now lets examine the outputs generated from above comands. As discussed above in first step stringtie genrates a gtf file for each sample with details of coverage, FPKM, TPM and other information on the transcripts based on sample bam file.
	

Now lets have a look at out merged GTF file stringtie_merged.gtf from the previous step:


This is our new reference GTF file we will be using to quantify the expression of dfferent genes and transcripts. If we look closer we can see that the file have information different features but exclude coverage, TPM and FPKM information. Thats how we want it to be for use as reference in subsequent analysis. Also note that the first two transcripts have known ENSEMBL transcrip-id,gene_name and ref_gene_id, however it is missing in transcript 3. This is because it represents a novel transcript identified in the study.

Now lets go ahead and do the transcript quantification using stringtie.

### Transcript quantification with StringTie

In this step we will use the stringtie_merged.gtf file as reference and measure the expression of exons, transcripts and other features present in the gtf file. The syntax of command we will be executing is,

	stringtie -e -B -p 4 sample.bam -G stringtie_merged.gtf -o output.count -A gene_abundance.out
Where:

	-B returns a Ballgown input table file
	-e only estimate the abundance of given reference transcripts
	-o output path/file name
	-A gene abundance estimation output file

Type in the below based on the above command to run all our samples.

	stringtie -e -B -p 4 con_Rep1/con_Rep1_sort.bam -G stringtie_merged.gtf -o con_Rep1/con_Rep1.count -A con_Rep1/con_Rep1_gene_abun.out
	stringtie -e -B -p 4 inf_Rep1/con_Rep1_sort.bam -G stringtie_merged.gtf -o inf_Rep1/con_Rep1.count -A inf_Rep1/inf_Rep1_gene_abun.out
	stringtie -e -B -p 4 con_Rep2/con_Rep2_sort.bam -G stringtie_merged.gtf -o con_Rep2/con_Rep2.count -A con_Rep2/con_Rep2_gene_abun.out
	stringtie -e -B -p 4 inf_Rep2/con_Rep2_sort.bam -G stringtie_merged.gtf -o inf_Rep2/con_Rep2.count -A inf_Rep2/inf_Rep2_gene_abun.out
      stringtie -e -B -p 4 con_Rep3/con_Rep3_sort.bam -G stringtie_merged.gtf -o con_Rep3/con_Rep3.count -A con_Rep3/con_Rep3_gene_abun.out
	stringtie -e -B -p 4 inf_Rep3/con_Rep3_sort.bam -G stringtie_merged.gtf -o inf_Rep3/con_Rep3.count -A inf_Rep3/inf_Rep3_gene_abun.out

DESeq2 and edgeR are two popular Bioconductor packages for analyzing differential expression, which take as input a matrix of read counts mapped to particular genomic features (e.g., genes). We provide a Python script (prepDE.py) to extract this read count information directly from the files generated by StringTie (run with the -e parameter).

Enter the following

	vim sample_lst.txt

Type the following 

	con1 <con_Rep1/transcripts.gtf.gtf>
	con2 <con_Rep2/transcripts.gtf.gtf>
	con3 <con_Rep3/transcripts.gtf.gtf>
	inf1 <inf_Rep1/transcripts.gtf.gtf>
	inf2 <inf_Rep2/transcripts.gtf.gtf>
	inf3 <inf_Rep3/transcripts.gtf>
	
	:wq

With this the information typed is saved.

Now,

	python prepDE.py -i sample_lst.txt
	
The script will produce two result files: gene_count_matrix.csv and transcript_count_matrix.csv.These are the input files for the differential expression analysis or the principal component analysis (PCA).

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


This step requires the R language and an IDE such as RStudio installed on a local machine. The R DESeq2 library also must be installed. To install this package, start the R console and enter:

	source("http://bioconductor.org/biocLite.R")
	biocLite("DESeq2")

If any dependencies fail, install them using the command: install.packages(PackageName, repos='http://cran.rstudio.com/')
Before running this script make sure to set the working directory and path to your Listeria_deseqFile (after copying it from the server to your local machine). The script was adapted slightly from Dave Wheeler's comprehensive tutorial on analysis with DESeq2. The only changes were a few bug fixes, adding an outputPrefix variable to allow easy modification of the output file names in the code for future use, and adding filtering by adjusted p value.
The most important information comes out as -replaceoutliers-results.csv. This file only contains the genes that have adjusted p values less than 0.05. These genes are the differentially expressed genes we are interested in. Depending on the experiment, this file can be adjusted to include p values less than 0.10 or a different value.

	# Import data from featureCounts
	countdata <- read.table("rawcounts.txt", header=TRUE, row.names=1)

	# Convert to matrix
	countdata <- as.matrix(countdata)
	head(countdata)

	# Assign condition (first two are WT, second two contain ABA, third two contain BAP, 
	#fourth two contain GAA, fifth two contain  IAA, sixth two contain S24)
	(condition <- factor(c(rep("WT", 2), rep("ABA", 2),rep("GAA", 2),rep("IAA", 2),rep("S24", 2),rep("BAP",2))))

	# Analysis with DESeq2 ----------------------------------------------------

	library(DESeq2)

	# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
	(coldata <- data.frame(row.names=rownames(countdata[,1]), condition))
	dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
	dds

	# Run the DESeq pipeline
	dds <- DESeq(dds)

	# Plot dispersions
	png("qc-dispersions.png", 1000, 1000, pointsize=20)
	plotDispEsts(dds, main="Dispersion plot")
	dev.off()

	# Regularized log transformation for clustering/heatmaps, etc
	rld <- rlogTransformation(dds)
	head(assay(rld))
	hist(assay(rld))


	# Colors for plots below
	## Ugly:
	## (mycols <- 1:length(unique(condition)))
	## Use RColorBrewer, better
	library(RColorBrewer)
	(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

	# Sample distance heatmap
	sampleDists <- as.matrix(dist(t(assay(rld))))
	library(gplots)
	png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
	heatmap.2(as.matrix(sampleDists), key=F, trace="none",
		  col=colorpanel(100, "black", "white"),
		  ColSideColors=mycols[condition], RowSideColors=mycols[condition],
		  margin=c(10, 10), main="Sample Distance Matrix")
	dev.off()

	# Principal components analysis
	## Could do with built-in DESeq2 function:
	## DESeq2::plotPCA(rld, intgroup="condition")
	## I like mine better:
	rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
	  require(genefilter)
	  require(calibrate)
	  require(RColorBrewer)
	  rv = rowVars(assay(rld))
	  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
	  pca = prcomp(t(assay(rld)[select, ]))
	  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
	  if (is.null(colors)) {
	    if (nlevels(fac) >= 3) {
	      colors = brewer.pal(nlevels(fac), "Paired")
	    }   else {
	      colors = c("black", "red")
	    }
	  }
	  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
	  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
	  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
	  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
	  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
	  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
	  legend(legendpos, legend=levels(fac), col=colors, pch=20)
	  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
	  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
	  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
	}
	png("qc-pca.png", 1000, 1000, pointsize=20)
	rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
	dev.off()

	# Get differential expression results
	res <- results(dds)
	table(res$padj<0.1)
	## Order by adjusted p-value
	res <- res[order(res$padj), ]
	## Merge with normalized count data
	resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
	names(resdata)[1] <- "Gene"
	head(resdata)
	## Write results
	write.csv(resdata, file="diffexpr-results.csv")

	## Examine plot of p-values
	hist(res$pvalue, breaks=50, col="red")

	## Examine independent filtering
	attr(res, "filterThreshold")
	#plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

	## MA plot
	## Could do with built-in DESeq2 function:
	## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
	## I like mine better:
	maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
	  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
	  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
	  if (labelsig) {
	    require(calibrate)
	    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
	  }
	}
	png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
	maplot(resdata, main="MA Plot")
	dev.off()

	## Volcano plot with "significant" genes labeled
	volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
	  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
	  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
	  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
	  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
	  if (labelsig) {
	    require(calibrate)
	    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
	  }
	  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
	}
	png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
	volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-5, 5))
	dev.off()


	#collapse technical replicates
	ddscol <- collapseReplicates(dds,dds$condition)

	ddscol

	res <-results(ddscol)
	res
	table(res$padj<0.05)
	rescol <- res[order(res$padj), ]
	## Merge with normalized count data
	rescoldata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
	names(rescoldata)[1] <- "Gene"
	head(rescoldata)
	## Write results
	write.csv(rescoldata, file="collapsed_diffexpr-results.csv")


