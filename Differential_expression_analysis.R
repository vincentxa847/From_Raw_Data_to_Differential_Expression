library(DESeq2)
library(vsn)
library(apeglm) # Approximate posterior estimation for GLM coefficients

#### Create dds object for gene and transcript and plot dispersion plot,PCA plot ####
# load gene-count matrix
countTable_gene = read.csv("2802815_gene_count_matrix_nondisc.csv",row.names=1)
colTable_gene = read.csv("Design-Exp.csv",row.names=1) # s10,s11 and s12 are placed behind s1 because the arrangement of gene_count_matrix

dim(countTable_gene)[1] # calculate how many gene : 2014

## Generate Deseq Data Set (dds) object
dds_gene = DESeqDataSetFromMatrix(countData=countTable_gene,colData=colTable_gene,design= ~ group)
notAllZero_gene = (rowSums(counts(dds_gene)) > 0) # keep only the sum of counts in a row bigger than zero
dds_gene = dds_gene[notAllZero_gene,]

# Run the differential expression (DE) analysis for gene
DE_gene = DESeq(dds_gene,test="Wald") # create dds object "DE_gene"

# make dispersion plot for gene
plotDispEsts(DE_gene, main="Gene Dispersion Estimates")

# make rlog transformation for gene
rld_gene = rlog(DE_gene,blind = TRUE) # blind means unbiased, for performing QC

# make rlog-based PCA plot for gene
plotPCA(rld_gene,intgroup=c("group")) # group is assigned in dds design

# using vst transformation to make PCA plot for gene
vst_gene = vst(DE_gene,blind = TRUE) # blind means unbiased, for performing QC
plotPCA(vst_gene,intgroup=c("group"))

## generate dds object for transcript counts
countTable_transcript = read.csv("2802815_transcript_count_matrix_nondisc.csv",row.names=1)
colTable_transcript = read.csv("Design-Exp.csv",row.names=1)

dim(countTable_transcript)[1] # calculate how many transcript 3063
        
dds_transcript = DESeqDataSetFromMatrix(countData=countTable_transcript,colData=colTable_transcript,design= ~ group)
notAllZero_transcript = (rowSums(counts(dds_transcript)) > 0) # keep only the sum of counts in a row bigger than zero
dds_transcript = dds_transcript[notAllZero_transcript,]

# Run the differential expression (DE) analysis for transcript
DE_transcript = DESeq(dds_transcript,test="Wald") # create dds object "DE_transcript"
        
# make dispersion plot for transcript
plotDispEsts(DE_transcript, main="Transcript Dispersion Estimates")
# make rlog transformation for transcript
rld_transcript = rlog(DE_transcript,blind = TRUE) # blind means unbiased, for performing QC
        
# make rlog-based PCA plot for transcript
plotPCA(rld_transcript,intgroup=c("group")) # group is assigned in dds design

#### Make SD verus mean plots for log and rlog ####
# Belong to the section of "Data transformations and visualization" 

# Perform transformation 
# generate log2(raw-counts), log2(normalized-counts) 
lgc.raw_gene = log2(counts(DE_gene,normalized=FALSE)+0.01) # +0.01 is to avoid log2 of zero
lgc.norm_gene = log2(counts(DE_gene,normalized=TRUE)+0.01)

# Visualization, need package "vsn"
meanSdPlot(lgc.raw_gene)  
meanSdPlot(lgc.norm_gene)
meanSdPlot(assay(rld_gene)) # need "assay" to extract rlog transformation data