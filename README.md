# From Raw Data to Differential Expression
## This project demonstrates fundamental RNA-Seq data processing skills by using the NewTuxedo cascade and DESeq2 to convert raw data into a gene count table and perform differential expression analysis.

## Method
Raw reads undergo Read pre-processing and new Tuxedo 2 cascade of Hisat2 and Stringtie to produce expression matrix that features by samples for differential expression analysis. (*From Raw Data to Gene Count Table.sh*)

[Scythe](https://github.com/vsbuffalo/scythe) is for removing 3’ end adaptor, which using poor quality bases to identify and remove adaptor. 
[Sickle](https://github.com/najoshi/sickle) is for quality trimming, which using window based algorithms to remove bases with low quality that accumulate in both ends of read. 
Scythe and Sickle compose the pre processing of the procedure to improve the quality of alignment subsequent to the pre-processing. 

[Hisat2](https://doi.org/10.1038/s41587-019-0201-4) is utilized to perform alignment, which is featured with the using of graph Ferragina Manzini (gFM) index for aligning reads to graph, 
Hisat2 also adopts hierarchical indexing from the HISAT to accelerate the alignment and reduce memory requirements. 
[Samtools](https://doi.org/10.1093/gigascience/giab008) is adopted to convert the sam file from hisat2 to bam file and sort the bam file for further analysis. 
[Stringtie](https://doi.org/10.1038/nprot.2016.095) that computing bases on the building of flow network is utilized to perform assembly and quantification of transcripts. Then, 
prepDE.py from Stringtie turns the GTF file producing by stringtie into expression matrix.

Differential expression (DE) analysis is performed using [DESeq2](https://doi.org/10.1186/s13059-014-0550-8), which is featured with the shrinkage for dispersion and fold-change estimation to deal 
with small sample sizes and heteroskedasticity. (*Differential_expression_analysis.R*)
### Dataset Used
Raw data files (not provided here)
- s1.c2.fq -> group “A” replicate 1 
- s2.c2.fq -> group “A” replicate 2 
- s3.c2.fq -> group “A” replicate 3 
- s4.c2.fq -> group “A” replicate 4 
- s5.c2.fq -> group “B” replicate 1 
- s6.c2.fq -> group “B” replicate 2 
- s7.c2.fq -> group “B” replicate 3 
- s8.c2.fq -> group “B” replicate 4 
- s9.c2.fq -> group “C” replicate 1 
- s10.c2.fq -> group “C” replicate 2 
- s11.c2.fq -> group “C” replicate 3 
- s12.c2.fq -> group “C” replicate 4

These data sets contain sequencing reads generated with single-end sequencing using NextSeq500 sequencer. The reads are preselected so that they originate from chromosome 2 of mouse genome. 

Reference genome : 
The reference genome in form of fasta file represents chromosome 2 of mm10 mouse genome

## Result 
From the result of error file of *From Raw Data to Gene Count Table.sh*, It can be seen from the error file that about 1500 adapters were identified and trimmed 
for each group, adapter contamination rates of twelve groups range from only 0.0002 to 0.0003 before trimming. It seems that reads had already trimmed before. 
Most of the contamination come from universal adapters, with the purpose of anchoring reads on flow cell. Also, no read was discarded due to the length of read below threshold after 
quality trimming, which can be deduced that the quality of raw read is good. The percentage of read aligned exactly one time of each group range from 84 to 89, demonstrating a good quality of alignment. 

Dispersion plots of gene and transcript show the shrinkage of gene-wise estimates toward the consensus represented by red line. It can be seen that the 
dispersion of transcript is generally variable than gene, which means noise is likely to dominate the biologically meaningful signal, therefore shrinkage is important for 
the transcript to ensure an accurate DE analysis. [However, the effect of shrinkage on  transcript is less extent when comparing with the effect of shrinkage on gene, this 
is because the ambiguity in mapping fragments to transcripts](https://doi.org/10.1038/nbt.2450) 
![image](https://github.com/vincentxa847/From_Raw_Data_to_Differential_Expression/assets/118545004/751ebdcd-a403-44ec-8d77-6e490e64704c)\
*Figure1: Dispersion plot of gene and transcript. Black dots are the maximum-likelihood estimate (MLE), Curve (red) is fit to the MLE to depict the overall trend of dispersion 
mean dependence. Blue dots are the maximum a posteriori that come from the second estimation using fit. For gene’s and transcript’s dispersions that show significant above 
the curve will be viewed as outliner (black dots circled in blue) and not undergo shrinkage.*

rlog-based PCA plot of gene shows a similarity between group A and B, group C show difference when comparing with group A and B. When comparing 
the rlog-based PCA plot of gene with transcript, it is observed that the variance of transcript within group is significant than gene, which is consistent with the result 
of Figure1 that dispersion of transcript is greater than gene. 
![image](https://github.com/vincentxa847/From_Raw_Data_to_Differential_Expression/assets/118545004/bff27090-c760-45ed-a146-d9804c9d1909)\
*Figure 2: rlog-based PCA plot of gene. Each group has three replicates. Similarity is observed between the group A and C.*

Regularized logarithm (rlog) transformation behaves similarly to standard logarithm transformation for genes with high counts but shrinking the genes with low counts. This 
can avoid noise dominate biologically meaningful signal in low read counts data, which is common when standard logarithm transformation is used to transform data. 
SD versus mean plots show variance of each gene across samples. It can be seen that the deviation of variance is stable across the range of mean after rlog 
transformation, which is the effect of regularized logarithm transformation on the variance that makes the data homoskedastic.  
![image](https://github.com/vincentxa847/From_Raw_Data_to_Differential_Expression/assets/118545004/ecea2690-ecdc-44f6-9708-6444e1ea24e2)
![image](https://github.com/vincentxa847/From_Raw_Data_to_Differential_Expression/assets/118545004/198766b2-8e27-41d5-9a18-b58c01db411f)
*Figure 3: SD versus mean plots of data after standard logarithm transformation (upper) and Regularized logarithm transformation (lower). SD on y-axis is the variance across all samples, therefore the SD also depicts the variance of experimental conditions that changes the shape of curve.* 

MA-plots of the DE analysis between group B and group A show that genes with low read counts have stronger difference between the compared groups than 
genes with high read counts. This phenomenon called “heteroskedasticity”, which is a common issue of count data that noise appears when counts are low. 
![image](https://github.com/vincentxa847/From_Raw_Data_to_Differential_Expression/assets/118545004/e8a440a1-b2ae-479a-8d7c-d8bc2e0f6639)
*Figure 4A:  MA-plots of the differential expression analysis between group B and group A*

Shrinking of LFC by DESeq2 can overcome the problem. Shrinking of LFC by DESeq2 is stronger to genes with less information that result from low read counts 
or high dispersion, this adjustment can correct the exaggeration of LFC of low read counts. 
![image](https://github.com/vincentxa847/From_Raw_Data_to_Differential_Expression/assets/118545004/ea66e275-da96-4416-9921-7dd292b95a12)
*Figure 4B:  MA-plots of the differential expression analysis between group B and group A, empirical Bayes shrinkage for 
fold-change estimation is adopted*

In order to deal with genes that have statistically significance but weak in effect strength, DESeq2 can set the LFC threshold to directly filter genes with log fold change above 
specified level (Figure 4). DE expression analysis of group B and A against the null hypothesis of LFC = 0 has 275 genes with significant upregulation and 259 genes with 
significant downregulation. In contrast, null hypothesis of LFC = 1 has only 45 genes with significant upregulation and 2 genes with significant downregulation. Specifying 
minimum effect size can not only filter out genes that poor in effect strength but also narrow down the number of target genes for downstream analysis. 
