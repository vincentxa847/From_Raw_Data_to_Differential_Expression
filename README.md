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

Differential expression (DE) analysis is conducted using [DESeq2](https://doi.org/10.1186/s13059-014-0550-8), which is featured with the shrinkage for dispersion and fold-change estimation to deal 
with small sample sizes and heteroskedasticity. (*Differential_expression_analysis.R*)
### Dataset Used
Raw data files 
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
The reference genome in form of fasta file represents chromosome 2 of mm10 mouse genome:
