
# Get to know the system
``` {bash}
1. commands: sinfo -lN
Number nodes: 4
CPUS: 12 each
RAM: 121000

2. command: df -h /vol/COMPEPIWS/.
Volume size: 3.9T

3. command: echo $TMPDIR
/vol/COMPEPIWS/tmp

4. /vol/COMPEPIWS/groups/atacseq_1

```


# Tips and Tricks
Performed all tasks

# Working with tables on the command line, or awk
``` {bash}
1. Working on worker 1 node

2. command: head /vol/COMPEPIWS/pipelines/references/genome_genes.gtf

3. wc -l /vol/COMPEPIWS/pipelines/references/genome_genes.gtf
Number of lines: 671462

4. two possible: cat  /vol/COMPEPIWS/pipelines/references/genome_genes.gtf |awk '{print $3}' | grep "exon" | wc -l
cut -f 3 genome_genes.gtf | grep "exon"|wc -l
               
Result: 324748

5.cat  /vol/COMPEPIWS/pipelines/references/genome_genes.gtf |awk '{  if($5-$4>1000) print $3, $5-$4}'|grep "exon" | wc -l
awk '{if ($5-$4>1000) print $3,$5-$4}' genome_genes.gtf|grep "exon"|wc -l
Result: 19593

6. cat  /vol/COMPEPIWS/pipelines/references/genome_genes.gtf |awk '{ print $3, $10}' | grep "exon" |grep "Sox17"|wc -l
Result= 20

7. cat /vol/COMPEPIWS/pipelines/references/genome_genes.gtf | awk '{if($1=="chr2" && $3 =="exon") print $1}'|wc -l
Result= 29373

8. cut -f 3 genome_genes.gtf |sort|uniq -c
Result= 286783 CDS
        324748 exon
        29978 start_codon
        29953 stop_codon
        
9. cat /vol/COMPEPIWS/pipelines/references/genome_genes.gtf |sort $1 |sort $3 > double_sorted.txt

```
# Conda and Bioconda

a) base loaded
b) conda create -p atacseq1_ap
c) conda activate /vol/COMPEPIWS/groups/atacseq1/conda/atacseq1
d) conda install -c conda-forge r-ggplot2
   conda install -c bioconda fastqc
   conda install -c bioconda bedtools=2.22
e) bedtools --help
   fastqc --help
   
f) conda update -c bioconda bedtools

g) conda deactivate
  

# Basics in R 
``` {r}
#Load required libraries
library(ggplot2)
library(reshape2)
library(GenomicRanges)

3. file_names=list.files(path = "/vol/COMPEPIWS/groups/shared",pattern = ".txt",full.names = TRUE)
 
  dfs=lapply(file_names,FUN = function(i){
  read.table(file = i,header=TRUE)
  })
  
4. dim(dfs[[1]]) results=  1362728        6
   dim(dfs[[2]]) results=  1362728        6

5. colnames(dfs[[1]]) results= "Chr"      "Start"    "End"      "H3K27me3" "H3K36me3" "H3K9me3"
   colnames(dfs[[2]]) results= "Chr"      "Start"    "End"      "H3K27me3" "H3K36me3" "H3K9me3"

6. ## re-do this 
 data1=dfs[[1]]
   data2=dfs[[2]]
   df1_bases=data1[,3]-data1[,2]
   sum(df1_bases)
[1] 2725456000

   df2_bases=data2[,3]-data2[,2]
   sum(df2_bases)
[1] 2725456000

   
7. data1=cbind(data1,rep("kidney",nrow(data1)))
   colnames(data1)[7]="cell_type"
   data2=cbind(data2,rep("liver",nrow(data2)))
   colnames(data2)[7]="cell_type"
   combined_df= rbind(data1,data2)
   
8. dim(combined_df) result= 2725456       7

9. melted_data=melt(combined_df[1:nrow(combined_df),4:7],id="cell_type")
   dim(melted_data)   result= 8176368       3

10.  p= ggplot(melted_data,aes(x=value,color=cell_type))+geom_density()+facet_wrap(~cell_type)
     p
     ggsave("cell_type.pdf", plot = p)
     
```
![](https://i.imgur.com/ZlZZ9En.png)

```{r}

11. p= ggplot(melted_data,aes(x=value,color=cell_type))+geom_density()+facet_wrap(~cell_type) + xlim(0,100)
    p
    ggsave("cell_type_edit.pdf", plot = p)
```
![](https://i.imgur.com/jKN6Lh8.png)

```{r}

12. kidney_epi=data1[1:nrow(data1),4:6]
    liver_epi=data2[1:nrow(data2),4:6]
    lk_df=cbind(kidney_epi,liver_epi)

13. new_cols=colnames(lk_df)

 new_cols=sapply(1:6, function(i){
  if(i<=3){
  paste(new_cols[i],"_kidney",sep = "")
  }
  else{
    paste(new_cols[i],"_liver",sep = "")
  }
  
})

colnames(lk_df)=new_cols


14. pca=prcomp(lk_df)


15. pca$sdev[1]^2   result for PC1=869.9908
    pca$sdev[2]^2   result for PC2=486.0194
    
16. var_explained = pca$sdev^2 / sum(pca$sdev^2)

#create scree plot

qplot(c(1:6), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)
  
```
![](https://i.imgur.com/DFBqaGj.png)

```{r}
 
 #Cummulative scree plot
 cum_exp=cumsum(var_explained)
 
 qplot(c(1:6), cum_exp) +
  geom_line() +
  xlab("Principal Component") +
  ylab("Variance Explained") +
  ggtitle("Cumulative Scree Plot") +
  ylim(0, 1)
  
```
![](https://i.imgur.com/wXc7oCo.png)

```{r}

17. num_loads=dim(pca$rotation)[1] *dim(pca$rotation)[2]
    num_loads   result= 36
    

18. rots=as.data.frame(pca$rotation) 
    p= ggplot(rots,aes(x=PC1,y=PC2,color=rownames(rots)))+geom_point()
    p= p + labs(color = "Epi Marks")
    p  
```
![](https://i.imgur.com/EHy1Fva.png)

The variables are distributed into clusteres as we can see on the top right the both H3K36me3 marks from both samples
are clustered together. Same for the H3K9me3 on the bottom left. However the H3K27me3 are not clustered together.
    
```{r}

##GenomicRanges

1. ##For kidney
   gr <- GRanges(seqnames = data1[,1], ranges = IRanges(start = data1[,2], end= data1[,3]))
   values(gr) <- data1[,4:6]
   
   ##For liver
   gr2 <- GRanges(seqnames = data2[,1], ranges = IRanges(start = data2[,2], end= data2[,3]))
   values(gr2) <- data2[,4:6]    
   
2. gr=as.data.frame(gr)
   gr2=as.data.frame(gr2)
   gr_bases=gr[,3]-gr[,2]
  
  sum(gr_bases)
[1] 2725456000

  gr2_bases=gr2[,3]-gr2[,2]
  
 sum(gr2_bases)
[1] 2725456000

3. subset_gr=gr[which(gr[,1]=="chr2"),]

4. gr2 <- GRanges(seqnames = data2[,1], ranges = IRanges(start = data2[,2], end= data2[,3]))
   values(gr2) <- data2[,4:6]
   shift(gr2, 100)
   
5. gr <- GRanges(seqnames = data1[,1], ranges = IRanges(start = data1[,2], end= data1[,3]))
   values(gr) <- data1[,4:6]
   
   gr2 <- GRanges(seqnames = data2[,1], ranges = IRanges(start = data2[,2], end= data2[,3]))
   values(gr2) <- data2[,4:6]
   
   findOverlaps(gr, gr2) 
   
6. grl <- GRangesList("txA" = gr, "txB" = gr2)
   grl  
```

# 1 Nextflow ATAC-seq pipeline

## 1.1 Running the pipeline

#### 1. Organize your raw data
a) Raw fastq files
kidney_14.5_ATAC_1_R1.fastq.gz  
kidney_15.5_ATAC_2_R1.fastq.gz  
liver_15.5_ATAC_1_R1.fastq.gz
kidney_14.5_ATAC_1_R2.fastq.gz  
kidney_15.5_ATAC_2_R2.fastq.gz  
liver_15.5_ATAC_1_R2.fastq.gz
kidney_14.5_ATAC_2_R1.fastq.gz  
liver_14.5_ATAC_1_R1.fastq.gz   
liver_15.5_ATAC_2_R1.fastq.gz
kidney_14.5_ATAC_2_R2.fastq.gz  
liver_14.5_ATAC_1_R2.fastq.gz   
liver_15.5_ATAC_2_R2.fastq.gz
kidney_15.5_ATAC_1_R1.fastq.gz  
liver_14.5_ATAC_2_R1.fastq.gz
kidney_15.5_ATAC_1_R2.fastq.gz  
liver_14.5_ATAC_2_R2.fastq.gz

b) 2 files per cell_type per time point per replicate are there.
*_R{1,2} indicates the replicate of the sequence, that replicate 1 and 2 are two copies of each other.
#### 2. Create a samplesheet
a) nextflow_run directory created
b) nano samplesheet.csv
Manually entered the data row wise for each column(group,
replicate, fastq_1, fastq_2)
#### 3. Run the ATAC-Seq pipeline
a) i) Paired end
  ii) skip_diff_analysis = 'True'
        skip_preseq = 'True'
        skip_plot_profile = 'True'
        skip_plot_fingerprint = 'True'
        skip_spp = 'True'
iii) Blacklist file is used to remove artefacts such as PCR duplicates. Significant background noise is removed on removing these regoins. Removal of blacklist regions act an important quality measure.
These regions are found at centromeres. telomeres and satellite repeats
## 1.2 Quality Control (QC)

#### 1. fastqc report
a) 51

b) The “Per base sequence quality” tells us the range of quality values for all bases at each position in the FastQ file are    displayed in this view. Moreover, all samples have passed this check
   
c) Per Base Sequence Content plots out the proportion of each base position in a file for which each of the four normal DNA      bases has been called. In all of our samples this check has been failed. To illustrate, at the start of each read position in all the samples the difference between A and T, or G and C is greater than 20%.

#### 2. Reads alignment

a) Mapped reads
kidney_14.5_R1.mLb.clN.sorted.bam.flagstat: 3122196 
   
kidney_14.5_R1.mLb.mkD.sorted.bam.flagstat: 3785633
kidney_14.5_R2.mLb.clN.sorted.bam.flagstat: 3608726
kidney_14.5_R2.mLb.mkD.sorted.bam.flagstat: 4061793 

kidney_15.5_R1.mLb.mkD.sorted.bam.flagstat: 3999685
kidney_15.5_R1.mLb.clN.sorted.bam.flagstat: 3542396
kidney_15.5_R2.mLb.mkD.sorted.bam.flagstat: 4576576
kidney_15.5_R2.mLb.clN.sorted.bam.flagstat: 3980732


#### 3. Peak calling


a)Number of peaks
  kidney_14.5_R1.mLb.clN_peaks.count_mqc.tsv: 2002
  kidney_14.5_R2.mLb.clN_peaks.count_mqc.tsv: 3470
  kidney_15.5_R1.mLb.clN_peaks.count_mqc.tsv: 3053
  kidney_15.5_R2.mLb.clN_peaks.count_mqc.tsv: 3329

b) *_peaks files are format files which contain the peak locations together with peak summit, pvalue and qvalue. However, *_summits files contain peak summits locations for every peak in order to find the motifs at the binding sites, this file is recommended


#### 4. MultiQC report
a) Highest Duplication rate
kidney_14.5_R1_T1_1: 15.7%

b) Most peaks called sample
kidney_14.5_R2: 3470

c) Lowesst FRiP score sample
kidney_14.5_R1: 0.1

d) Percentage range of peaks overlapping with promoters
15% - 27.3% (kidney_14.5_R2 - kidney_14.5_R1)

e) MERGED LIB is the section of the report that shows results after merging libraries and before filtering. However, MERGED REP is the section of the report that shows results after merging replicates and filtering.



# 2. Integrative analysis using ChrAccR



### 2.1 Running the pipeline
#### 3. GenomicRanges objects
a) These regions represent the accessibility regions. These are important for chromatin accesibility since they are regulatory elements which can modulate chormatinn accessibility

b)
```
library(ChrAccR)                                                                                               

setConfigElement("annotationColumns", c("sampleId", "tissue", "time","replicate"))

setConfigElement("filteringCovgCount", 1L)
setConfigElement("filteringCovgReqSamples", 0.25)
setConfigElement("filteringSexChroms", TRUE)
setConfigElement("normalizationMethod", "quantile")

diffCompNames <- c(
  "liver vs kidney [tissue]"
)
setConfigElement("differentialCompNames", diffCompNames)

rds_file=readRDS("/vol/COMPEPIWS/data/annotation/regionSetList.rds")   

bedfile <- read.table("/vol/COMPEPIWS/groups/atacseq1/tasks/nextflow_run/results/bwa/mergedLibrary/macs/narrowPeak/consensus/consensus_peaks.mLb.clN.bed", sep="\t", header=FALSE, comment.char="", quote="", stringsAsFactors=FALSE)

gr <- GenomicRanges::GRanges(bedfile[,"V1"], IRanges::IRanges(start=bedfile[,"V2"], end=bedfile[,"V3"], names=bed
gr <- sort(gr)

regionslist <- c(rds_file, list(gr))
names(regionslist)[4]<-"consensusPeaks"

temp=read.table("/vol/COMPEPIWS/groups/atacseq1/tasks/sample_annotation_file.tsv", header=TRUE)

run_atac("ChrAccR_analysis_final_final2", "bamFile", temp, genome="mm10", sampleIdCol="sampleId", regionSets=regionslist)

```


 

### 2.2 Interpreting the output
#### 1. Summary
a)
![](https://i.imgur.com/WjRekgn.png)

b) low numbers of fragments -2 
![](https://i.imgur.com/qvzOz1G.png)

c) It corresponds to +1 nucleosome
d) Yes, kidney_14.5_R1 has the lowest TSS enrichment score of 4.078.

#### 2. Filtering
a) ![](https://i.imgur.com/Co7OkLn.png)
 
These regions represent were removed because they have a coverage of less than 1 and represent regions and fragments on chromosomes chrM, chrX, chrY

#### 3. Normalization
a) quantileNorm
b) Normalization increased the number of insertion counts in these samples as can be seen in the quantile heatmap plot when comparing the heatmap of unnormalized and the normalized heat maps. This is a typical aftermath of normalization as it makes the data more uniformly distributed and easier to interpret.

#### 4. Exploratory analysis
a)best discriminated region 
dimension reduction-  peaks_nf: since the samples are better seperated due to large variance explained by PC1(80.6) and PC2(5.96%) 

clustering plots-gr region: since the distances between the clusters are better well defined while using gr as a discriminater rather than tilling500bp. Furthermore, the colors in the heatmap as well are more differntial and distinguishable while comparing the tissue types.

b) Yes, samples have been separated using PCA based on developmental time(gr region best separeates it based on time)

c)  Gata4, GATA::TAL1, and RAX


#### 5. Differential analysis
a) Peaks in the liver are generally more accessible


b)less accessible in kidney compared to liver       
MA0514.1_Sox3                                        
MA0718.1_RAX                                         


more accessible in kidney compared to liver
MA0634.1_ALX3
MA0521.1_Tcf12
# 3 Exploratory and differential analysis
1.
```
a) 8 samples
b) Fragments per sample
kidney_14.5_R1 kidney_14.5_R2 kidney_15.5_R1 kidney_15.5_R2  liver_14.5_R1
       1560477        1804344        1771188        1990352         907666
 liver_14.5_R2  liver_15.5_R1  liver_15.5_R2
       2261596         875366         875366
`

c) Region per sample
4 region types: tiling, promoter, peaks_nf, consensusPeaks
```
2) 
```
filter_data=loadDsAcc("/vol/COMPEPIWS/groups/atacseq1/tasks/ChrAccR_analysis_final_final2/data/dsATAC_filtered")
transformed_data=transformCounts(filter_data,method="CPM")
peakCounts <- getCounts(transformed_data, "peaks_nf")
isVar <- rank(-matrixStats::rowVars(peakCounts)) <= 100
dsa_varPeaks <- removeRegions(transformed_data, !isVar, "peaks_nf")
getNRegions(dsa_varPeaks, "peaks_nf")
cvRes <- getChromVarDev(dsa_varPeaks, "peaks_nf")
devZ <- chromVAR::deviationScores(cvRes)
pheatmap::pheatmap(devZ[1:100,])
```
![](https://i.imgur.com/NBeScN6.png)




3) peak count= 6380
   Promoter counts =1241
   


```{r}
peakCounts <- getCounts(filter_data, "peaks_nf")
peakCountspro <- getCounts(filter_data, "promoter")
write.table(peakCountspro, file="promoterCounts.tsv", sep="\t")
write.table(peakCounts, file="peakCounts.tsv", sep="\t")    
```

4.
```{r}
cvRes <- getChromVarDev(filter_data, "consensusPeaks", motifs="jaspar_vert")
    
devZ <- chromVAR::deviationScores(cvRes) 
str(devZ)
```
5.
```{r}
library(chromVAR)
varia=computeVariability(cvRes)   
plotVariability(varia, n = 5, use_plotly = FALSE)
```
![](https://i.imgur.com/U6ZwODj.png)
 
Top 5: "JUND(var.2)" "RBPJ"        "SOX10"       "EWSR1-FLI1"  "MZF1(var.2)"

6.

```{r}
pheatmap::pheatmap(devZ[1:20,])
```
![](https://i.imgur.com/1ttj1Bn.png)

7.

```{r}
## a
daTab <- getDiffAcc(filter_data, "peaks_nf", "tissue", grp1Name="kidney", grp2Name="liver")

## b
daTab$diffAcc <- "NO"
daTab$diffAcc[daTab$log2FoldChange > 3 & daTab$padj < 0.01] <- "kidney"
daTab$diffAcc[daTab$log2FoldChange < -3 & daTab$padj < 0.01] <- "liver"
table(daTab$diffAcc)


table(daTab$diffAcc)
more accessibile in liver          not differential    more accessibile in kidney
                     513                      5713                             154    


## c
library(ggplot2)
p=ggplot(data=daTab, aes(x=log2FoldChange, y=-log10(padj),color=diffAcc)) + geom_point()
ggsave("volcano.png",p)

```
![](https://i.imgur.com/CcpZvJC.png)


```{r}
## d
genome_coord=getCoord(filter_data,"peaks_nf")

## e
library(GenomicRanges)
TSS_cord=ChrAccRAnnotationMm10::getGeneAnnotation(anno="gencode_coding",type="tssGr")
TSS_close=distanceToNearest(genome_coord,TSS_cord)

## f
difPeak_ano =as.data.frame(TSS_close)
TSS_merged=cbind(daTab,difPeak_ano)
write.table(TSS_merged, file="difPeak_ano.tsv", sep="\t")

## g
#### Check with tutorsss
p=ggplot(TSS_merged, aes(x=diffAcc, y=distance,fill=diffAcc)) + geom_violin(trim=FALSE) + ylim(NA,100000)
ggsave("violin.png",p)
```
![](https://i.imgur.com/taZnreB.png)




```{r}
## h is optional so skipped for time

## i
genome_coord=as.data.frame(genome_coord)
TSS_merged=cbind(genome_coord,daTab)
kidney_index=vector()
countk=0
liver_index=vector()
countl=0
for(i in 1:nrow(TSS_merged)){
  if(TSS_merged[i,21]=="kidney"){
    countk=countk+1
    kidney_index[countk]=i
  }
  else if(TSS_merged[i,21]=="liver"){
    countl=countl+1
    liver_index[countl]=i
  }
}
kidney_df=TSS_merged[kidney_index,]
liver_df=TSS_merged[liver_index,]
kidney_df=kidney_df[,-c(4:7)]
liver_df=liver_df[,-c(4:7)]
rownames(kidney_df)=NULL
colnames(kidney_df)=NULL
rownames(liver_df)=NULL
colnames(liver_df)=NULL
write.table(kidney_df, "kidney.bed",sep= "\t", row.names=FALSE,quote=FALSE)
write.table(liver_df, "liver.bed", sep= "\t",row.names=FALSE,quote=FALSE)

```


8.
```{r}
#a
dsa_merged <- mergeSamples(filter_data, "tissue")

#b
tiling200bp = muRtools::getTilingRegions("mm10", width=200L, onlyMainChrs=TRUE)
dsa_erb <- regionAggregation(dsa_merged, tiling200bp, "tiling200bp", signal="insertions")

#c
exportCountTracks(dsa_erb,type="tiling",outDir=".",format="igv")
```

9.

```{r}
motifNames <- c("JUND(var.2)","RBPJ","SOX10","EWSR1-FLI1")
samples <- getSamples(dsa_merged)
fps <- getMotifFootprints(dsa_merged, motifNames, samples)
```

10.d) Firstly all the top differential regions are in the liver sample compared to the kidney. Furthermore, there is completely distinct peak signals between the kidney and liver samples when checking the most differential regions which is expected since these were found to be most significant. Finally there are more changes in the distal peaks compared to the promoter regions as was clear when checking in the IGV browser.


![](https://i.imgur.com/QrwZDaa.png)
![](https://i.imgur.com/Makid5W.png)



# 4 Prepare data for downstream integrative analysis


# 1 Share assay-specific results with other groups


# 2 Integrative data exploration in IGV
1) In general the promoter regions are low methylated regions and the gene body regions are highly methylated regions
2) UMR HMR

# 3 Integrative analysis

## 1. Summarize signals over regions of interest
### a) Defining regions of interest (ROI)
```{bash}

 awk '{ OFS="\t"; if ($6=="+"){ print $1,$2-1500,$2+500} else {print $1,$3+1500,$3-500} }' /vol/COMPEPIWS/pipelines/references/mm10_genome_genes.bed> promoter_roi.bed

```


 b) PLOTTING
```{bash}
computeMatrix scale-regions -S  /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K27ac_R1.bigWig\
					  /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K36me3_R1.bigWig\
					  /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K27me3_R1.bigWig\
					  /vol/COMPEPIWS/groups/shared/ChIP-seq/chipseq1/signals/kidney_14.5_H3K4me1_R1.bigWig\
                              -R promoter_roi.bed \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o matrix.mat.gz

plotHeatmap -m matrix.mat.gz \
      -out Heatmap1.png 
      
plotProfile -m matrix.mat.gz \
     -out Profile.png 
   ```   
 ![](https://i.imgur.com/2Q4pS0L.jpg)

 ![](https://i.imgur.com/46FLG2X.png)
     
we can generate more if required for presentation(with different bigwig fies)


## 2. Chromatin states and DNA methylation
```{bash}
a)  
bedtools intersect -wa -wb -a /vol/COMPEPIWS/groups/shared/WGBS/wgbs1/segmentation/kidney_final.bed -b kidney_15_track_dense.bed> /vol/COMPEPIWS/groups/atacseq1/tasks/average_Meth.bed

bedtools intersect -wa -wb -a /vol/COMPEPIWS/groups/shared/WGBS/wgbs1/segmentation/liver_final.bed -b liver_15_track_dense.bed> /vol/COMPEPIWS/groups/atacseq1/tasks/average_LMeth.bed
```

```{r}
b)

avgM= read.table("average_Meth.bed", header=FALSE) 
 p=ggplot(avgM,aes(y=V5,color=V13))+geom_boxplot()+facet_wrap(~V13)

```
![](https://i.imgur.com/mbvO6Xr.png)

## 3. DNA methylation and chromatin accessibility




## 4. Differentially Methylated Regions (DMRs)

## 5. Differential methylation and gene expression
```












