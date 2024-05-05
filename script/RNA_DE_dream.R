##https://www.bioconductor.org/packages/devel/bioc/vignettes/variancePartition/inst/doc/dream.html
library(variancePartition)
library(edgeR)
library(BiocParallel)
library(tidyverse)
library(Matrix)
param = MulticoreParam(workers = 16)

genemeta=read_tsv('/cndd2/jchien/iGenome/gencodev37/gencode.v37.annotation.intragenic.bed.gz',
                  col_names = c('chr','start','end','id','strand','description','name','func'))%>% 
  filter(chr!=c('chrY') & chr!=c('chrM'))%>% 
  filter(str_detect(func, 'lncRNA|protein_coding'))
samples=c('YM1B', 'YM2A', 'YM2B', 'YM3A', 'YM3B', 'YF1A', 'YF1B', 'YF2A', 'YF2B', 
 'AM1A', 'AM1B', 'AF1A', 'AF1B', 'AF2B', 'AF3A', 'AF3B')# 'AF2A',
files=list.files('/cndd2/jchien/project/CZI_human/pseudobulk_rna', pattern='.tsv.gz' ,full.names=TRUE)
# celltypes=c(
#   'CGE_ADARB2_ADAM33'
# )
for (fileoi in files){
#for (celltype in celltypes){
  #fileoi=paste0('/cndd2/jchien/project/CZI_human/pseudobulk_rna/',celltype,'.tsv.gz')    	
  print(fileoi)
  counts=read.delim(fileoi,row.names=1)
  print(dim(counts))
  goi=Reduce(intersect,list(genemeta$id,rownames(counts)))
  soi=Reduce(intersect,list(samples,colnames(counts)))
  counts=counts[goi,soi]
  print(dim(counts))
  keep.exprs <-rowSums(cpm(counts)>5) >= 2
  dge = DGEList( counts[keep.exprs,] )
  dge = calcNormFactors( dge )
  print(dim(dge))
  
  # metadata
  samplenames <- colnames(dge)
  Donor <- substr(samplenames,1,3)
  Age <- substr(samplenames,1,1)
  Sex <- substr(samplenames,2,2)
  metadata=data.frame(sample=samplenames,donor=Donor,age=Age,sex=Sex)
  rownames(metadata)=metadata$sample
  
  
  # # age
  design <- ~ 0 + age + sex + (1|donor)
  L = makeContrastsDream( design, metadata, 
                          contrasts = c("ageA - ageY"))
  plotContrasts(L) 
  vobjDream = voomWithDreamWeights( dge, design, metadata, BPPARAM=param )
  fit = dream( vobjDream, design, metadata, L)
  fit = eBayes(fit)
  topTable( fit, coef='ageA - ageY', number=3 )
  top.table <- topTable(fit, coef='ageA - ageY',sort.by = "P", n = Inf)
  print(length(which(top.table$adj.P.Val < 0.05)))
  top.table$Gene <- rownames(top.table)
  top.table <- top.table[,c("Gene", names(top.table)[1:6])]
  outfile=str_replace(fileoi,'pseudobulk_rna/','pseudobulk_rna/dream_result/rm_CZI1_AM2_AM3_AF2A/ageDE_rmRINPMI_')
  write_tsv(top.table,outfile)
  
  # sex
  design <- ~ 0 + sex + age + (1|donor)
  L = makeContrastsDream( design, metadata, 
                          contrasts = c("sexF - sexM"))
  plotContrasts(L) 
  vobjDream = voomWithDreamWeights( dge, design, metadata, BPPARAM=param )
  fit = dream( vobjDream, design, metadata, L)
  fit = eBayes(fit)
  topTable( fit, coef='sexF - sexM', number=3 )
  top.table <- topTable(fit, coef='sexF - sexM',sort.by = "P", n = Inf)
  print(length(which(top.table$adj.P.Val < 0.05)))
  top.table$Gene <- rownames(top.table)
  top.table <- top.table[,c("Gene", names(top.table)[1:6])]
  print(top.table['ENSG00000229807.12',])
  outfile=str_replace(fileoi,'pseudobulk_rna/','pseudobulk_rna/dream_result/rm_CZI1_AM2_AM3_AF2A/sexDE_rmRINPMI_')
  write_tsv(top.table,outfile)

  ## variance partition
  #design <- ~ sex + age + (1|donor)
  #vobjDream = voomWithDreamWeights( dge, design, metadata, BPPARAM=param )
  #form = ~ (1|age) + (1|sex) + (1|donor)
  #vp = fitExtractVarPartModel( vobjDream, form, metadata)
  ## plotVarPart( sortCols(vp))
  #vp$Gene <- rownames(vp)
  #vp <- vp[,c("Gene", names(vp)[1:4])]
  #outfile=str_replace(fileoi,'pseudobulk_rna/','pseudobulk_rna/dream_result/rm_CZI1_AM2_AM3/vp_')
  #write_tsv(vp,outfile)
}
