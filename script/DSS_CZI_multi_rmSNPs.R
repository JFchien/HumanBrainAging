library(tidyverse)
library(DSS)
require(bsseq)


datadir="/cndd2/jchien/project/CZI_human/pseudobulk_allc/merge_CG/tab/"
dmrdir="/cndd2/jchien/project/CZI_human/pseudobulk_allc/merge_CG/DMR/"
agedir="/cndd2/jchien/project/CZI_human/pseudobulk_allc/merge_CG/DMR/ageDMR/"
sexdir="/cndd2/jchien/project/CZI_human/pseudobulk_allc/merge_CG/DMR/sexDMR/"
Donor=c('YF1','YF2','YM1','YM2','YM3','AF1','AF2','AF3','AM1','AM2','AM3')
Age = substr(Donor,1,1)
Sex = substr(Donor,2,2)
design = data.frame(Age,Sex)
celltypes=c(
  'L2-4IT_CUX2',
  'L3-5IT_RORB_PLCH1','L4-5IT_RORB_TSHZ2','L4-5IT_RORB_ARHGAP15','L4-5IT_RORB_LRRK1', 
  'L6IT_THEMIS_LINC00343','L6IT_THEMIS_CUX1', 'L56NP_TLE4_TSHZ2','L6b_TLE4_NXPH4','L6CT_TLE4_FAM95C',
  'CGE_LAMP5','CGE_LAMP5_LHX6','CGE_PAX6','CGE_VIP','CGE_ADARB2_ADAM33',
  'MGE_SST','MGE_SST_CLMP','MGE_PVALB'
)
for (celltype in celltypes){
  rm(BSobj)
  message(sprintf("Handling : %s", paste0(dmrdir,"BSobj_",celltype,".rda")))
  load(paste0(dmrdir,"BSobj_",celltype,".rda"))
  # data.list <- vector("list", length(Donor))
  # for (i in 1:length(Donor)) {
  #   message(sprintf("Handling : %s", paste0("tab_",celltype, "-", Donor[i],".CGN-Merge.allc.tsv.gz")))
  #   toDoFile <- paste0(datadir, "tab_",celltype, "-", Donor[i],".CGN-Merge.allc.tsv.gz")
  #   data.list[[i]] <- read.table(toDoFile, header=TRUE)
  # }

  # BSobj <- makeBSseqData(data.list,Donor)
  # save(BSobj, file = paste0(dmrdir,"BSobj_",celltype,".rda"))
  
  # Age
  message(sprintf("Age DMLfit.multiFactor..."))
  DMLfit = DMLfit.multiFactor(BSobj, design=design, formula=~Age+Sex, smoothing=TRUE, smoothing.span=500)
  dml_multi = DMLtest.multiFactor(DMLfit, coef='AgeY') #DMLfit$X
  dml_multi$stat=-dml_multi$stat
  save(dml_multi, file = paste0(agedir,"DMLfit-multiFactor-Age+Sex_AvsY_",celltype,"_smoothed.rda"))
  
  pthres=0.0005
  message(sprintf("Age dmr_multi..."))
  dmr_multi = callDMR(dml_multi, p.threshold=pthres)
  dmr_multi.df <- as_data_frame(dmr_multi)
  write_tsv(dmr_multi.df, paste0(agedir,"DMR-Age+Sex_AvsY_",celltype,"_smoothed_p",pthres,".tsv.gz"))
  rm(dml_multi)
  rm(dmr_multi)
  rm(dmr_multi.df)
  
  # Sex
  message(sprintf("Sex DMLfit.multiFactor..."))
  DMLfit = DMLfit.multiFactor(BSobj, design=design, formula= ~ Sex + Age, smoothing=TRUE, smoothing.span=500)
  dml_multi = DMLtest.multiFactor(DMLfit, coef='SexM') #DMLfit$X
  dml_multi$stat=-dml_multi$stat
  save(dml_multi, file = paste0(sexdir,"DMLfit-multiFactor-Sex+Age_FvsM_",celltype,"_smoothed.rda"))
  
  pthres=0.0005
  message(sprintf("Sex dmr_multi..."))
  dmr_multi = callDMR(dml_multi, p.threshold=pthres)
  dmr_multi.df <- as_data_frame(dmr_multi)
  write_tsv(dmr_multi.df, paste0(sexdir,"DMR-Sex+Age_FvsM_",celltype,"_smoothed_p",pthres,".tsv.gz"))
  rm(dml_multi)
  rm(dmr_multi)
  rm(dmr_multi.df)
}


