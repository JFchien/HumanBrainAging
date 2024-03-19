library(tidyverse)
library(clusterProfiler)
library(magrittr)
library(viridis)
library(data.table)
library(org.Hs.eg.db)

celltypes=c(
            'Exc_up','Inh_up',
            'Exc_down','Inh_down'
            'all_up','all_down'
            )

for (celltype in celltypes){
  bg <-read_tsv(paste0("/cndd2/jchien/project/CZI_human/pseudobulk_rna/dream_result/rm_CZI1_AM2_AM3_AF2A/GO/ageDE_",celltype,"_bg.txt"),col_names='gene')
  de <-read_tsv(paste0("/cndd2/jchien/project/CZI_human/pseudobulk_rna/dream_result/rm_CZI1_AM2_AM3_AF2A/GO/ageDE_",celltype,"_up.txt"),col_names='gene')

  bg_list <- bg$gene
  de_list <- de$gene

  print(celltype)
  print(length(bg_list))
  print(length(de_list))

  set.seed(667)
  ego <- enrichGO(de_list,
                        OrgDb = org.Hs.eg.db,
                        ont = "BP",
                        universe = bg_list,
                        pAdjustMethod = "BH",
                        keyType = "SYMBOL",
                        # pvalueCutoff = 0.01,
                        # qvalueCutoff = 0.05,
                        minGSSize = 25,
                        maxGSSize = 500,
                        readable = TRUE)
  final_res <- ego@result %>% as_tibble()
  write_tsv(final_res, paste0("/cndd2/jchien/project/CZI_human/pseudobulk_rna/dream_result/rm_CZI1_AM2_AM3_AF2A/GO/clusterProfiler/GO_BP_ageDE_",celltype,"_up.tsv.gz"))
  ego_simp <- clusterProfiler::simplify(ego, cutoff=0.5,by="p.adjust", select_fun=min,measure="Wang")
  final_res2 <- ego_simp@result %>% as_tibble()
  write_tsv(final_res2, paste0("/cndd2/jchien/project/CZI_human/pseudobulk_rna/dream_result/rm_CZI1_AM2_AM3_AF2A/GO/clusterProfiler/simp05/GO_BP_ageDE_",celltype,"_up_simp05.tsv.gz"))
}
