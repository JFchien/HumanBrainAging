library(tidyverse)
library(Seurat)
library(Augur)
library(tictoc)
Sys.setenv(VROOM_CONNECTION_SIZE=131072*8)

meta=read_csv('/cndd2/jchien/project/CZI_human/Augur/metadata.csv.gz') %>%
        column_to_rownames('...1')
expr=read_csv('/cndd2/jchien/project/CZI_human/Augur/cellxgene_log2cpm.csv.gz') %>%
        column_to_rownames('...1') %>% data.matrix()

selected_n_subsamples <- 50
selected_subsample_size <- 400
selected_folds <- 3
selected_var_quantile <- 0.1
selected_feature_perc <- 0.1

auc_out_1 <- calculate_auc(
  expr,
  meta =meta,
  cell_type_col = "level2",
  label_col = "age",
  n_subsamples = selected_n_subsamples, # default 50, the number of random subsamples of fixed size for each cell type
  subsample_size = selected_subsample_size, # default 20, the number of cells per type from each experimental condition
  folds = selected_folds, # default 3
  min_cells = NULL,
  var_quantile = selected_var_quantile,
  feature_perc = selected_feature_perc,
  n_threads = 8,
  show_progress = TRUE,
  classifier = "rf",
  rf_params = list(trees = 100, mtry = 20, min_n = NULL, importance = "accuracy") # mtry = 2
)
saveRDS(auc_out_1, paste0("/cndd2/jchien/project/CZI_human/Augur/nSubsamples_", selected_n_subsamples, "_subsampleSize_", selected_subsample_size, "_folds_", selected_folds,'_',selected_feature_perc,".rds"))
write_tsv(auc_out_1$AUC, paste0("/cndd2/jchien/project/CZI_human/Augur/nSubsamples_", selected_n_subsamples, "_subsampleSize_", selected_subsample_size, "_folds_", selected_folds,'_',selected_feature_perc,".txt"))
