library(tidyverse)
library(MethylSeekR)
library(BSgenome.Hsapiens.UCSC.hg38)


genome_length <- seqlengths(Hsapiens)
CpGislands_gr <-import('/cndd2/jchien/project/CZI_human/DMV/cpgIslandExt_hg38.bed')
CpGislands_gr <- suppressWarnings(resize(CpGislands_gr, 5000, fix="center"))
CpGislands_gr <- keepStandardChromosomes(CpGislands_gr, pruning.mode = "coarse")

find_UMR_LMR <- function(sample_name){
  set.seed(666)
  message(paste0("Analazing sample:", sample_name))
  ## read in data as tibble, and convert into Granges
  meth_df <- read_tsv(paste0("cg_merged_strands_", sample_name, ".txt.gz"),
                      col_types = "cicii",
                      col_names = c("seqnames", "start", "strand", "M", "T"))

  meth_gr <- makeGRangesFromDataFrame(meth_df %>% mutate(end = start, strand = "*") %>% select(seqnames, start, end, strand, T, M),
                                      keep.extra.columns = TRUE,
                                      starts.in.df.are.0based = FALSE) %>%
    dropSeqlevels("L", pruning.mode = "coarse")
  seqlevelsStyle(meth_gr) <- "UCSC"
  seqlevels(meth_gr) <- paste0("chr", c(1:22))
  seqinfo(meth_gr) <- seqinfo(Hsapiens)
  meth_gr <- keepStandardChromosomes(meth_gr, pruning.mode = "coarse")


  stats <- calculateFDRs(m = meth_gr,
                         CGIs = CpGislands_gr,
                         # PMDs = PMDsegments_gr,
                         num.cores = 8,
                         nCpG.cutoffs = seq(1, 8, by = 1),
                         pdfFilename = paste0("meth_threshold_vs_FDR_and_CpG_counts_", sample_name, ".pdf"))
  print(stats)

  # smooth methylation levels over 3 consecutive CpGs to reduce the sampling noise and then identify hypo-methylated regions as
  #  regions of consecutive CpGs with smoothed methylation levels below a fixed cut-off m, containing a minimal number of CpGs n
  FDR_cutoff <- 5 #5%
  m_sel <- 0.5 # mCG < 0.5
  n_sel <- 10 # num. of CG sites in the region
  # n_sel <- as.integer(names(stats$FDRs[as.character(m_sel), ] [stats$FDRs[as.character(m_sel), ]<FDR_cutoff])[1])
  message(paste0("selected n: ", n_sel))

  UMRLMRsegments_gr <- segmentUMRsLMRs(m = meth_gr,
                                       meth.cutoff = m_sel,
                                       nCpG.cutoff = n_sel,
                                       # PMDs = PMDsegments_gr,
                                       num.cores = 8,
                                       myGenomeSeq = Hsapiens,
                                       seqLengths = genome_length,
                                       pdfFilename = paste0("num_CpG_vs_median_meth_", sample_name,
                                                            "_fdr", FDR_cutoff,  "_m", m_sel, "_n", n_sel, ".pdf")
                                       )
  head(UMRLMRsegments_gr)

  plotFinalSegmentation(m = meth_gr,
                        segs = UMRLMRsegments_gr,
                        # PMDs = PMDsegments_gr,
                        meth.cutoff = m_sel,
                        pdfFilename = paste0("example_segmentation_", sample_name,
                                             "_fdr", FDR_cutoff,  "_m", m_sel, "_n", n_sel, ".pdf")
                        )

  UMRLMRsegments_df <- tibble(
    chrom = seqnames(UMRLMRsegments_gr) %>% as.vector(),
    start = start(UMRLMRsegments_gr) - 1L,
    end = end(UMRLMRsegments_gr),
    strand = strand(UMRLMRsegments_gr) %>% as.vector(),
    nCG_covered_by_5more_reads = UMRLMRsegments_gr$nCG.segmentation %>% as.vector(),
    nCG_in_region = UMRLMRsegments_gr$nCG,
    mc = UMRLMRsegments_gr$M %>% as.vector(),
    c = UMRLMRsegments_gr$T %>% as.vector(),
    mean_mc_level = UMRLMRsegments_gr$pmeth  %>% as.vector(),
    median_mc_level_by_each_site = UMRLMRsegments_gr$median.meth,
    region_type = UMRLMRsegments_gr$type
  )

  write_tsv(UMRLMRsegments_df, paste0("UMR_LMR_segments_", sample_name, "_fdr", FDR_cutoff,  "_m", m_sel, "_n", n_sel, ".tsv.gz"))
}

celltypes=c(
  'L2-4IT_CUX2',
  'L3-5IT_RORB_PLCH1',
  'L4-5IT_RORB_TSHZ2',
  'L4-5IT_RORB_ARHGAP15','L4-5IT_RORB_LRRK1',
  'L6IT_THEMIS_LINC00343','L6IT_THEMIS_CUX1', 'L56NP_TLE4_TSHZ2','L6b_TLE4_NXPH4','L6CT_TLE4_FAM95C',
  'CGE_LAMP5','CGE_LAMP5_LHX6','CGE_PAX6','CGE_VIP','CGE_ADARB2_ADAM33',
  'MGE_SST','MGE_SST_CLMP','MGE_PVALB'
)
walk(celltypes, find_UMR_LMR)
