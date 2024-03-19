#!/bin/bash

allcpath=$1
out_file=$(basename ${allcpath%.allc_rmSNPs.tsv.gz})
echo $out_file

allcools standardize-allc --allc_path $allcpath \
                          --chrom_size_path /cndd2/jchien/iGenome/hg38_noALT/withchrL/forallc/genome_withchrL_forallc.chrom.sizes

allcools allc-to-region-count --allc_path $allcpath \
                              --output_prefix $out_file \
                              --chrom_size_path /cndd2/jchien/iGenome/hg38_noALT/withchrL/forallc/genome_withchrL_forallc.chrom.sizes \
                              --mc_contexts CGN CHN \
                              --region_bed_paths /cndd2/jchien/iGenome/gencodev37/gencode.v37.annotation.intragenic_filter.bed.gz \
                              --region_bed_names intragenic
