#!/bin/bash

allcpath=$1
out_file=${allcpath#allclist_}
out_file=${out_file/.txt/.allc.tsv.gz}
echo $out_file
allcools merge-allc --allc_paths $allcpath \
                    --output_path $out_file \
                    --chrom_size_path 'hg38.chrom_1-X.sizes' \
                    --cpu 20
