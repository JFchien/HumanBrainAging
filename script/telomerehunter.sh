#!/bin/bash

## test telomerehunter
## activate the `telomerehunter` conda environment before running this script
# conda activate telomerehunter

## logging which server and which conda environment this scirpt is running in
echo "Running this on ### $(hostname) ### \(* w *)/"
echo "Currently in conda environment: ### $CONDA_DEFAULT_ENV ###, located in: $CONDA_PREFIX"
echo "Starting at $(date)"
echo "telomerehunter path: $(command -v telomerehunter)"

for cell_id in $(cat cell_id_list_batch_00)
do
	echo "Processing ${cell_id}..."
	telomerehunter \
		-ibc "./bam/${cell_id}.hisat3n_dna.rm_rna_reads.deduped.bam" \
		-b cytoBand_hg38_filtered.tsv \
		-o ./telomerehunter_results \
		-p "${cell_id}" \
		-d \
		-pl
done

echo "Ending at $(date)"
