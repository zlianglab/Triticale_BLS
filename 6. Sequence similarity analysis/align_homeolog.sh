#!/bin/bash

# Read the homeolog triads
i=1
while read geneA geneB geneR; do
	# Combined fasta file for triplet
	outfile="triplet_${i}.fasta"

	# Extract sequence using seqkit
	seqkit grep -n -p "$geneA" fasta/1000bp/ABR_diff_A.fa > mafft/1000bp/ABR-diff/$outfile
	seqkit grep -n -p "$geneB" fasta/1000bp/ABR_diff_B.fa >> mafft/1000bp/ABR-diff/$outfile
	seqkit grep -n -p "$geneR" fasta/1000bp/ABR_diff_R.fa >> mafft/1000bp/ABR-diff/$outfile

	# Run mafft alignment
	mafft --localpair --maxiterate 1000 --thread 8 \
	mafft/1000bp/ABR-diff/$outfile > mafft/1000bp/ABR-diff/triplet_${i}.aln

	# Clean up the temp file
	rm mafft/1000bp/ABR-diff/$outfile

	i=$((i+1))

done < genelist/ABR-diff.txt
