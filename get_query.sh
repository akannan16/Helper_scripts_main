#!/bin/bash

# Parse command line arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <reference.fasta> <input.bed>"
    exit 1
fi

ref=$1
bed=$2

# Create a temporary file to store the output of the awk command
tmpfile=$(mktemp)

# Run the awk command and store its output in the temporary file
#awk '{print $1,$2,$3}' $bed > $tmpfile

# Use the temporary file as input to the while loop to extract query sequences for each interval in the BED file
st=/Users/abhiks/Desktop/samtools-1.15.1/samtools
while read -r chrom start end name score strand; do
    $st faidx -n 120 $ref ${chrom}:${start}-${end} | grep -v '>' | sed 's/$/\n/' | sed '/^$/d'
done < $bed > sequences.txt

# Remove the temporary file
#rm "$tmpfile"

paste $bed sequences.txt | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5}' > output.txt
