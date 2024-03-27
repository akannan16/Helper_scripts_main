#!/bin/bash

# Store the list of fastq.gz files
files=$(ls *.gz)

# Open a log file for writing
log_file=wc_results.log
exec 3>&1 1>>${log_file}

# Loop through the files and run wc -l on each one
for file in $files; do
    echo "$file"
    gunzip -c $file | awk '{if(NR%4==2) print}' | wc -l
done

# Close the log file
exec 1>&3 3>&-

