#!/bin/bash

# v2.py is the modified version of Oxford's OPIG group
for file in ./pdb_structures/*.pdb; do
    python v2.py --f1 ${file} --f2 ${file} --c1 AB --c2 C --c 4.5 --i 10.0 --jobid .\\$(basename ${file%.*})\\
done

for folder in ./out_/* ; do
    for file in ./${folder}/*.txt; do
        sed 's/[[:space:]]/,/g' ${file} > ${file}.csv
    done
done

for folder in ./out_/* ; do
    for file in ./${folder}/*.txt; do
        find . -type f -name "parameters.txt.csv*" -print -delete

    done
done
