#!/bin/bash


#1 Check contacts
#for file in ./pdb_str_2/*.pdb; do
   #python v2.py --f1 ${file} --f2 ${file} --c1 AB --c2 C --c 4.5 --i 10.0 --jobid .\\$(basename ${file%.*})\\
#done


#2 delete empty folders for which #1 did not work
#for folder in ./out_/*; do
#    find . -empty -type d -delete
#done

#3 delete pdb files in pdb_Str_2 for which #1 did now work

#4 Extract sequence in FASTA format from pdb
#for file in ./pdb_str_2/*.pdb; do
#    pdb_fetch $(basename ${file%.*}) | pdb_selchain -C | pdb_tofasta> ${file}.fasta
#done

#5 Merge all FASTAs together
#cd pdb_str_2
#tail -n *.fasta > merged.fasta


#6 Get csvs with coordinates
#for folder in ./out_/* ; do
 #   for file in ./${folder}/*.txt; do
  #      sed 's/[[:space:]]/,/g' ${file} > ${file}.csv
    #done
#done

#7 Delete parameters.txt.csv as we do not need them
#for folder in ./out_/* ; do
    #for file in ./${folder}/*.txt; do
    #    find . -type f -name "parameters.txt.csv*" -print -delete
 #   done
#done

#run ls -d */ > pdb_names.txt


#run ls -d */ > pdb_names.txt

#sed  's/\\\\//g' pdb_names.txt > _pdb_names.txt




