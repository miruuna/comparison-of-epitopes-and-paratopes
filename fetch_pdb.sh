#!/usr/bin/env bash


for pdb_file in 1fbl 1fdl 1g7h 1ic4 1jhl 1jto 1kip 1mlc 1nby 1ndg 1op9 1pdc 1ri8 1sq2 1t6v 1ua6 1vfb 1xgr 1yqv 1ua6 1xgq 1yqv 1zmy 1zv5 2dqc 2hfm 2i25 2i26 2iff 2yss 2znw 2znx 3eba 3hfm 4i0c 4n1e 4pgj 4ttd 5vjo; do
    pdb_fetch ${pdb_file} > ${pdb_file}.pdb
done
