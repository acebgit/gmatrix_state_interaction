#!/bin/bash

#################################################
##   Obtain g-tensor of several files
#################################################

for file in "$@"
do
    python test.py $file
done

#module load Python/3.7.4-Anaconda3-2019.10
#molecule="cucl4_2-"
##reference="d10"
##basis="6-31Gd def2-TZVP"
#
#for file_ms_notnull in $molecule/$molecule*/*.out; do
# echo $file_ms_notnull | cut -f3 -d "/"
# python g_main.py $file_ms_notnull
#done

