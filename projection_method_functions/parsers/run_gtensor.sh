#!/bin/bash

# module load Python/3.7.6-Anaconda3-2020.02
# export MODULEPATH=/scratch/abel/SOFTWARE/privatemodules:$MODULEPATH
# export PYTHONPATH=/dipc/acebg/0_scripts/gtensor:$PYTHONPATH
# module load qchem_group
# module load Python

# Comprueba si se proporcionaron los dos argumentos
if [ $# -ne 1 ]; then
    echo "Check arguments: $0 <archivo>"
    exit 1
fi

archivo=$1

if [ ! -e "$archivo.json" ]; then
    python ~/Desktop/my_programs/gtensor/projection_method/parsers/parser_read_data.py $archivo
fi
    python ~/Desktop/my_programs/gtensor/projection_method/parsers/input.py $archivo