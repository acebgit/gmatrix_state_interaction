#!/bin/bash

# module load Python/3.7.6-Anaconda3-2020.02
# export MODULEPATH=/scratch/abel/SOFTWARE/privatemodules:$MODULEPATH
# export PYTHONPATH=/dipc/acebg/0_scripts/gtensor:$PYTHONPATH
# module load qchem_group
# module load Python

# Comprueba si se proporcionaron los dos argumentos
if [ $# -ne 2 ]; then
    echo "Check arguments: $0 <method> <archivo>"
    exit 1
fi

method=$1
archivo=$2

# Utiliza la opción en una declaración case
case $method in
    "rasci")
        python ~/Desktop/my_programs/gtensor/projection_method/parsers/parser_qchem_rasci.py $archivo
        python ~/Desktop/my_programs/gtensor/projection_method/parsers/gtensor_calculation.py $archivo
        ;;
    "rascifull")
        python ~/Desktop/my_programs/gtensor/projection_method/parsers/parser_qchem_rasci_detailed.py $archivo
        ;;
    "tddft")
        python ~/Desktop/my_programs/gtensor/projection_method/parsers/parser_qchem_tddft.py $archivo
        python ~/Desktop/my_programs/gtensor/projection_method/parsers/gtensor_calculation.py $archivo
        ;;
    "eomcc")
        python ~/Desktop/my_programs/gtensor/projection_method/parsers/parser_qchem_eom.py $archivo
        python ~/Desktop/my_programs/gtensor/projection_method/parsers/gtensor_calculation.py $archivo
        ;;
    *)
        echo "Not valid. Methods available: rasci, rascifull, tddft, eomcc"
        ;;
esac
