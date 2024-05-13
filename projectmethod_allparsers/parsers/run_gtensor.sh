#!/bin/bash

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
        python parser_qchem_rasci.py $archivo
        python gtensor_calculation.py $archivo
        ;;
    "tddft")
        python parser_qchem_tddft.py $archivo
        python gtensor_calculation.py $archivo
        ;;
    "eomcc")
        python parser_qchem_eom.py $archivo
        python gtensor_calculation.py $archivo
        ;;
    *)
        echo "Not valid"
        ;;
esac
