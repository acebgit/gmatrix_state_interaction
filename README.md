[//]: # ([![Build, Test & Upload]&#40;https://github.com/abelcarreras/PyQchem/actions/workflows/test-publish.yaml/badge.svg&#41;]&#40;https://github.com/abelcarreras/PyQchem/actions/workflows/test-publish.yaml&#41;)

Qchem
=======
Python wrapper for Q-Chem (https://www.q-chem.com)

Main feature of the code
--------------------------
- Calculation of the g-tensor by using Q-Chem outputs. 

Installation requirements
-------------------------
- requests
- operator
- warnings
- numpy
- hashlib
- json
- matlib
- sys
- pymatgen
- decimal

Files
-------------------------
- [parser_eom.py](parser_eom_antiguo.py): generates "json" general input from EOM qchem outputs to be used by the [gtensor.py](gtensor.py) script.
- [parser_rasci.py](parser_rasci.py): generates "json" general input from RAS-CI qchem outputs to be used by the [gtensor.py](gtensor.py) script. 
- [gtensor.py](gtensor.py): calculates g-tensor using "json" input.

Examples 
--------
**Simple command line to obtain the g-tensor:**
```console
# To calculate the EOM output gtensor:
python parser_eom.py example_eom.out 
python gtensor.py example_eom.out > example_eom_results.out

# To calculate the RAS-CI output gtensor:
python parser_eom.py example_ras.out 
python gtensor.py example_ras.out > example_ras_results.out
```

Contact info
------------
Anna Krylov anna.i.krylov@gmail.com, University of Southern California (USC)

David Casanova david.casanova@ehu.eus, Donostia International Physics Center (DIPC)