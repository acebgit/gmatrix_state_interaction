[//]: # ([![Build, Test & Upload]&#40;https://github.com/abelcarreras/PyQchem/actions/workflows/test-publish.yaml/badge.svg&#41;]&#40;https://github.com/abelcarreras/PyQchem/actions/workflows/test-publish.yaml&#41;)

Qchem
=======
Python wrapper for g-tensor calculation in PyQChem (https://pyqchem.readthedocs.io/en/master/)

Main feature of the code
--------------------------
- Calculation of the g-tensor by using Q-Chem outputs. 

Installation requirements
-------------------------
- PyQChem
- numpy
- sys

Files
-------------------------
- [parser_eom.py](parser_eom_antiguo.py): generates "json" general input from EOM qchem outputs to be used by the [gtensor.py](gtensor.py) script.
- [parser_rasci.py](parser_rasci.py): generates "json" general input from RAS-CI qchem outputs to be used by the [gtensor.py](gtensor.py) script. 
- [gtensor.py](gtensor.py): calculates g-tensor using "json" input.