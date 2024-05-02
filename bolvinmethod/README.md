[//]: # ([![Build, Test & Upload]&#40;https://github.com/abelcarreras/PyQchem/actions/workflows/test-publish.yaml/badge.svg&#41;]&#40;https://github.com/abelcarreras/PyQchem/actions/workflows/test-publish.yaml&#41;)

Main feature of the code
=======
Python program to calculate the g-tensor using QChem (https://www.q-chem.com) outputs and 
PyQChem parsers (https://pyqchem.readthedocs.io/en/master/).

Installation requirements
-------------------------
- PyQChem
- numpy
- sys

Files
-------------------------
- [gtensor_doublets.py](gtensor_doublets.py): functions to be used in the [test_doublets.py](test_doublets.py) file. 
- [test_doublets.py](test_doublets.py): file to be run selection of the input file and input values. 

Example 
--------
```console
python test_doublets.py example_doublets.out
```
