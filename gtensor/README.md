[//]: # ([![Build, Test & Upload]&#40;https://github.com/abelcarreras/PyQchem/actions/workflows/test-publish.yaml/badge.svg&#41;]&#40;https://github.com/abelcarreras/PyQchem/actions/workflows/test-publish.yaml&#41;)

Qchem
=======
Python wrapper for Q-Chem (https://www.q-chem.com)

Main features of the code
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

Examples 
--------
**Simple API to define the input:** An example is shown below (pending)

[//]: # ()
[//]: # (```python)

[//]: # (from pyqchem import Structure, QchemInput, get_output_from_qchem)

[//]: # (from pyqchem.parsers.basic import basic_parser_qchem)

[//]: # ()
[//]: # (molecule = Structure&#40;coordinates=[[0.0, 0.0, 0.0],)

[//]: # (                                  [0.0, 0.0, 0.9]],)

[//]: # (                     symbols=['H', 'H'],)

[//]: # (                     charge=0,)

[//]: # (                     multiplicity=1&#41;)

[//]: # ()
[//]: # (qc_input = QchemInput&#40;molecule,)

[//]: # (                      jobtype='sp',)

[//]: # (                      exchange='hf',)

[//]: # (                      basis='6-31G'&#41;)

[//]: # ()
[//]: # (data = get_output_from_qchem&#40;qc_input,)

[//]: # (                             processors=4,)

[//]: # (                             parser=basic_parser_qchem&#41;)

[//]: # ()
[//]: # (# obtain a python dictionary)

[//]: # (print&#40;'Energy: ', data['scf_energy']&#41;)

[//]: # (```)

Files
-------------------------
- [gtensor.py](gtensor.py): Calculation of the g-tensor from "json" file.
- [parser_eom.py](parser_eom.py): generates the "json" and "xml" file from EOM Q-Chem output.
- [parser_rasci.py](parser_rasci.py): generates the "json" and "xml" file from RAS-CI Q-Chem output.

Contact info
------------
Anna Krylov anna.i.krylov@gmail.com, University of Southern California (USC)

David Casanova david.casanova@ehu.eus, Donostia International Physics Center (DIPC)