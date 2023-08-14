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

Examples 
--------
**Simple API to define the input:** An example is shown in file [test.py](test.py)

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
- [test.py](test.py): file to test programs. At the moment, only RAS results. 

- [parser_eom.py](parser_eom.py): generate the input from eom output and calculate eom g-tensor
- [gnt_fe_pyms2_+_eomip_def2-SVPD.in.out.json](gnt_fe_pyms2_+_eomip_def2-SVPD.in.out.json): general input created from eom output

- [fe_pyms2_def2tzvp_ras.out](fe_pyms2_def2tzvp_ras.out): rasci output
- [parser_rasci_gtensor.py](parser_rasci_gtensor.py): calculate rasci g-tensor
- [parser_rasci_pyqchem.py](parser_rasci_pyqchem.py): pyqchem parser used in rasci

Contact info
------------
Anna Krylov anna.i.krylov@gmail.com, University of Southern California (USC)

David Casanova david.casanova@ehu.eus, Donostia International Physics Center (DIPC)