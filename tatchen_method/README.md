# g-matrix analysis toolkit

Python toolkit for analyzing **g-matrix** using Quasi-degenerate perturbation theory within **PyQchem** (https://github.com/abelcarreras/PyQchem), a Python wrapper for Q-Chem (https://www.q-chem.com).

## 📁 Main features
Perform:
* Calculation of the g-shift 
* Sum-over-states (SOS) plots: g-matrix calculation between ground and the different excited states separately
* Scaling analysis: evaluation of the g-matrix when increasing the spin-orbit coupling as a scaling parameter
* Spin contamination analysis for each excited state

## 📁 Project Structure

```yaml 
tatchen_method/
│
├── gmatrix_program/
│ ├── gmatrix_functions.py # Functions used in the different classes
│ ├── gmatrix_classes.py # Classes wrapping the functions into workflows
│
├── examples/
│ └── gmatrix_example.py # Example of how to use the pipeline
│ └── test_1.out # QChem output example to be used in the gmatrix_example.py
│ └── test_2.out # QChem output example to be used in the gmatrix_example.py
│ └── test_2.out # QChem output example to be used in the gmatrix_example.py│
|
├── README.md
└── requirements.txt
``` 

## 🚀 Usage

### 1. Installation

Check https://github.com/abelcarreras/PyQchem. 

### 2. Example
Simple Python API to define your pipeline, defined in *gmatrix_example.py*: 

```python
import sys
from gmatrix_program.gmatrix_classes import OutToJsonConverter, GTensorConfig, GTensorPipeline 

# Path to your Q-Chem output file
qchemout = sys.argv[1] 

# Obtain the json file from the QChem .out file
converter = OutToJsonConverter(qchemout)
converter.read_file()
converter.parse()
converter.save_json()

# Create a configuration object
config = GTensorConfig()

# Create and run the pipeline using this configuration
pipeline = GTensorPipeline(qchemout, config=config)
pipeline.load_data()
pipeline.select_states()
pipeline.run_gshift_calculation()
```
Then in your command line:
```bash
python gmatrix_example.py test1.out
```

### Dependencies
See requirements.txt.
matplotlib for plotting
(Standard libraries like json, sys, and os are used internally and do not need installation.)
