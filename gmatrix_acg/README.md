# g-Tensor Analysis Toolkit

This project provides a Python toolkit for analyzing **g-tensors** and related spin–orbit coupling (SOC) properties from quantum chemistry output files (e.g., Q-Chem).  
It includes:

- **Utility functions** for reading and processing SOC, angular momentum, and transitions data.
- **Object-oriented classes** that organize workflows (g-shift, S² analysis, sum-over-states (SOS), scaling).
- **Example scripts** demonstrating how to use the pipeline.

---

## 📁 Project Structure

myproject/
│
├── myproject/
│ ├── utils.py # Functions for SOC matrices, angular momentum, transitions, etc.
│ ├── models.py # Classes wrapping the functions into workflows
│ ├── init.py # Expose main classes/functions
│
├── examples/
│ └── example_basic.py # Example of how to use the pipeline
│
├── README.md
└── requirements.txt

yaml
Copy code

---

## 🚀 Usage

### 1. Installation

Clone the repository and install the requirements:

```bash
git clone https://github.com/yourusername/g-tensor-toolkit.git
cd g-tensor-toolkit
pip install -r requirements.txt
```

### 2. Example: Run the Pipeline
```python
from myproject.models import GTensorPipeline, GTensorConfig

# Path to your Q-Chem output file
qchemout = "o2_4_3.out"

# Optional: configure state selection
config = GTensorConfig()
# config.state_selection = 1
# config.initial_states = [1, 2, 3, 4]

# Create and run the pipeline
pipeline = GTensorPipeline(qchemout, config=config)
pipeline.load_data()
pipeline.select_states()
pipeline.run_gshift_calculation()
pipeline.run_sos_analysis()
```

### Features
Parse Q-Chem output to extract:
Spin–orbit coupling matrices

Occupation numbers

Angular momentum elements

Transitions data

Perform:

g-shift analysis

S² analysis

Sum-over-states (SOS) plots

Scaling analysis

Object-oriented design for clean workflows

Easy plotting with matplotlib

### Dependencies
See requirements.txt.
Core libraries:
numpy for linear algebra
matplotlib for plotting
(Standard libraries like json, sys, and os are used internally and do not need installation.)
