# Magnetic property calculation: g-matrix with state-interaction procedure

**Development of computer codes for calculating the molecular g-matrix using the state-interaction (SI) procedure.**

---

## 🧠 Overview

This repository contains Python programs developed to compute **magnetic g-matrices** through the **state-interaction (SI)** approach.  
The workflow consists of two main post-processing stages:

1. **Parsing Q-Chem outputs** to extract relevant electronic structure data (state energies, spin–orbit couplings, and angular momentum matrix elements).  
2. **Computing the g-matrix** based on the SI effective Hamiltonian formalism.

Electronic structure calculations are performed with a developer version of **Q-Chem 6.0**, which provides all the quantities required by the SI scheme. The evaluation of the g-matrix is carried out with in-house codes integrated within **PyQChem** and **ezMagnet** frameworks.

---

## ⚙️ Methodology Summary

The SI procedure constructs and diagonalizes an **effective Hamiltonian** whose:
- diagonal elements correspond to the **non-relativistic state energies**, and  
- off-diagonal elements contain the **spin–orbit coupling (SOC)** matrix elements.

The resulting **spin–orbit-coupled states** are linear combinations of the non-relativistic states.  
Two implementations are available:
- **Bolvin’s approach**, for systems with doublet multiplicity)  
- **Tatchen’s method**, for arbitrary spin multiplicities

This framework enables a detailed **rationalization of g-shifts**, analyzing how excitation energies, SOCs, and orbital angular momenta jointly determine magnetic anisotropy.

---

## 🚀 Main Features

- Computation of the **molecular g-tensor** from Q-Chem output files.  
- Support for **RAS-CI**, **TDDFT** and **EOM** calculations.  
- JSON-based intermediate representation of extracted electronic data.  
- Modular structure for integration into larger workflows (e.g., PyQChem, ezMagnet).  

---

## 🧰 Installation Requirements

Python ≥ 3.8 and the following packages:
numpy
json
decimal
pymatgen


(Other built-in modules such as `warnings`, `sys`, `hashlib`, and `operator` are used internally.)

---

## 📁 Repository Structure

| File | Description |
|------|--------------|
| `parser_eom.py` | Parses Q-Chem EOM-CC outputs and generates a `.json` input for g-tensor calculation. |
| `parser_rasci.py` | Parses Q-Chem RAS-CI outputs and generates a `.json` input for g-tensor calculation. |
| `gtensor.py` | Computes the g-tensor from the generated `.json` data file. |

---

## 🧪 Example Usage

**Simple command-line workflow:**

```bash
# For EOM outputs
python parser_eom.py example_eom.out
python gtensor.py example_eom.json > example_eom_results.out

# For RAS-CI outputs
python parser_rasci.py example_ras.out
python gtensor.py example_ras.json > example_ras_results.out
```

## 📖 References

* A. Krylov et al., Q-Chem 6.0 (developer version)
* Bolvin, Chem. Phys. Lett. 2006, 428, 142–147.
* Tatchen et al., J. Chem. Phys. 2007, 126, 034303.

## 📬 Contact

Antonio Cebreiro

Donostia International Physics Centre

📧 antonio.cebreiro@dipc.org

## 🧾 License

This project is distributed for research and educational use.
Please cite appropriately if used in scientific publications.
