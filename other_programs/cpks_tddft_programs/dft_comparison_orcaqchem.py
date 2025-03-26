import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re


def take_reference_energy(file_path, search_word, line_position):
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if search_word in line:
                    energy = float((line.split()[line_position]))
                    return energy
        print(f"'{search_word}' not found in {file_path}")
    except FileNotFoundError:
        print(f"File {file_path} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")
    exit()


def comparison_dft_references(qchem_file, orca_file):
    # Take reference energy differences
    molecule = qchem_file.split("_", 1)[0]
    tddft_ener = take_reference_energy(qchem_file, "Total energy in the final basis set", -1)
    cpks_ener = take_reference_energy(orca_file, "Total Energy       :", -4)

    print("--------------------------")
    print("DFT REFERENCE ENERGIES")
    print("--------------------------")
    print("molecule & TDDFT & CPKS & Difference & Relative Difference \\\\ \\hline")
    print(f"\ce{{{molecule}}} & {tddft_ener:.5f} & {cpks_ener:.5f} & {tddft_ener - cpks_ener:.5f} & {(tddft_ener-cpks_ener)*100/tddft_ener:.5f} \\\\")
    compare_s2_value(qchem_file, orca_file)
    print("")
   

def take_qchem_orbitals(file_path, alpha_or_beta):
    """
    Extract even lines between "Alpha MOs, Unrestricted" and "-- Virtual --" from a given file.
    """
    if alpha_or_beta == 'alpha':
        search_word = "Alpha MOs"
    elif alpha_or_beta == 'beta':
        search_word = "Beta MOs"

    lines_to_store = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        # Initialize flags
        found_mo_section = False
        capture = False
        relevant_lines = []

        for line in lines:
            if search_word in line:
                found_mo_section = True
                continue  # Skip this line

            if found_mo_section and "-- Occupied --" in line:
                capture = True
                continue  # Skip this line

            if capture:
                if "-- Virtual --" in line:
                    break  # Stop capturing
                relevant_lines.append(line.strip())  # Store the line

    # Convert list to a single string
    data_str = ' '.join(relevant_lines).replace('********', '0')

    # Use regex to split correctly, preserving negative signs
    numbers = re.findall(r'-?\d+\.\d+|0', data_str)

    # Convert to a list of floats
    occupied_orbital_energ = [float(x) for x in numbers]

    print("If the orbital energy is 0, it means there has been an error in the output format.")
    return occupied_orbital_energ


def extract_spin_up_orbitals(file_path, alpha_or_beta):
    """
    Extract lines from a file starting after the "SPIN UP ORBITALS" line
    until a blank line is encountered.
    
    Args:
        file_path (str): Path to the file.
    
    Returns:
        list: A list of extracted lines.
    """
    if alpha_or_beta == 'alpha':
        search_word = "SPIN UP ORBITALS"
    elif alpha_or_beta == 'beta':
        search_word = "SPIN DOWN ORBITALS"

    lines_to_return = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        # Locate "SPIN UP ORBITALS"
        start_index = None
        for i, line in enumerate(lines):
            if "SPIN UP ORBITALS" in line:
                start_index = i + 2  # Jump one line after "SPIN UP ORBITALS"
                break
        
        # If "SPIN UP ORBITALS" not found, return an empty list
        if start_index is None:
            return lines_to_return

        # Extract lines until a blank line
        for line in lines[start_index:]:
            if line.strip() == "":  # Stop at the first blank line
                break
            lines_to_return.append(line.strip())

    # Extract the third number from each element
    orbital_ener = [float(line.split()[2]) for line in lines_to_return]
    return orbital_ener


def obtain_energ_difference(orca, qchem):
    # Truncate orca_alpha to the length of qchem_alpha
    orca_truncated = orca[:len(qchem)]

    # Create a pandas DataFrame
    df = pd.DataFrame({
        "QChem": qchem,
        "ORCA": orca_truncated
    })

    # Calculate the difference
    df["difference"] = df["QChem"] - df["ORCA"]

    # Create a new column 'Orbitals' with row numbers starting from 1
    df['Orbitals'] = df.index + 1

    # Reorder columns to have 'Orbitals' as the first column
    df = df[['Orbitals'] + [col for col in df.columns if col != 'Orbitals']]
    return df


def plot_difference_vs_orbitals(df, spin='alpha', save_as_png=False, filename='plot.png'):
    titles = {
    'alpha': 'Energy difference in ALPHA orbitals (QChem - ORCA)',
    'beta': 'Energy difference in BETA orbitals (QChem - ORCA)'
    }

    # Ensure the 'Orbitals' column exists and is at the correct position
    if 'Orbitals' not in df.columns or 'difference' not in df.columns:
        raise ValueError("DataFrame must contain 'Orbitals' and 'difference' columns")
    
    # Plotting as a bar plot
    plt.figure(figsize=(8, 6))
    plt.bar(df['Orbitals'], df['difference'], color='b', edgecolor='black')
    
    # Adding labels and title
    plt.xlabel('Orbital Number')
    plt.ylabel(r'$\Delta E$ (QChem - ORCA)')
    plt.title(titles.get(spin, 'Energy difference in UNKNOWN spin orbitals'))

    # Display grid and show the plot
    plt.grid(True)

    # Optionally save the plot as PNG
    if save_as_png:
        plt.savefig(filename, format='png')
    else:
        plt.show()


def comparison_dft_orbitals(qchem_file, orca_file):
    molecule = qchem_file.split("_", 1)[0]

    # Extract relevant even lines from the first input file
    qchem_alpha = take_qchem_orbitals(qchem_file, 'alpha')
    qchem_beta = take_qchem_orbitals(qchem_file, 'beta')

    # Handle second input (example usage)
    orca_alpha = extract_spin_up_orbitals(orca_file, 'alpha')
    orca_beta = extract_spin_up_orbitals(orca_file, 'beta')

    # Print results
    print("--------------------------")
    print("ALPHA ORBITAL ENERGIES")
    print("--------------------------")
    alpha_results = obtain_energ_difference(orca_alpha, qchem_alpha)
    print(alpha_results.to_string(index=False))
    print("")

    print("--------------------------")
    print("BETA ORBITAL ENERGIES")
    print("--------------------------")
    beta_results = obtain_energ_difference(orca_beta, qchem_beta)
    print(beta_results.to_string(index=False))
    print("")
    
    # Call the function to plot
    plot_difference_vs_orbitals(alpha_results, spin='alpha', save_as_png=True, filename=molecule+'_alpha_orbitals.png')
    plot_difference_vs_orbitals(beta_results, spin='beta', save_as_png=True, filename=molecule+'_beta_orbitals.png')


def compare_s2_value(qchem_output, orca_output):
    """
    Extracts the <S^2> value from an output file.

    Parameters:
    file_path (str): Path to the output file.

    Returns:
    float: The extracted <S^2> value, or None if not found.
    """
    with open(qchem_output, 'r') as file:
        for line in file:
            if "<S^2>" in line:  # Find the line containing "<S^2>"
                qchem_s2 = float(line.split()[-1])  # Extract last value as float
    
    with open(orca_output, 'r') as file:
        for line in file:
            if "Expectation value of <S**2>" in line:  # Find the line containing "<S^2>"
                orca_s2 = float(line.split()[-1])  # Extract last value as float
    print("QChem S^2:", qchem_s2, ", ORCA S^2:", orca_s2)


def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <qchem_output> <orca_output> <option>")
        print("Please provide an option: --energy or --orbital")
        sys.exit(1)
    
    qchem_file = sys.argv[1]
    orca_file = sys.argv[2]
    print(f"QChem input: '{qchem_file}'")
    print(f"ORCA input: '{orca_file}'")
    print("")

    # Two options: compare only DFT energies or also orbital energies
    args = sys.argv[1:]

    if '--energy' in args:
        comparison_dft_references(qchem_file, orca_file)
    if '--orbital' in args:
        comparison_dft_references(qchem_file, orca_file)
        comparison_dft_orbitals(qchem_file, orca_file)
    if not args:
        print("Please provide an option: --energy or --orbital")


if __name__ == "__main__":
    main()
