import sys
import re

def correct_excitation_energies(filename):
    hartree_to_ev = 27.2114

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Step i: Find all SCF energies and use the last one
    scf_energies = [
        float(match.group(1))
        for line in lines
        if (match := re.search(r"SCF\s+energy in the final basis set\s*=\s*(-?\d+\.\d+)", line))
    ]

    if not scf_energies:
        print("No SCF energy found in the file.")
        return

    scf_energy = scf_energies[-1]  # Use the second (or last) SCF energy

    output_lines = []
    i = 0
    while i < len(lines):
        line = lines[i]

        # Step ii: Look for TDDFT/TDA block
        if "TDDFT/TDA Excitation Energies" in line:
            output_lines.append(line)
            i += 1
            while i < len(lines):
                line = lines[i]

                # Step iii: Replace excitation energy if it is *********
                if "excitation energy (eV) = *********" in line:
                    print(line)
                    exit()
                    prefix = line.split("= *********")[0] + "= "

                    if i + 1 < len(lines):
                        next_line = lines[i + 1]
                        match = re.search(r"Total energy for state\s+\d+:\s+(-?\d+\.\d+)", next_line)
                        if match:
                            excited_energy = float(match.group(1))
                            delta_e_ev = (excited_energy - scf_energy) * hartree_to_ev
                            line = f"{prefix}{delta_e_ev:.6f}\n"

                output_lines.append(line)
                i += 1

                # Optional: stop if another block starts
                if re.match(r"\s*-{5,}", line):
                    break
            continue

        output_lines.append(line)
        i += 1

    # Save output
    output_file = f"corrected_{filename}"
    with open(output_file, "w") as f:
        f.writelines(output_lines)

    print(f"Corrected file saved as: {output_file}")

# Entry point
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py input_file")
        sys.exit(1)

    correct_excitation_energies(sys.argv[1])
