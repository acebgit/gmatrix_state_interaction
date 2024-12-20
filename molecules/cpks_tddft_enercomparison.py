import sys

def search_and_print_line(file_path, search_word, line_position):
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

if len(sys.argv) != 3:
    print("Usage: python script.py TDDFT.out CPKS.out")
else:
    file1 = sys.argv[1]
    file2 = sys.argv[2]

    tddft_ener = search_and_print_line(file1, "Total energy in the final basis set", -1)
    cpks_ener = search_and_print_line(file2, "Total energy after final integration", -2)
    print("TDDFT energy: {}, CPKS: {}, Difference: {}.".format(tddft_ener, cpks_ener, tddft_ener-cpks_ener))

