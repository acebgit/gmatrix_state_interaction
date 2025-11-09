import sys

archivo = 'phenalenyl_opt.out'  # str(sys.argv[1])

output_name = archivo[:-4] + '.molecule'
# sys.stdout = open(output_name, "w")

with open(archivo, encoding="utf8") as file:
    for line in file:
        if '**  OPTIMIZATION CONVERGED  **' in line:

            # To take the Standard Nuclear Orientation:
            # while '    1  ' not in line:
            #     line = next(file)
            #
            # while ' -------------------------' not in line:
            #     line_nueva = line.split()
            #     print(line[10:],end='')
            #     line = next(file)

            # To take the z-matrix:
            while 'Z-matrix Print:' not in line:
                line = next(file)

            line = next(file)
            while '$end' not in line:
                print(line, end='')
                line = next(file)
            print(line, end='')
