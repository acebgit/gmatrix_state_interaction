import sys 

# Name of the file you want to read
input_name = str(sys.argv[1])
output_name = input_name.replace(".out", ".inp")

with open(input_name, 'r') as file:
    with open(output_name, 'w') as output_file:
        # Flag to indicate if "User input:" has been found
        user_input_found = False

        # Read the file line by line
        for line in file:
            # Check if the line contains "User input:"
            if "User input:" in line:
                user_input_found = True
                next(file)

            # Check if "User input:" has been found and print subsequent lines
            elif user_input_found:
                if "-----------------------------------------" in line:
                    break
                else:
                    output_file.write(line)
