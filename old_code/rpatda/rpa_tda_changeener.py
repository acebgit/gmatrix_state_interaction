import sys 
from projection_method.parsers.gtensor_calculation import extract_data_from_json
from projection_method.parsers.parser_read_data import output_json

try:
    json_file = str(sys.argv[1])
    rpa_file = str(sys.argv[2])
except (IndexError):
      print("INPUT ERROR!! Two files are required: TDA json and RPA qchem output")
      exit()

output_file = json_file.replace('tda', 'rpa').replace('.json','')

# Take the json dictionaries 
output_dict = extract_data_from_json(json_file)

# Take the RPA energies and put them in an energy dictionary
start_search = False
energy_dict = {'1': output_dict['energy_dict']['1']} # Ground state energy (HF energy)
state = 1

with open(rpa_file, "r") as infile:
    lines = infile.readlines()

for nline, line in enumerate(lines):
    if "TDDFT Excitation Energies" in line:
            start_search = True
    
    if start_search and "Excited state" in line:
            excit_ener = float(line.split()[-1])
            total_ener = float(lines[nline+1].split()[-2])
            state += 1
            energy_dict.update({state: [total_ener, excit_ener]})

# Replace the TDA energy dictionary by RPA energy dictionary and put them in a new file
output_dict['energy_dict'] = energy_dict
output_json(output_file, output_dict)

