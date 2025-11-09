##############################################
#   Preparing the run_qchem.sh only
# with 'python run_general.py files partition'
##############################################
# web: https://pythonexamples.org/python-replace-string-in-file/

import sys

# Replace partition, archivo and time required in the same run files
file = 'run_qchem.sh'
partition = str(sys.argv[2])
archivo = str(sys.argv[1])

archivo = archivo.replace('.inp', '')
time = {'regular': '1-00:00:00',
              'test': '00:10:00',
              'long': '2-00:00:00',
              'xlong': '8-00:00:00',
              'large': '2-00:00:00',
              'xlarge': '2-00:00:00'}

partition_line = '#SBATCH --partition='
archivo_line = '#SBATCH --job-name='
error_line = '#SBATCH --error='
time_line = '#SBATCH --time='

with open(file, "r") as f: # Read all lines
    lines = f.readlines()
with open(file, "w") as f: # Rewrite all lines except partition_line
    for line in lines:
        # if line.strip("\n") != partition_line:
        #     f.write(line)
        if partition_line in line:
            new_line = partition_line + partition + "\n"
            f.write(new_line)
        elif archivo_line in line:
            new_line = archivo_line + archivo + "\n"
            f.write(new_line)
        elif error_line in line:
            new_line = error_line + archivo + '.err' + "\n"
            f.write(new_line)
        elif time_line in line:
            new_line = time_line + time[partition] + "\n"
            f.write(new_line)

        else:
            f.write(line)

# Replace partition, archivo and time required in another files
# file_old = 'run_qchem.sh'
# file_new = 'run_qchem.sh'
#
# partition = 'regular' # str(sys.argv[2])
# time = {'regular': '1-00:00:00',
#               'test': '00:10:00',
#               'long': '2-00:00:00',
#               'xlong': '8-00:00:00',
#               'large': '2-00:00:00',
#               'xlarge': '2-00:00:00'}
# archivo = 'nuevo.inp' #str(sys.argv[1])
#
# fin = open(file_old, "rt")
# data = fin.read()
# data = data.replace('mypartition', partition)
# data = data.replace('mytime', time[partition])
# data = data.replace('myname', archivo)
# fin.close()
#
# fin = open(file_new, "wt")
# fin.write(data)
# fin.close()