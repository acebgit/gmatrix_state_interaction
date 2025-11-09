#########################################################
# PROGRAM TO ADD REQUIRED TEXT TO A
# FILE IN THE REQUIRED PORTION
#########################################################
import sys

inputs = sys.argv[1:]
search_line = '!PLOTS = 1'
lines_to_include = [
' ',
'!**SR-DFT**',
'RAS_SRDFT           TRUE !Perform short-range density functional RAS-CI',
'RAS_SRDFT_EXC       0 !EX functional. 0=HF.',
'RAS_SRDFT_ONLY_COR  TRUE !No use of EXC functional',
'RAS_SRDFT_COR       srPBE !short-range correlation functional',
'ras_srdft_0         TRUE !true: no iteration (one-shot); false: iteration',
'RAS_OMEGA           400 !Coulomb range-separation parameter. def=400',
'!RAS_SRDFT_DAMP     0.5 !dumping factor. def=0.5'
]

for file in inputs:
	with open(file, "r") as f:  # Read all lines
		lines = f.readlines()
	with open(file, "w") as f:  # Rewrite all lines except partition_line
		for line in lines:
			if search_line in line:
				f.write(line)
				for i in range(0, len(lines_to_include)):
					f.write(lines_to_include[i] + "\n")
			else:
				f.write(line)
