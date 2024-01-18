import parser_qchem_tddft, parser_qchem_rasci, parser_qchem_eom
import subprocess

state_selection = 1  # 0: use "state_ras" ; 1: use all states_selected
initial_states = [1,2]

file = '../test/qchem_eomcc.out'  # str(sys.argv[1])


def check_parser(filee):
    with open(filee, encoding="utf8") as f:
        for line in f:
            if 'EOM-MP2' in line:
                parser = parser_qchem_eom
            if 'CORRELATION        RASCI' in line:
                parser = parser_qchem_rasci
            if 'TDDFT/TDA Excitation Energies' in line:
                parser = parser_qchem_tddft
    return parser

parser_selected = check_parser(file)
subprocess.run(parser_selected)