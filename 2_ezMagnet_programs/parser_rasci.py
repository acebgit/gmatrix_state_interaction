#########################################################
# THIS PROGRAM GENERATES A "JSON" FILE                  #
# (WITH SOCS, ENERGIES, SPIN AND ORBITAL ANGULAR        #
# MOMENTUMS OF ALL THE STATES) FROM THE RAS-CI Q-CHEM   #
# OUTPUT.                                               #
#########################################################

__author__ = 'Antonio Cebreiro-Gallardo, Abel Carreras'

# from dicttoxml import dicttoxml
import sys

import re
import operator
import warnings
import numpy as np
import hashlib, json
import pymatgen
from numpy import linalg, sqrt


class ParserError(Exception):
    def __init__(self, parser_name, message, output):
        self.parser_name = parser_name
        self.message = message
        self.full_output = output

    def __str__(self):
        return 'Error found while parsing output using "{}" parser: {}'.format(self.parser_name, self.message)


class OutputError(Exception):
    def __init__(self, output, error_output):
        self.full_output = output
        self.error_lines = error_output + '\n'.join(output.split('\n')[-20:])

    def __str__(self):
        return 'Error in Q-Chem calculation:\n{}'.format(self.error_lines)


class StructureError(Exception):
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return 'Error in Structure:\n{}'.format(self._message)


class QchemInputWarning(UserWarning):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


class QchemInputError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


class Structure:
    """
    Structure object containing all the geometric data of the molecule
    """
    def __init__(self,
                 coordinates=None,
                 symbols=None,
                 atomic_numbers=None,
                 connectivity=None,
                 charge=0,
                 multiplicity=1,
                 name=None):
        """
        :param coordinates: List containing the cartesian coordinates of each atom in Angstrom
        :param symbols: Symbols of the atoms within the molecule
        :param atomic_numbers: Atomic numbers of the atoms within the molecule
        :param charge: charge of the molecule
        :param multiplicity: multiplicity of the molecule
        """

        self._coordinates = np.array(coordinates)
        self._atomic_numbers = atomic_numbers
        self._connectivity = connectivity
        self._symbols = symbols
        self._charge = charge
        self._multiplicity = multiplicity
        self._name = name
        self._atomic_masses = None
        self._number_of_atoms = None

        # check input data
        if symbols is not None and coordinates is not None:
            if len(coordinates) != len(symbols):
                raise StructureError('coordinates and symbols do not match')

        if atomic_numbers is not None:
            self._symbols = [atom_data[i][1] for i in atomic_numbers]

    def __str__(self):
        return self.get_xyz()

    def __hash__(self):
        digest = hashlib.md5(json.dumps((self.get_xyz(), self.alpha_electrons, self.beta_electrons),
                                        sort_keys=True).encode()).hexdigest()
        return int(digest, 16)

    def get_coordinates(self, fragment=None):
        """
        gets the cartesian coordinates

        :param fragment: list of atoms that are part of the fragment

        :return: coordinates list
        """
        if self._coordinates is None:
            return None

        if fragment is None:
            return np.array(self._coordinates).tolist()
        else:
            return np.array(self._coordinates)[fragment].tolist()


    def set_coordinates(self, coordinates):
        """
        sets the cartessian coordinates

        :param coordinates: cartesian coordinates matrix
        """

        self._coordinates = np.array(coordinates)
        self._number_of_atoms = None

    @property
    def name(self):
        """
        returns the name
        :return: structure name
        """
        return self._name

    @property
    def file_name(self):
        return self._file_name

    @file_name.setter
    def file_name(self, file_name):
        self._file_name = file_name

    @property
    def charge(self):
        """
        returns the charge
        :return: the charge
        """
        return self._charge

    @charge.setter
    def charge(self, charge):
        self._charge = charge

    @property
    def multiplicity(self):
        """
        returns the multiplicity

        :return: the multiplicity
        """
        return self._multiplicity

    @multiplicity.setter
    def multiplicity(self, multiplicity):
        self._multiplicity = multiplicity

    @property
    def number_of_electrons(self):
        """
        returns the total number of electrons

        :return: number of total electrons
        """
        return int(np.sum(self.get_atomic_numbers()) + self.charge)

    @property
    def alpha_electrons(self):
        """
        returns the alpha electrons

        :return: number of alpha electrons
        """
        alpha_unpaired = self.multiplicity // 2
        return self.number_of_electrons // 2 + alpha_unpaired

    @property
    def beta_electrons(self):
        """
        returns the number of beta electrons

        :return: number of beta electrons
        """
        return self.number_of_electrons - self.alpha_electrons

    def get_atomic_numbers(self):
        """
        get the atomic numbers of the atoms of the molecule

        :return: list with the atomic numbers
        """
        if self._atomic_numbers is None:
            self._atomic_numbers = [[data[1].upper() for data in atom_data].index(element.upper())
                                    for element in self.get_symbols()]
        return self._atomic_numbers

    def set_atomic_numbers(self, atomic_numbers):
        self._atomic_numbers = atomic_numbers

    def get_symbols(self):
        """
        get the  atomic element symbols of the atoms of the molecule

        :return: list of symbols
        """
        if self._symbols is None:
            self._symbols = np.array(atom_data)[self.get_atomic_numbers()].T[1]
        return np.array([i for i in self._symbols if i != "X"], dtype=str)

    def set_symbols(self, atomic_elements):
        self._symbols = atomic_elements

    def _get_connectivity(self):
        if self._connectivity is None:
            print('No atom connectivity available')
            exit()

        return self._connectivity

    def _set_connectivity(self, connectivity):
        self._connectivity = connectivity

#   Real methods
    def get_number_of_atoms(self):
        """
        get the number of atoms

        :return: number of atoms
        """
        if self._number_of_atoms is None:
            self._number_of_atoms = np.array(self.get_coordinates()).shape[0]

        return self._number_of_atoms

    def get_atomic_masses(self):
        """
        get the atomic masses of the atoms of the molecule

        :return: list of atomic masses
        """
        if self._atomic_masses is None:

            try:
                masses_string = np.array(atom_data)[:, 3:4][[np.where(np.array(atom_data)==element)[0][0]
                                                             for element in self.get_symbols()]]
                self._atomic_masses = np.array(masses_string, dtype=float).T[0]
            except TypeError:
                print('Error reading element labels')
                exit()
        return self._atomic_masses

    def get_valence_electrons(self):
        """
        get number of valence electrons

        :return: number of valence electrons
        """
        valence_electrons = 0
        for number in self.get_atomic_numbers():
            if 2 >= number > 0:
                valence_electrons += np.mod(number, 2)
            if 18 >= number > 2:
                valence_electrons += np.mod(number-2, 8)
            if 54 >= number > 18:
                valence_electrons += np.mod(number-18, 18)
            if 118 >= number > 54:
                valence_electrons += np.mod(number-54, 32)
            if number > 118:
                raise Exception('Atomic number size not implemented')

        valence_electrons -= self.charge

        return valence_electrons

    def get_xyz(self, title=''):
        """
        generates a XYZ formatted file

        :param title: title of the molecule
        :return: string with the formatted XYZ file
        """
        txt = '{}\n{}'.format(self.get_number_of_atoms(), title)
        for s, c in zip(self.get_symbols(), self.get_coordinates()):
            txt += '\n{:2} '.format(s) + '{:15.10f} {:15.10f} {:15.10f}'.format(*c)

        return txt

    def get_connectivity(self, thresh=1.2):
        """
        get the connectivity as a list of pairs of indices of atoms
        from atomic radii

        :param thresh: radii threshold used to determine the connectivity
        :return:
        """
        from scipy.spatial import distance_matrix

        try:
            radius = [atom_data[sym][4] for sym in self.get_atomic_numbers()]
        except KeyError:
            warnings.warn('failed to generate connectivity, no connectivity will be used')
            return None

        distances_matrix = distance_matrix(self.get_coordinates(), self.get_coordinates())

        radii_matrix = np.array([radius] * len(radius))
        radii_matrix = radii_matrix + radii_matrix.T

        try:
            relative_differences = np.abs(radii_matrix - distances_matrix) / radii_matrix
        except ValueError:
            warnings.warn('failed to generate connectivity')
            return None

        if not (np.array(np.where(relative_differences < thresh - 1)).T + 1).tolist():
            return None
        else:
            return (np.array(np.where(relative_differences < thresh - 1)).T + 1).tolist()

    def get_point_symmetry(self):
        """
        Returns the point group of the molecule using pymatgen

        :return: point symmetry label
        """
        from pymatgen.core import Molecule
        from pymatgen.symmetry.analyzer import PointGroupAnalyzer

        pymatgen_mol = Molecule(self.get_symbols(), self.get_coordinates())
        symm_group = PointGroupAnalyzer(pymatgen_mol, tolerance=0.1)

        return symm_group.sch_symbol


atom_data = [
    # atomic number, symbols, names, masses, bohr radius
    [  0, "X", "X",            0.000000, 0.000],  # 0
    [  1, "H", "Hydrogen",     1.007940, 0.324],  # 1
    [  2, "He", "Helium",      4.002602, 0.000],  # 2
    [  3, "Li", "Lithium",     6.941000, 1.271],  # 3
    [  4, "Be", "Beryllium",   9.012182, 0.927],  # 4
    [  5, "B", "Boron",       10.811000, 0.874],  # 5
    [  6, "C", "Carbon",      12.010700, 0.759],  # 6
    [  7, "N", "Nitrogen",    14.006700, 0.706],  # 7
    [  8, "O", "Oxygen",      15.999400, 0.678],  # 8
    [  9, "F", "Fluorine",    18.998403, 0.568],  # 9
    [ 10, "Ne", "Neon",       20.179700, 0.000],  # 10
    [ 11, "Na", "Sodium",     22.989769, 1.672],  # 11
    [ 12, "Mg", "Magnesium",  24.305000, 1.358],  # 12
    [ 13, "Al", "Aluminium",  26.981539, 1.218],  # 13
    [ 14, "Si", "Silicon",    28.085500, 1.187],  # 14
    [ 15, "P", "Phosphorus",  30.973762, 1.105],  # 15
    [ 16, "S", "Sulfur",      32.065000, 1.045],  # 16
    [ 17, "Cl", "Chlorine",   35.453000, 1.006],  # 17
    [ 18, "Ar", "Argon",      39.948000, 0.000],  # 18
    [ 19, "K", "Potassium",   39.098300, 2.247],  # 19
    [ 20, "Ca", "Calcium",    40.078000, 1.748],  # 20
    [ 21, "Sc", "Scandium",   44.955912, 1.664],  # 21
    [ 22, "Ti", "Titanium",   47.867000, 1.620],  # 22
    [ 23, "V", "Vanadium",    50.941500, 1.543],  # 23
    [ 24, "Cr", "Chromium",   51.996100, 1.418],  # 24
    [ 25, "Mn", "Manganese",  54.938045, 1.569],  # 25
    [ 26, "Fe", "Iron",       55.845000, 1.514],  # 26
    [ 27, "Co", "Cobalt",     58.933195, 1.385],  # 27
    [ 28, "Ni", "Nickel",     58.693400, 1.390],  # 28
    [ 29, "Cu", "Copper",     63.546000, 1.382],  # 29
    [ 30, "Zn", "Zinc",       65.380000, 1.416],  # 30
    [ 31, "Ga", "Gallium",    69.723000, 1.235],  # 31
    [ 32, "Ge", "Germanium",  72.640000, 1.201],  # 32
    [ 33, "As", "Arsenic",    74.921600, 1.232],  # 33
    [ 34, "Se", "Selenium",   78.960000, 1.210],  # 34
    [ 35, "Br", "Bromine",    79.904000, 1.190],  # 35
    [ 36, "Kr", "Krypton",    83.798000, 0.000],  # 36
    [ 37, "Rb", "Rubidium",   85.467800, 2.284],  # 37
    [ 38, "Sr", "Strontium",  87.620000, 1.942],  # 38
    [ 39, "Y", "Yttrium",     88.905850, 1.993],  # 39
    [ 40, "Zr", "Zirconium",  91.224000, 1.758],  # 40
    [ 41, "Nb", "Niobium",    92.906380, 1.610],  # 41
    [ 42, "Mo", "Molybdenum", 95.960000, 1.639],  # 42
    [ 43, "Tc", "Technetium",  0.000000, 1.493],  # 43
    [ 44, "Ru", "Ruthenium",  101.07000, 1.467],  # 44
    [ 45, "Rh", "Rhodium",    102.90550, 1.437],  # 45
    [ 46, "Pd", "Palladium",  106.42000, 1.422],  # 46
    [ 47, "Ag", "Silver",     107.86820, 1.466],  # 47
    [ 48, "Cd", "Cadmium",    112.41100, 1.441],  # 48
    [ 49, "In", "Indium",     114.81800, 1.421],  # 49
    [ 50, "Sn", "Tin",        118.71000, 1.408],  # 50
    [ 51, "Sb", "Antimony",   121.76000, 1.397],  # 51
    [ 52, "Te", "Tellurium",  127.60000, 1.395],  # 52
    [ 53, "I", "Iodine",      126.90447, 1.396],  # 53
    [ 54, "Xe", "Xenon",      131.29300, 1.336],  # 54
    [ 55, "Cs", "Caesium",    132.90545, 2.470],  # 55
    [ 56, "Ba", "Barium",     137.32700, 2.219],  # 56
    [ 57, "La", "Lanthanum",  138.90547, 2.089],  # 57
    [ 58, "Ce", "Cerium",     140.11600, 2.054],  # 58
    [ 59, "Pr", "Praseodymium", 140.90765, 1.979],  # 59
    [ 60, "Nd", "Neodymium",  144.24200, 0.000],  # 60
    [ 61, "Pm", "Promethium",   0.00000, 0.000],  # 61
    [ 62, "Sm", "Samarium",   150.36000, 2.535],  # 62
    [ 63, "Eu", "Europium",   151.96400, 0.000],  # 63
    [ 64, "Gd", "Gadolinium", 157.25000, 0.000],  # 64
    [ 65, "Tb", "Terbium",    158.92535, 0.000],  # 65
    [ 66, "Dy", "Dysprosium", 162.50000, 0.000],  # 66
    [ 67, "Ho", "Holmium",    164.93032, 0.000],  # 67
    [ 68, "Er", "Erbium",     167.25900, 0.000],  # 68
    [ 69, "Tm", "Thulium",    168.93421, 0.000],  # 69
    [ 70, "Yb", "Ytterbium",  173.05400, 0.000],  # 70
    [ 71, "Lu", "Lutetium",   174.96680, 0.000],  # 71
    [ 72, "Hf", "Hafnium",    178.49000, 1.779],  # 72
    [ 73, "Ta", "Tantalum",   180.94788, 1.723],  # 73
    [ 74, "W", "Tungsten",    183.84000, 1.627],  # 74
    [ 75, "Re", "Rhenium",    186.20700, 1.536],  # 75
    [ 76, "Os", "Osmium",     190.23000, 1.521],  # 76
    [ 77, "Ir", "Iridium",    192.21700, 1.456],  # 77
    [ 78, "Pt", "Platinum",   195.08400, 1.390],  # 78
    [ 79, "Au", "Gold",       196.96657, 1.402],  # 79
    [ 80, "Hg", "Mercury",    200.59000, 1.371],  # 80
    [ 81, "Tl", "Thallium",   204.38330, 1.384],  # 81
    [ 82, "Pb", "Lead",       207.20000, 1.820],  # 82
    [ 83, "Bi", "Bismuth",    208.98040, 1.507],  # 83
    [ 84, "Po", "Polonium",     0.00000, 0.000],  # 84
    [ 85, "At", "Astatine",     0.00000, 0.000],  # 85
    [ 86, "Rn", "Radon",        0.00000, 0.000],  # 86
    [ 87, "Fr", "Francium",     0.00000, 0.000],  # 87
    [ 88, "Ra", "Radium",       0.00000, 0.000],  # 88
    [ 89, "Ac", "Actinium",     0.00000, 0.000],  # 89
    [ 90, "Th", "Thorium",    232.03806, 0.000],  # 90
    [ 91, "Pa", "Protactinium",231.03588, 0.000], # 91
    [ 92, "U", "Uranium",     238.02891, 0.000],  # 92
    [ 93, "Np", "Neptunium",    0.00000, 0.000],  # 93
    [ 94, "Pu", "Plutonium",    0.00000, 0.000],  # 94
    [ 95, "Am", "Americium",    0.00000, 0.000],  # 95
    [ 96, "Cm", "Curium",       0.00000, 0.000],  # 96
    [ 97, "Bk", "Berkelium",    0.00000, 0.000],  # 97
    [ 98, "Cf", "Californium",  0.00000, 0.000],  # 98
    [ 99, "Es", "Einsteinium",  0.00000, 0.000],  # 99
    [100, "Fm", "Fermium",      0.00000, 0.000],  # 100
    [101, "Md", "Mendelevium",  0.00000, 0.000],  # 101
    [102, "No", "Nobelium",     0.00000, 0.000],  # 102
    [103, "Lr", "Lawrencium",   0.00000, 0.000],  # 103
    [104, "Rf", "Rutherfordium",0.00000, 0.000],  # 104
    [105, "Db", "Dubnium",      0.00000, 0.000],  # 105
    [106, "Sg", "Seaborgium",   0.00000, 0.000],  # 106
    [107, "Bh", "Bohrium",      0.00000, 0.000],  # 107
    [108, "Hs", "Hassium",      0.00000, 0.000],  # 108
    [109, "Mt", "Meitnerium",   0.00000, 0.000],  # 109
    [110, "Ds", "Darmstadtium", 0.00000, 0.000],  # 110
    [111, "Rg", "Roentgenium",  0.00000, 0.000],  # 111
    [112, "Cn", "Copernicium",  0.00000, 0.000],  # 112
    [113, "Uut", "Ununtrium",   0.00000, 0.000],  # 113
    [114, "Uuq", "Ununquadium", 0.00000, 0.000],  # 114
    [115, "Uup", "Ununpentium", 0.00000, 0.000],  # 115
    [116, "Uuh", "Ununhexium",  0.00000, 0.000],  # 116
    [117, "Uus", "Ununseptium", 0.00000, 0.000],  # 117
    [118, "Uuo", "Ununoctium",  0.00000, 0.000],  # 118
    ]


def read_basic_info(output):

    there_vector = [m.start() for m in re.finditer('There are ', output)]
    n_alpha = int(output[there_vector[0]:there_vector[0]+100].split()[2])
    n_beta = int(output[there_vector[0]:there_vector[0]+100].split()[5])

    nshell = int(output[there_vector[1]:there_vector[1]+100].split()[2])
    nbas = int(output[there_vector[1]:there_vector[1]+100].split()[5])

    return {'n_alpha': n_alpha,
            'n_beta': n_beta,
            'n_shells': nshell,
            'n_basis_functions': nbas}


def get_rasci_occupations_list(configuration, occupied_orbitals, total_orbitals):
    # occupied_orbitals = get_occupied_electrons(configuration, structure)
    n_extra = total_orbitals - occupied_orbitals - len(configuration['alpha'])
    vector_alpha = [1] * occupied_orbitals + [int(c) for c in configuration['alpha']] + [0] * n_extra

    n_extra = total_orbitals - occupied_orbitals - len(configuration['beta'])
    vector_beta = [1] * occupied_orbitals + [int(c) for c in configuration['beta']] + [0] * n_extra

    if configuration['hole'] != '':
        if np.sum(vector_alpha) > np.sum(vector_beta):
            vector_alpha[int(configuration['hole']) - 1] = 0
        else:
            vector_beta[int(configuration['hole']) - 1] = 0

    if configuration['part'] != '':
        if np.sum(vector_alpha) < np.sum(vector_beta):
            vector_alpha[int(configuration['part']) - 1] = 1
        else:
            vector_beta[int(configuration['part']) - 1] = 1

    return {'alpha': vector_alpha, 'beta': vector_beta}


def search_bars(output, from_position=0, bar_type='---'):
    output = output[from_position:]
    positions = []
    previous = 0
    for m in re.finditer(bar_type, output):
        if m.start() > previous + 1:
            positions.append(m.start() + from_position)
        previous = m.end()

    return positions


def standardize_vector(vector):
    import numpy as np
    if vector[0] != 0:
        if vector[0] < 0:
            vector = np.array(vector) * -1
            vector = vector.tolist()
    elif vector[1] != 0:
        if vector[1] < 0:
            vector = np.array(vector) * -1
            vector = vector.tolist()
    else:
        if vector[2] < 0:
            vector = np.array(vector) * -1
            vector = vector.tolist()

    for i in range(3):
        vector[i] = vector[i] + 0

    return vector


def _read_soc_matrix(lines, dimensions):
    # for line in lines:
    #     print(line)

    col_per_line = 5
    matrix = []
    for ib in range(dimensions[0]):
        real = []
        complex = []
        for j in range((dimensions[1] - 1) // col_per_line + 1):
            real += lines[j*dimensions[0] + 1 * (j+1) + ib][11:].split()[0::2]
            complex += lines[j*dimensions[0] + 1 * (j+1) +ib][11:].split()[1::2]

        row = [float(r) + float(c[:-1]) * 1j for r, c in zip(real, complex)]
        matrix.append(row)

    return matrix


def read_input_structure(output):

    enum = output.find('Standard Nuclear Orientation')
    end_section = search_bars(output, from_position=enum, bar_type=r'-----')
    section_structure = output[end_section[0]: end_section[1]].split('\n')

    symbols = []
    coordinates = []
    for line in section_structure[1:-1]:
        symbols.append(line.split()[1])
        coordinates.append(line.split()[2:5])

    coordinates = np.array(coordinates, dtype=float).tolist()

    # basic info
    enum = output.find('Nuclear Repulsion Energy')
    basic_data = read_basic_info(output[enum:enum + 5000])

    n_nucleus = 0
    for s in symbols:
        for i, row in enumerate(atom_data):
            if row[1] == s:
                n_nucleus += i

    charge = n_nucleus - (basic_data['n_alpha'] + basic_data['n_beta'])
    multiplicity = abs(basic_data['n_alpha'] - basic_data['n_beta']) + 1

    return Structure(coordinates=coordinates,
                     symbols=symbols,
                     charge=charge,
                     multiplicity=multiplicity)


def _complete_interstate_pairs(interstate_dict):
    additional_items = {}
    for key, value in interstate_dict.items():
        if not key[::-1] in interstate_dict:
            dict_entry = {}
            for k, v in interstate_dict[key].items():
                if k == 'state_a':
                    dict_entry.update({k: key[1]})
                elif k == 'state_b':
                    dict_entry.update({k: key[0]})
                elif k == 'angular_momentum':
                    dict_entry.update({k: np.conjugate(v).tolist()})
                elif k == '1e_soc_mat':
                    dict_entry.update({k: np.conjugate(v).T.tolist()})
                elif k == 'hso_l-':
                    dict_entry.update({k: (-np.array(v)).tolist()})
                elif k == 'hso_l+':
                    dict_entry.update({k: (-np.array(v)).tolist()})
                elif k == '2e_soc_mat':
                    dict_entry.update({k: np.conjugate(v).T.tolist()})
                elif k == 'total_soc_mat':
                    dict_entry.update({k: np.conjugate(v).T.tolist()})
                else:
                    dict_entry.update({k: v})
            additional_items.update({key[::-1]: dict_entry})

    interstate_dict.update(additional_items)


def parser_rasci(output):
    """
    Parser for RAS-CI calculations
    Include:
    - Diabatization scheme data
    - Structure
    - Adiabatic states
    - SOC

    :param output:
    :return:
    """

    data_dict = {}

    # Molecule
    data_dict['structure'] = read_input_structure(output)
    n_atoms = data_dict['structure'].get_number_of_atoms()

    # basic info
    enum = output.find('Nuclear Repulsion Energy')
    basic_data = read_basic_info(output[enum:enum + 5000])

    # scf_energy
    enum = output.find('SCF   energy in the final basis set')
    scf_energy = float(output[enum:enum+100].split()[8])

    data_dict['scf_energy'] = scf_energy
    # total energy
    # enum = output.find('Total energy in the final basis set')
    # total_energy = float(output[enum:enum+100].split()[8])

    # RASCI dimensions
    ini_section = output.find('RAS-CI Dimensions')
    end_section = search_bars(output, from_position=ini_section, bar_type=r'\*\*\*')[0]
    dimension_section = output[ini_section: end_section]

    enum = dimension_section.find('Doubly Occ')
    doubly_occ = int(dimension_section[enum: enum+50].split()[3])
    enum = dimension_section.find('Doubly Vir')
    doubly_vir = int(dimension_section[enum: enum+50].split()[2])
    enum = dimension_section.find('Frozen Occ')
    frozen_occ = int(dimension_section[enum: enum+50].split()[3])
    enum = dimension_section.find('Frozen Vir')
    frozen_vir = int(dimension_section[enum: enum+50].split()[2])

    enum = dimension_section.find('Total CI configurations')
    total_conf = int(dimension_section[enum: enum+50].split()[3])
    enum = dimension_section.find('Active configurations')
    active_conf = int(dimension_section[enum: enum+50].split()[2])
    enum = dimension_section.find('Hole configurations')
    hole_conf = int(dimension_section[enum: enum+50].split()[2])
    enum = dimension_section.find('Particle configurations')
    particle_conf = int(dimension_section[enum: enum+50].split()[2])

    rasci_dimensions = {'doubly_occupied': doubly_occ,
                        'doubly_virtual': doubly_vir,
                        'frozen_occupied': frozen_occ,
                        'frozen_virtual': frozen_vir,
                        'total_configurations': total_conf,
                        'active_configurations': active_conf,
                        'hole_configurations': hole_conf,
                        'particle_configurations': particle_conf}

    data_dict.update({'rasci_dimensions': rasci_dimensions})

    # excited states data
    excited_states = []
    for m in re.finditer('RAS-CI total energy for state', output):
        # print('ll found', m.start(), m.end())

        end_section = search_bars(output, from_position=m.start(), bar_type=r'\*\*\*\*\*\*\*')[0]
        section_state = output[m.start():end_section]

        # energies
        enum = section_state.find('total energy for state')
        tot_energy = float(section_state[enum: enum + 50].split()[5])
        enum = section_state.find('Excitation energy')
        exc_energy_units = section_state[enum: enum + 30].split()[2].strip('(').strip(')')
        exc_energy = float(section_state[enum: enum + 50].split()[4])

        # multiplicity
        n_multi = section_state.find('<S^2>')
        multi_data = section_state[n_multi:n_multi + 30].split(':')[1]
        state_multiplicity = float(multi_data.split()[0])

        # dipole moment
        enum = section_state.find('Dipole Moment')
        dipole_mom = [float(section_state[enum:].split()[2]) + 0.0,
                      float(section_state[enum:].split()[4]) + 0.0,
                      float(section_state[enum:].split()[6]) + 0.0]

        # Transition moment
        enum = section_state.find('Trans. Moment')
        trans_mom = strength = None
        if enum > -1:
            trans_mom = [float(section_state[enum:].split()[2]) + 0.0,
                         float(section_state[enum:].split()[4]) + 0.0,
                         float(section_state[enum:].split()[6]) + 0.0]
            trans_mom = standardize_vector(trans_mom)
            strength = float(section_state[enum:].split()[10])

        # Mulliken population analysis
        mulliken_population = None
        enum = section_state.find('Mulliken population analysis')
        if enum > -1:
            max = search_bars(section_state, from_position=enum, bar_type=r'\-'*30)
            pop_charges = []
            pop_spin = []
            for line in section_state[max[1]: max[2]].split('\n')[1:n_atoms+1]:
                c, s = line.split()[2:4]
                pop_charges.append(float(c))
                pop_spin.append(float(s))

            mulliken_population = {'charge': pop_charges, 'spin': pop_spin, 'units': 'au'}

        # Natural orbitals
        nato_occ = None
        enum = section_state.find('NATURAL OCCUPATION NUMBERS')
        if enum > -1:
            lines = []
            for line in section_state[enum:].split('\n')[2::2]:
                if len(line) == 0:
                    break
                if line.split()[0].isnumeric():
                    lines += line.split()[1:]
                else:
                    break
            nato_occ = [float(num) for num in lines]

        # configurations table
        enum = section_state.find('AMPLITUDE')
        enum2 = section_state.find('Contributions')
        section_table = section_state[enum: enum2].split('\n')[2:-2]

        # ' HOLE  | ALPHA | BETA  | PART | AMPLITUDE'
        table = []
        for row in section_table:
            table.append({'hole': row.split('|')[1].strip(),
                          'alpha': row.split('|')[2].strip(),
                          'beta': row.split('|')[3].strip(),
                          'part': row.split('|')[4].strip(),
                          'amplitude': float(row.split('|')[5]) + 0.0})
            table[-1]['occupations'] = get_rasci_occupations_list(table[-1],
                                                                  doubly_occ,
                                                                  basic_data['n_basis_functions'])

        table = sorted(table, key=operator.itemgetter('hole', 'alpha', 'beta', 'part'))
        table = sorted(table, key=lambda x: abs(x['amplitude']), reverse=True)

        # Contributions RASCI wfn
        contributions_section = section_state[enum2:]
        contributions = {'active' : float(contributions_section.split()[4]),
                         'hole': float(contributions_section.split()[6]),
                         'part': float(contributions_section.split()[8])}

        # complete dictionary
        tot_energy_units = 'au'

        state_dict = {'total_energy': tot_energy,
                      'total_energy_units': tot_energy_units,
                      'excitation_energy': exc_energy,
                      'excitation_energy_units': exc_energy_units,
                      'multiplicity': state_multiplicity,
                      'dipole_moment': dipole_mom,
                      'transition_moment': trans_mom,
                      'dipole_moment_units': 'ua',
                      'oscillator_strength': strength,
                      'configurations': table,
                      'contributions_wfn': contributions}

        if nato_occ is not None:
            state_dict.update({'natural_occupation_numbers': nato_occ})

        if mulliken_population is not None:
            state_dict.update({'mulliken_population': mulliken_population})

        excited_states.append(state_dict)

    data_dict.update({'excited_states': excited_states})

    # Interstate transition properties
    done_interstate = bool(output.find('Interstate Transition Properties')+1)
    if done_interstate:
        ini_section = output.find('Interstate Transition Properties')
        end_section = search_bars(output, from_position=ini_section)[1]
        interstate_section = output[ini_section: end_section]

        interstate_dict = {}
        for m in re.finditer('State A: Root', interstate_section):
            section_pair = interstate_section[m.start():m.start() + 10000]
            end_section = search_bars(section_pair, bar_type=r'\*'*20)[0]
            section_pair = section_pair[:end_section]

            lines = section_pair.split('\n')

            state_a = int(lines[0].split()[-1])
            state_b = int(lines[1].split()[-1])

            pair_dict = {'state_a': state_a,
                         'state_b': state_b}

            s_a = s_b = 0
            for i, line in enumerate(lines):
                # RAS-CI SOC section
                if 'Angular momentum components' in line:
                    pair_dict['angular_momentum'] = [complex(lines[i+1+k].split()[-1].replace('i', 'j')) for k in range(3)]

                if '||gamma^AB||_total' in line:
                    pair_dict['gamma_total'] = float(lines[i+0].split()[-1])
                    pair_dict['gamma_sym'] = float(lines[i+1].split()[-1])
                    pair_dict['gamma_anti_sym'] = float(lines[i+2].split()[-1])

                if "KET: S',Sz'" in line:
                    s_a = float(lines[i].split('=')[1].split()[0])
                    s_b = float(lines[i+1].split('=')[1].split()[0])

                na = int(2 * s_a + 1)
                nb = int(2 * s_b + 1)

                if 'Spin Matrices' in line:
                    warnings.warn('Spin Matrices parsing is deprecated')
                    spinmat_x = _read_soc_matrix(lines[i + 2:], [nb, na])
                    spinmat_y = _read_soc_matrix(lines[i + 4 + nb:], [nb, na])
                    spinmat_z = _read_soc_matrix(lines[i + 6 + 2*nb:], [nb, na])
                    pair_dict['spin_matrices'] = [spinmat_x, spinmat_y, spinmat_z]

                if 'Spin matrices Sx, Sy and Sz for states' in line:
                    warnings.warn('Spin Matrices parsing is deprecated')
                    pair_dict['spin_matrices'] = [np.zeros((nb, na)).tolist()]*3

                if '1-elec SOC matrix (cm-1)' in line:
                    pair_dict['1e_soc_mat'] = _read_soc_matrix(lines[i+1:], [nb, na])

                if '1-elec SOCC' in line:
                    pair_dict['1e_socc'] = float(line.split()[3])

                if '2e-SOMF Reduced matrix elements (cm-1)' in line:
                    r, c = lines[i+1].split()[-2:]
                    pair_dict['hso_l-'] = float(r) + float(c) * 1j
                    r, c = lines[i+2].split()[-2:]
                    pair_dict['hso_l0'] = float(r) + float(c) * 1j
                    r, c = lines[i+3].split()[-2:]
                    pair_dict['hso_l+'] = float(r) + float(c) * 1j

                if '2-elec mean-field SOC matrix (cm-1)' in line:
                    pair_dict['2e_soc_mat'] = _read_soc_matrix(lines[i + 1:], [nb, na])
                if 'Total mean-field SOC matrix (cm-1)' in line:
                    pair_dict['total_soc_mat'] = _read_soc_matrix(lines[i + 1:], [nb, na])
                if 'Mean-Field SOCC' in line:
                    pair_dict['mf_socc'] = float(line.split()[-2])
                    pair_dict['units'] = line.split()[-1]

                if 'Skipping SOCs between states' in line:
                    pair_dict['1e_soc_mat'] = np.zeros((nb, na)).tolist()
                    pair_dict['1e_socc'] = 0.0

                    pair_dict['hso_l-'] = complex(0.0)
                    pair_dict['hso_l0'] = complex(0.0)
                    pair_dict['hso_l+'] = complex(0.0)

                    pair_dict['2e_soc_mat'] = np.zeros((nb, na)).tolist()
                    pair_dict['total_soc_mat'] = np.zeros((nb, na)).tolist()

                    pair_dict['mf_socc'] = 0.0
                    pair_dict['units'] = 'cm-1'

            interstate_dict[(state_a, state_b)] = pair_dict

        _complete_interstate_pairs(interstate_dict)
        data_dict.update({'interstate_properties': interstate_dict})

    return data_dict


def get_symmetry_states(file, totalstates):
    """
    Create two lists: first with the symmetry of each state (A1,A2,A3...) and second with
    the order of these symmetries (1A1,2A2,3A3...)
    :param: file_ras, nstates
    :return: all_state_symmetries, ordered_state_symmetries
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    word_search = ['State symmetry: ']
    elements = []
    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[2]
            element = element.replace('*', '')  # '*' means mix of symmetries
            elements.append(element)

        if len(elements) == totalstates:
            break  # All symmetries appended
    all_state_symmetries = np.array(elements)

    ordered_state_symmetries = []  # State symmetries with the order (1,2,3..)
    for i in range(0, len(all_state_symmetries)):
        number = 1
        for j in range(0, i):
            if all_state_symmetries[i] == all_state_symmetries[j]:
                number += 1
        element = str(number) + all_state_symmetries[i]
        ordered_state_symmetries.append(element)
    return all_state_symmetries, ordered_state_symmetries


def get_number_of_states(file):
    """
    Obtain the total number of states in ras
    :param: file
    :return: nstates
    """
    with open(file, encoding="utf8") as f:
        data = f.readlines()

    element = 0
    word_search = ['Requested states: ']

    for line in data:
        if any(i in line for i in word_search):
            line = line.split()
            element = line[3]
            break
    totalstates = int(element)
    states = list(range(1, totalstates + 1))
    return totalstates, states


class MyEncoder(json.JSONEncoder):
    # https://itsourcecode.com/typeerror/typeerror-object-of-type-is-not-json-serializable-solved/?expand_article=1
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, complex):
            # return {'real': obj.real, 'image': obj.imag}
            a = obj.real
            b = obj.imag
            return {a, b}
        else:
            return super().default(obj)


def get_socs(outpuut, nstates, selected_state):
    """
    Obtain a dictionary with SOCs between ground state and all the rest, written as strings.
    :param outpuut:
    :param nstates:
    :param selected_state:
    :return: socs
    """
    data = outpuut['interstate_properties']
    socs = {}
    for i in range(2, nstates + 1):
        list_1 = data[(1, i)]['total_soc_mat']
        a = str(list_1[0][0]), str(list_1[0][1])
        b = str(list_1[1][0]), str(list_1[1][1])
        list_2 = [a, b]
        socs.update({selected_state[i - 1]: list_2})
    return socs


def get_energies(outpuut, nstates, selected_state):
    """
    Obtain a dictionary with eigenenergies of all the states, written as floats.
    :param outpuut:
    :param nstates:
    :param selected_state:
    :return: socs
    """
    data = outpuut['excited_states']
    energy = {}
    for i in range(0, nstates):
        lista = [data[i]['total_energy'], data[i]['excitation_energy']]
        energy.update({selected_state[i]: lista})
    return energy


def get_orbital_angmoment(outpuut, nstates, selected_state):
    """
    Obtain a dictionary with orbital angular momentum between ground state and all the rest, written as strings.
    :param outpuut:
    :param nstates:
    :param selected_state:
    :return: socs
    """
    data = outpuut['interstate_properties']
    momentums = {}
    for i in range(2, nstates + 1):
        list_1 = data[(1, i)]['angular_momentum']
        a = str(list_1[0])
        b = str(list_1[1])
        c = str(list_1[2])
        list_2 = [a, b, c]
        momentums.update({selected_state[i - 1]: list_2})
    return momentums


def get_spin_momentum(outpuut, nstates):
    """
    get s2 of each state from Q-Chem otuput
    :param: outpuut
    :param: nstates
    :return: s2
    """
    search = ['  <S^2>      : ']
    s2_each_states = []
    with open(outpuut, encoding="utf8") as outpuut:
        for line in outpuut:
            if any(ii in line for ii in search):
                element = float(line[16:])
                s2_each_states.append(element)
    # s2_each_states = np.array(elements, dtype=float)

    momentum = {}
    for i in range(0, totalstates):
        momentum.update({nstates[i]: s2_each_states[i]})
    return momentum


def output_json(outpuut):
    outfile_name = outpuut + ".json"
    with open(outfile_name, 'w') as ff:
        json.dump(output_dict, ff, separators=(',', ':'), sort_keys=True, indent=4, cls=MyEncoder)


# def output_xml(outpuut, dictionary):
#     xmloutfile_name = outpuut + ".xml"
#     xml = dicttoxml(dictionary)
#     with open(xmloutfile_name, 'w') as ff:
#         ff.write(str(xml))


file = str(sys.argv[1])
with open(file, encoding="utf8") as f:
    output = f.read()
output = parser_rasci(output)

totalstates, states_ras = get_number_of_states(file)
all_state_symmetry, ordered_state_symmetry = get_symmetry_states(file, totalstates)
ordered_pair_symmetry = [ordered_state_symmetry[0] + "_" + stateSymm for stateSymm in ordered_state_symmetry]
aveTransSOCListDict = get_socs(output, totalstates, ordered_pair_symmetry)
stateEnergiesDict = get_energies(output, totalstates, ordered_state_symmetry)
stateTotalAngMomDict = get_spin_momentum(file, ordered_state_symmetry)
transAngMomListDict = get_orbital_angmoment(output, totalstates, ordered_pair_symmetry)

output_dict = {
    "stateEnergiesDict" : stateEnergiesDict,
    "stateTotalAngMomDict" : stateTotalAngMomDict,
    "transAngMomListDict" : transAngMomListDict,
    "aveTransSOCListDict" : aveTransSOCListDict
}

output_json(file)
# output_xml(file, output_dict)
