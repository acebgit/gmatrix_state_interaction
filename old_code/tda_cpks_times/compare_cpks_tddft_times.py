import os
import re
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from adjustText import adjust_text

folder_path = sys.argv[1]  # PUT THE FOLDER NAME
saveplot = 1
tda_fit_degree = 3
cpks_fit_degree = 2


def extract_cpks_time_basisfunct(file_path):
    """Extracts time in seconds from a CPKS output file."""
    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(r"TOTAL RUN TIME: (\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec", line)
            if match:
                days, hours, minutes, seconds, msec = map(int, match.groups())
                total_seconds = (days * 86400) + (hours * 3600) + (minutes * 60) + seconds + (msec / 1000)
            
            match2 = re.search(r"Number of basis functions", line)
            if match2:
                basis_functions = int(line.split()[-1])            
    return total_seconds, basis_functions


def extract_tda_time(file_path):
    """Extracts CPU time in seconds from a TDA output file."""
    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(r"Total job time:.*?([\d.]+)s\(cpu\)", line)
            if match:
                return float(match.group(1))
    return None


def take_times_and_basis(folder_path):
    """Processes all CPKS and TDA files in the folder and returns a dictionary with times."""
    times = {}

    for filename in os.listdir(folder_path):
        if filename.endswith("cpks.out"):
            molecule = filename.replace("_cpks.out", "")
            file_path = os.path.join(folder_path, filename)
            cpks_time, basis_functions = extract_cpks_time_basisfunct(file_path)
            if cpks_time is not None:
                if molecule in times:
                    times[molecule]["cpks_time"] = cpks_time
                else:
                    times[molecule] = {"cpks_time": cpks_time}
                times[molecule]["basis_functions"] = basis_functions

        elif filename.endswith("tda.out"):
            molecule = filename.replace("_tda.out", "")
            file_path = os.path.join(folder_path, filename)
            tda_time = extract_tda_time(file_path)
            if tda_time is not None:
                if molecule in times:
                    times[molecule]["tda_time"] = tda_time
                else:
                    times[molecule] = {"tda_time": tda_time}

    return times


def compute_ratio(data):
    result = {}
    for molecule, values in data.items():
        ratio = round(values['tda_time'] / values['cpks_time'], 3)
        result[molecule] = {'ratio': ratio, 'basis_functions': values['basis_functions']}
    return result


def split_molecule_data(data):
    cpks_dict = {}
    tda_dict = {}

    for molecule, values in data.items():
        cpks_dict[molecule] = {'CPKS': values['cpks_time'], 'basis_functions': values['basis_functions']}
        tda_dict[molecule] = {'TDA': values['tda_time'], 'basis_functions': values['basis_functions']}

    return cpks_dict, tda_dict


def fit_cpks_vs_basis(data, method, degree):
    # Extract basis functions and CPKS values
    basis_functions = np.array([values['basis_functions'] for values in data.values()])
    cpks_values = np.array([values[method] for values in data.values()])

    # Perform polynomial fitting
    coefficients = np.polyfit(basis_functions, cpks_values, degree)
    polynomial = np.poly1d(coefficients)

    # Calculate R² (coefficient of determination)
    y_pred = polynomial(basis_functions)
    ss_total = np.sum((cpks_values - np.mean(cpks_values)) ** 2)
    ss_residual = np.sum((cpks_values - y_pred) ** 2)
    r_squared = 1 - (ss_residual / ss_total)

    # Generate the equation as a string
    equation_terms = [f"{coeff:.6f} * x^{i}" if i > 0 else f"{coeff:.6f}" 
                      for i, coeff in enumerate(coefficients[::-1])]
    equation = " + ".join(equation_terms)

    return coefficients, equation, r_squared


def convert_data_to_minutes(data, method):
    """
    Converts the 'CPKS' values in the given dictionary from seconds to minutes.

    Parameters:
        data (dict): Dictionary containing molecule names as keys and 
                     their properties (including 'CPKS' in seconds) as values.

    Returns:
        dict: New dictionary with 'CPKS' values converted to minutes.
    """
    # Create a new dictionary with 'CPKS' values converted to minutes
    converted_data = {}
    for molecule, values in data.items():
        # Convert CPKS from seconds to minutes
        converted_data[molecule] = values.copy()
        converted_data[molecule][method] = values[method] / 60
    
    return converted_data


def plot_scatter(data, method, degree=1, save_plot=False, font_size=14, file_name='scatter_plot.png'):
    """
    Plots a scatter plot of Basis Functions vs. a selected method value and fits a polynomial regression line.

    Parameters:
        data (dict): Dictionary containing 'basis_functions' and method values.
        method (str): The key to extract y-values (e.g., 'CPKS', 'TDA/CPKS').
        degree (int): The degree of the polynomial regression line.
        save_plot (bool): Whether to save the plot as a file.
        font_size (int): Font size for labels and title.
        file_name (str): Name of the file to save the plot.

    Returns:
        tuple: (coefficients, equation, R²)
    """
    # Extract x and y values
    x_values = np.array([values['basis_functions'] for values in data.values()])
    y_values = np.array([values[method] for values in data.values()])
    labels = list(data.keys())

    # Create scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(x_values, y_values, color='blue', label='Data Points')

    # Fit a polynomial regression line
    coefficients = np.polyfit(x_values, y_values, degree)
    polynomial = np.poly1d(coefficients)

    # Generate R² value
    y_pred = polynomial(x_values)
    ss_total = np.sum((y_values - np.mean(y_values)) ** 2)
    ss_residual = np.sum((y_values - y_pred) ** 2)
    r_squared = 1 - (ss_residual / ss_total)

    # Generate equation as a string with six decimals for coefficients
    equation_terms = [f"{coeff:.6f} * x^{i}" if i > 0 else f"{coeff:.4f}" 
                      for i, coeff in enumerate(coefficients[::-1])]
    equation = " + ".join(equation_terms)

    # Plot regression line (without adding it to the legend)
    x_range = np.linspace(min(x_values), max(x_values), 100)
    y_fit = polynomial(x_range)
    plt.plot(x_range, y_fit, color='red', linestyle='--', alpha=0.7)

    # Annotate each point with the molecule name
    texts = []
    for i, label in enumerate(labels):
        text = plt.text(x_values[i], y_values[i], label, fontsize=font_size, ha='right')
        texts.append(text)

    # Adjust the text to avoid overlap
    adjust_text(texts, only_move={'points': 'xy', 'texts': 'xy'}, arrowprops=dict(arrowstyle="->", color='red'))

    # Labels and title
    plt.xlabel('Basis Functions', fontsize=font_size+4)
    plt.ylabel(method+' time (min)', fontsize=font_size+4)
    # plt.title(f'Scatter Plot of Basis Functions vs. {method} (Degree {degree} Fit)', fontsize=font_size + 8)

    # Set axis tick font size
    plt.tick_params(axis='both', labelsize=font_size)

    # Show grid
    # plt.grid(True)

    # Add equation text to the plot
    plt.text(0.05, 0.95, f'Equation: {equation}\nR² = {r_squared:.2f}', 
             transform=plt.gca().transAxes, fontsize=font_size, 
             verticalalignment='top', horizontalalignment='left', 
             bbox=dict(facecolor='white', alpha=0.8))

    # Save the plot if requested
    if save_plot:
        plt.savefig(method+"_"+file_name, bbox_inches='tight')

    # Show the plot
    # plt.legend()
    plt.show()

    return coefficients, equation, r_squared


def exponential_function(x, a, b):
    """Exponential function for curve fitting."""
    return a * np.exp(b * x)


def plot_scatter_exponential(data, method, save_plot=False, font_size=14, file_name='scatter_plot.png'):
    """
    Plots a scatter plot of Basis Functions vs. a selected method value and fits an exponential regression line.

    Parameters:
        data (dict): Dictionary containing 'basis_functions' and method values.
        method (str): The key to extract y-values (e.g., 'cpks_time', 'tda_time').
        save_plot (bool): Whether to save the plot as a file.
        font_size (int): Font size for labels and title.
        file_name (str): Name of the file to save the plot.

    Returns:
        tuple: (fitted parameters, equation, R²)
    """
    # Extract x and y values
    x_values = np.array([values['basis_functions'] for values in data.values()])
    y_values = np.array([values[method] for values in data.values()])
    labels = list(data.keys())

    # Ensure y-values are positive for log transformation
    if np.any(y_values <= 0):
        raise ValueError("All y-values must be positive for an exponential fit.")

    # Fit the data to an exponential function
    popt, _ = curve_fit(exponential_function, x_values, y_values, p0=(1, 0.01))
    a, b = popt

    # Compute R² value
    y_pred = exponential_function(x_values, *popt)
    ss_total = np.sum((y_values - np.mean(y_values)) ** 2)
    ss_residual = np.sum((y_values - y_pred) ** 2)
    r_squared = 1 - (ss_residual / ss_total)

    # Generate equation string
    equation = f"y = {a:.6f} * exp({b:.6f} * x)"

    # Create scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(x_values, y_values, color='blue', label='Data Points')

    # Plot exponential regression line
    x_range = np.linspace(min(x_values), max(x_values), 100)
    y_fit = exponential_function(x_range, *popt)
    plt.plot(x_range, y_fit, color='red', linestyle='--', alpha=0.7)

    # Annotate points with molecule names
    texts = []
    for i, label in enumerate(labels):
        text = plt.text(x_values[i], y_values[i], label, fontsize=font_size, ha='right')
        texts.append(text)

    # Adjust text to prevent overlap
    adjust_text(texts, only_move={'points': 'xy', 'texts': 'xy'}, arrowprops=dict(arrowstyle="->", color='red'))

    # Labels and formatting
    plt.xlabel('Basis Functions', fontsize=font_size + 4)
    plt.ylabel(method + ' time (min)', fontsize=font_size + 4)
    plt.tick_params(axis='both', labelsize=font_size)

    # Add equation text
    plt.text(0.05, 0.95, f'Equation: {equation}\nR² = {r_squared:.2f}', 
             transform=plt.gca().transAxes, fontsize=font_size,
             verticalalignment='top', horizontalalignment='left',
             bbox=dict(facecolor='white', alpha=0.8))

    # Save the plot if requested
    if save_plot:
        plt.savefig(method + "_" + file_name, bbox_inches='tight')

    # Show the plot
    plt.show()

    return popt, equation, r_squared


# Take the CPU times in ORCA and QChem
molecule_times_basis = take_times_and_basis(folder_path)

# Creating the new dictionary
cpks_times_dict, tda_times_dict = split_molecule_data(molecule_times_basis)
cpks_minutes_basis = convert_data_to_minutes(cpks_times_dict, 'CPKS')
tda_minutes_basis = convert_data_to_minutes(tda_times_dict, 'TDA')

# Fit and plot for TDA
coefficients_tda, equation_tda, r_squared_tda = plot_scatter(tda_minutes_basis, 'TDA', degree=tda_fit_degree, save_plot=saveplot) 
print('TDA fit:', equation_tda, 'R^2:', r_squared_tda)

# Fit and plot for CPKS
coefficients_cpks, equation_cpks, r_squared_cpks= plot_scatter(cpks_minutes_basis, 'CPKS', degree=cpks_fit_degree, save_plot=saveplot) 
print('CPKS fit:', equation_cpks, 'R^2:', r_squared_cpks)

# Take TDA/CPKS ratio
ratio_dict = compute_ratio(molecule_times_basis)
sorted_molecule_data = dict(sorted(ratio_dict.items(), key=lambda item: item[1]['ratio']))
print(sorted_molecule_data)
