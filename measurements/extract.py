#
# Course: High Performance Computing 2021/2022
#
# Lecturer: Francesco Moscato	fmoscato@unisa.it
#
# Group:
# Marseglia	Mattia		0622701697	    m.marseglia1@studenti.unisa.it
# Spingola     Camilla		0622701698  	c.spingola@studenti.unisa.it
# Turi		    Vito		0622701795  	v.turi3@studenti.unisa.it
# Sica 		Ferdinando	0622701794	    f.sica24@studenti.unisa.it
#
# Copyright (C) 2021 - All Rights Reserved
#
# This file is part of Contest-OMP: RadixSort.
#
# Contest-OMP: RadixSort is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Contest-OMP: RadixSort is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Contest-OMP: RadixSort.  If not, see <http://www.gnu.org/licenses/>.
#

import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import seaborn as sns
from prettytable import PrettyTable
from prettytable import MARKDOWN
from prettytable import MSWORD_FRIENDLY
import re

# set_up of what we want in addiction from every categories of measure
set_up = {
    'init': {

        'jpg': False,
        'speedup': False

    },
    'funct': {

        'jpg': False,
        'speedup': False

    },
    'user': {

        'jpg': False,
        'speedup': False

    },
    'sys': {

        'jpg': False,
        'speedup': False

    },
    'elapsed': {

        'jpg': True,
        'speedup': True

    }
}

# @brief main function that call other to extract all data in different folders
# @param root_algorithm_folders = the root were save all measurements extracted in different configuration


def extraction(root_algorithm_folders):
    # obtain all directory that respect folders_condition for the measurement

    algorithm_folders = [f for f in os.listdir(root_algorithm_folders) if (
        folders_condition(root_algorithm_folders, f, "ALGORITHM-[0-9]"))]
    for algorithm_folder in algorithm_folders:
        print(f"Found Folder : {algorithm_folder}")
        # obtain path were find optimization
        init_path = os.path.join(
            root_algorithm_folders, algorithm_folder)
        extract_from_init_folders(init_path)


# @brief call from main function to extract folders for different sizes
# @param root = the root were save measurements extracted organized for size


def extract_from_init_folders(root_init_folders):
    init_folders = [f for f in os.listdir(root_init_folders) if (
        folders_condition(root_init_folders, f, "INIT_MODE-[0-9]"))]
    for init_folder in init_folders:
        print(f"Found Folder : {init_folder}")
        # obtain path were find optimization
        size_path = os.path.join(root_init_folders, init_folder)
        extract_from_size_folders(size_path)


# @brief call from main function to extract folders for different sizes
# @param root = the root were save measurements extracted organized for size


def extract_from_size_folders(root_size_folders):
    size_folders = [f for f in os.listdir(root_size_folders) if (
        folders_condition(root_size_folders, f, "SIZE-[0-9]"))]
    for size_folder in size_folders:
        print(f"Found Folder : {size_folder}")
        # obtain path were find optimization
        max_digit_path = os.path.join(root_size_folders, size_folder)
        extract_from_max_digit_folders(max_digit_path)

# @brief call from main function to extract folders for different max digit
# @param root = the root were save measurements extracted organized for max digit


def extract_from_max_digit_folders(root_max_digit_folders):
    max_digit_folders = [f for f in os.listdir(root_max_digit_folders) if (
        folders_condition(root_max_digit_folders, f, "MAX_DIGIT-[0-9]"))]
    for max_digit_folder in max_digit_folders:
        print(f"Found Folder : {max_digit_folder}")
        # obtain path were find optimization
        optimization_path = os.path.join(
            root_max_digit_folders, max_digit_folder)
        extract_from_optimization_folders(optimization_path)

# @brief call from main function to extract folders for different optimization
# @param root = the root were save measurements extracted organized for optimization


def extract_from_optimization_folders(root_opt_folders):
    optimization_folders = [f for f in os.listdir(root_opt_folders) if
                            (folders_condition(root_opt_folders, f, "OPTIMIZATION-O[0-9]"))]
    for optimization_folder in optimization_folders:
        print(f"Found Sub-Folder: {optimization_folder}")
        complete_path = os.path.join(root_opt_folders, optimization_folder)
        # call function that extract measurements
        create_table_from_measure(_extract(complete_path), complete_path)

# @brief return a table containing values analized for different configurations
# @param subdir_where_measure = the directories in which found values obtained
# @param complete_path = the path in which save different tables
# CAMBIA QUA!!!!


def create_table_from_measure(subdir_where_measure, complete_path):
    # save header
    header = {
        'values': ['Version', 'Processes', 'Init', 'Real_Func', 'User', 'Sys', 'Elapsed', 'Speedup', 'Efficiency']
    }

    # values in each row
    cells = {'values': []}
    nt = -1
    for filename in subdir_where_measure:
        cell = []
        # obtain all data from filename name
        splitted_filename = filename.split("-")
        splitted_filename = splitted_filename[1].split(".")
        # if was serial
        if "SERIAL" in filename:
            seq = subdir_where_measure[filename]['elapsed']
            nt = 1
            cell.append('Serial')
            cell.append(nt)
        else:
            nt = int(splitted_filename[0])
            cell.append('Parallel')
            cell.append(nt)

        for col in set_up:
            cell.append(subdir_where_measure[filename][col])
            if set_up[col]['speedup']:
                speedup, efficiency = _compute_speedup(
                    seq, subdir_where_measure[filename][col], nt)
                cell.append(speedup)
                cell.append(efficiency)
        cells['values'].append(cell)

    table_filename = complete_path + "/summary_table.csv"
    plot_filename = complete_path + "/speedup-plot.jpg"

    table = _make_table(header['values'], cells['values'], name=table_filename)
    _plot_from_table(header["values"], cells["values"], name=plot_filename)

# @brief return all measures for all configurations obtained by all data captured
# @param path_were_extract = the root were extract all measures


def _extract(path_were_extract):
    # save the path were i'm now
    prev = os.getcwd()
    # change my directory into path_were_extract
    os.chdir(path_were_extract)

    # List of file from directory who are csv files
    filenames = [f for f in os.listdir('.') if file_condition(f)]
    if not os.path.exists("jpg"):
        os.mkdir("jpg")

    # order files in the list filenames
    filenames = sorted(filenames)
    all_measure_for_all_configuration = {}

    for filename in filenames:
        effective_time_with_this_configuration = {}
        print('Processing : ' + filename)
        # i read the "filename" that contains all measures
        ds = pd.read_csv(filename)

        # I scroll through all the categories
        for cat in set_up.keys():
            print('Processing : ' + filename + ", Col : " + cat)
            # extract the measure for the specific category
            x_data = ds[cat]

            # calculate the mean and standard variation from the Gaussian distribution of data
            mean, std = stats.norm.fit(x_data)

            # consider in x_data only data that are not busted
            # 68,3% = P{ μ − 1,00 σ < X < μ + 1,00 σ }
            x_data = ds[(ds[cat] <= (mean + std)) &
                        (ds[cat] >= (mean - std))][cat]

            # calculate more accurate mean and standard variation from the Gaussian distribution of data
            # minus data busted
            mean, std = stats.norm.fit(x_data)

            # insert the mean in the Set that contains measure for that file
            effective_time_with_this_configuration[cat] = mean
            # if category have set 'jpg' as elapsed it save an histogram of value in a directory called 'jpg'

            if set_up[cat]['jpg']:

                sns.histplot(x_data, kde=False)
                plt.savefig("jpg/" + str(cat) + "_" +
                            filename.split('.')[0] + ".jpg")
                plt.close()
        # save alla measure of file in this variable at the end of the process it have all measure of all files
        all_measure_for_all_configuration[filename] = effective_time_with_this_configuration

    # return to the folder from before
    os.chdir(prev)
    # return all correct measure of all files in this folder
    return all_measure_for_all_configuration

# @brief main function to compute speedup for values obtained
# @param t = time for serial algorithm
# @param tp = time for parallel algorithm
# @param nt = number of threads


def _compute_speedup(t, tp, nt):
    speedup = t / tp
    efficiency = t / (tp * float(nt))
    return speedup, efficiency


# @brief folders_condition verify condition to find folders were save measure files
# @param root = the root were check subdirectories
# @param subdirectories = specific subdirectories to check
# @param string = specific string to check
def folders_condition(root, subdirectories, string):
    if os.path.isdir(os.path.join(root, subdirectories)) and re.match(string, subdirectories):
        return True
    return False


# @brief file_condition verify condition to find file were save measure
# @param file = specific subdirectories to check
def file_condition(file):
    if os.path.isfile(file) and file.endswith(".csv") and re.match("N_PROCESS-[0-9]", file):
        return True
    return False

# @brief produce a table from data
# @param header = name of field
# @param rows = number of elements to table
# @param print_table = boolean that explain if print the table
# @param save = boolean that explain if save the table
# @param name = pathe where save the table


def _make_table(header, rows, print_table=False, save=True, name=""):
    if save and not name:
        raise Exception("No filename to save file")
    x = PrettyTable()
    x.field_names = header
    rows.sort(key=lambda x: x[1])
    x.add_rows(rows)
    if save:
        _save_table(x, name)
    if print_table:
        print(x)
    return x

# @brief produce a plot from data
# @param header = name of axis
# @param rows = number of elements to plot
# @param save = boolean that explain if save the plot
# @param name = pathe where save the plot
# @param show_plot = boolean that explain if show the plot


def _plot_from_table(header, rows, save=True, name="", show_plot=False):
    if save and not name:
        raise Exception("No filename to save file")

    x = [0]
    y = [0]
    speedup_pos = header.index("Speedup")
    thread_pos = header.index("Processes")
    for row in rows[1:]:
        x.append(row[thread_pos])
        y.append(row[speedup_pos])

    x_th = np.array(x)
    fig, ax = plt.subplots(figsize=(15, 8))
    ax.plot(x_th, y, 'ro-', label='Experimental')
    ax.plot(x_th, x_th, color='blue', label='Ideal')
    # same as y_th, bisection
    plt.style.use('seaborn-whitegrid')

    plt.xticks(x)
    for i, j in zip(x, y):
        ax.annotate(str(round(j, 3)), xy=(i, j))
    # ax.set_aspect(1/2)

    plt.autoscale(enable=True, axis='x', tight=True)
    plt.autoscale(enable=True, axis='y', tight=True)
    tableName = name[name.find("ALGORITHM"):-17].replace("/", "  ")

    plt.title(tableName, fontdict={
              'family': 'serif', 'color': 'black', 'weight': 'bold', 'size': 15})
    plt.legend()
    plt.xlabel("Processors", labelpad=10)
    plt.ylabel("Speedup")
    if show_plot:
        plt.show()
    if save:
        plt.savefig(name)

    plt.close()

# @brief save a table in a file
# @param table = reference to table
# @param filename = reference to file


def _save_table(table, filename):
    with open(filename, "w") as table_file:
        # table.set_style(MARKDOWN)
        table.set_style(MSWORD_FRIENDLY)
        data = table.get_string()
        table_file.write(data)


if __name__ == "__main__":
    # root_for_extraction = the name of directory where is this file by adding a subdirectory measure
    root_for_extraction = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "measures/")
    extraction(root_for_extraction)
