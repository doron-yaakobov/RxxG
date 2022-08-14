# index{R} = -4; index{x} = -3; index{x} = -2 ; index{G} = -1
import csv
from collections import OrderedDict
from datetime import datetime

import matplotlib.pylab as plt
import numpy as np
import pandas as pd


def calculate_aa_freq(data: pd.DataFrame) -> dict:
    first_x_population, sec_x_population, pairs_incidence = {}, {}, {}
    first_x_aa, first_x_incidences, sec_x_aa, sec_x_incidences, first_x_aa_freq, sec_x_aa_freq = (), (), (), (), (), ()
    is_relevant_aa_seq, is_found_rxxg_pattern = bool(), bool()
    rxxg_incidence: int = 0
    first_x, sec_x, pair = "", "", ""

    for _, row in data.iterrows():
        is_found_rxxg_pattern = False
        row = tuple(row)  # fixes pd error
        for idx, letter in enumerate(row, start=0):
            is_relevant_aa_seq = letter in ["R", "r"] and (idx + 3) < len(row) and row[idx + 3] in ["G", "g"]
            if is_relevant_aa_seq:
                is_found_rxxg_pattern = True
                first_x = row[idx + 1]
                sec_x = row[idx + 2]
                pair = f"{first_x}{sec_x}"
        # Makes sure we are updating only the last RxxG in a row:
        if is_found_rxxg_pattern:
            first_x_population[first_x] = first_x_population.get(first_x, 0) + 1
            sec_x_population[sec_x] = sec_x_population.get(sec_x, 0) + 1
            pairs_incidence[pair] = pairs_incidence.get(pair, 0) + 1

    # extract incidences:
    first_x_aa, first_x_incidences = zip(*sorted(first_x_population.items()))
    sec_x_aa, sec_x_incidences = zip(*sorted(sec_x_population.items()))
    for pair in pairs_incidence.items():
        rxxg_incidence += pair[1]

    # calculate experiment's freqs:
    first_x_aa_freq = dict(zip(first_x_aa,
                               tuple(map(lambda aa_incidence: aa_incidence / rxxg_incidence, first_x_incidences))))
    sec_x_aa_freq = dict(zip(sec_x_aa,
                             tuple(map(lambda aa_incidence: aa_incidence / rxxg_incidence, sec_x_incidences))))
    return {-3: first_x_aa_freq, -2: sec_x_aa_freq}


def update_missing_aa(data: dict, delighted_keys: list) -> dict:
    data_aa = list(data.keys())
    data_missing_aa = list(set(data_aa).symmetric_difference((set(delighted_keys))))
    for aa in data_missing_aa:
        data[aa] = float(0)
    return dict(OrderedDict(sorted(data.items())))


def plot_frequencies(experiment_frequencies: dict, population_frequencies: dict, index: int) -> None:
    # get \ arrange plot data:
    x_aa = list(population_frequencies[index].keys())
    y_pop = list(population_frequencies[index].values())
    exp_aa_freq[index] = update_missing_aa(experiment_frequencies[index], x_aa)
    z_exp = list(experiment_frequencies[index].values())

    # plot
    x_axis = np.arange(len(x_aa))
    plt.bar(x_axis - 0.2, y_pop, 0.4, label='Population Frequency')
    plt.bar(x_axis + 0.2, z_exp, 0.4, label='Experiment Frequency')
    plt.xticks(x_axis, x_aa)
    plt.xlabel("Amino Acids")
    plt.ylabel("Frequency")
    plt.suptitle("Amino Acids Population VS Experiment Frequencies")
    plt.title(f"Index = {index}")
    plt.legend()

    # log results on csv:
    date = datetime.now().strftime("%Y.%m.%d_%H-%M-%S")
    data = [["AminoAcid", "Population", "Experiment"]]
    for row in zip(x_aa, y_pop, z_exp):
        data.append(list(row))
    file_name = f"{date}Amino Acids Population VS Experiment Frequencies, index = {index}.csv"

    with open(file_name, 'w', newline='') as results_csv:
        writer = csv.writer(results_csv)
        writer.writerows(data)

    csv_data = pd.read_csv(file_name)
    csv_data.to_excel(file_name.replace("csv", "xlsx"), index=None, header=True)
    plt.savefig(f'{date}Amino Acids Population VS Experiment Frequencies, index = {index}.png')
    plt.show()


# region setup
EXPERIMENTS_OUTPUT = "Experiment's Results (Candidates).xlsx"
EXPERIMENTS_OUTPUT_SHEET_NAME = "23mer_CRL_sub"
POPULATION_PROTEINS = "Population Proteins (C23 Library).xlsx"
POPULATION_PROTEINS_SHEET_NAME = 'גיליון1'
ACHILLES_DATA = "Ubiquitome_AchillesData.xlsx"

exp_out: pd.DataFrame = pd.read_excel(EXPERIMENTS_OUTPUT, None)[EXPERIMENTS_OUTPUT_SHEET_NAME]
pop_proteins: pd.DataFrame = pd.read_excel(POPULATION_PROTEINS, None)[POPULATION_PROTEINS_SHEET_NAME]
# endregion

exp_aa_freq = calculate_aa_freq(exp_out)
pop_aa_freq = calculate_aa_freq(pop_proteins)

# debugging values:
# exp_aa_freq = {
#     -3: {'A': 0.09090909090909091, 'C': 0.030303030303030304, 'E': 0.09090909090909091, 'F': 0.06060606060606061,
#          'G': 0.12121212121212122, 'H': 0.030303030303030304, 'L': 0.12121212121212122, 'M': 0.09090909090909091,
#          'R': 0.15151515151515152, 'S': 0.18181818181818182, 'Y': 0.030303030303030304},
#     -2: {'A': 0.09090909090909091, 'C': 0.09090909090909091, 'E': 0.030303030303030304, 'G': 0.030303030303030304,
#          'H': 0.030303030303030304, 'K': 0.030303030303030304, 'L': 0.18181818181818182, 'M': 0.06060606060606061,
#          'N': 0.06060606060606061, 'Q': 0.030303030303030304, 'R': 0.06060606060606061, 'S': 0.09090909090909091,
#          'T': 0.09090909090909091, 'W': 0.06060606060606061, 'Y': 0.06060606060606061}}
# pop_aa_freq = {
#     -3: {'A': 0.08857938718662953, 'C': 0.019498607242339833, 'D': 0.03565459610027855, 'E': 0.05403899721448468,
#          'F': 0.02395543175487465, 'G': 0.10250696378830083, 'H': 0.025069637883008356, 'I': 0.040668523676880224,
#          'K': 0.06852367688022284, 'L': 0.09637883008356546, 'M': 0.017270194986072424, 'N': 0.027855153203342618,
#          'P': 0.07910863509749304, 'Q': 0.03231197771587744, 'R': 0.08022284122562674, 'S': 0.07855153203342619,
#          'T': 0.05236768802228412, 'V': 0.05125348189415042, 'W': 0.010584958217270195, 'Y': 0.015598885793871866},
#     -2: {'A': 0.07409470752089137, 'C': 0.028969359331476322, 'D': 0.03454038997214485, 'E': 0.05515320334261838,
#          'F': 0.025069637883008356, 'G': 0.09637883008356546, 'H': 0.022841225626740947, 'I': 0.024512534818941504,
#          'K': 0.05181058495821727, 'L': 0.10362116991643454, 'M': 0.017270194986072424, 'N': 0.029526462395543174,
#          'P': 0.08913649025069638, 'Q': 0.05181058495821727, 'R': 0.08969359331476323, 'S': 0.08022284122562674,
#          'T': 0.03955431754874652, 'V': 0.04233983286908078, 'W': 0.020612813370473538, 'Y': 0.022841225626740947}}

plot_frequencies(exp_aa_freq, pop_aa_freq, index=-2)
plot_frequencies(exp_aa_freq, pop_aa_freq, index=-3)
