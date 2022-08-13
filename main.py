# index{R} = -4; index{x} = -3; index{x} = -2 ; index{G} = -1

import gc

import matplotlib.pylab as plt
import pandas as pd

EXPERIMENTS_OUTPUT = "Experiment's Results (Candidates).xlsx"
POPULATION_PROTEINS = "Population Proteins (C23 Library).xlsx"
ACHILLES_DATA = "Ubiquitome_AchillesData.xlsx"

# region find experimantal aa's frequencies, with respect to their index:
exp_out: pd.DataFrame = pd.read_excel(EXPERIMENTS_OUTPUT, None)["23mer_CRL_sub"]

first_x_population, sec_x_population, pairs_incidence = {}, {}, {}
# aa stands for "Amino Acid"
first_x_aa, first_x_incidences, sec_x_aa, sec_x_incidences, first_x_aa_freq, sec_x_aa_freq = (), (), (), (), (), ()
is_relevant_aa_seq, is_found_rxxg_pattern = bool(), bool()
rxxg_incidence: int = 0
first_x, sec_x, pair = "", "", ""

for _, row in exp_out.iterrows():
    is_found_rxxg_pattern = False
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
# cleanup namespace:
del exp_out, first_x, first_x_incidences, first_x_population, idx, is_relevant_aa_seq, letter, pair, pairs_incidence, \
    row, rxxg_incidence, sec_x, sec_x_incidences, sec_x_population, sec_x_aa, first_x_aa, _
gc.collect()
# endregion


# col_vals = plt.figure(figsize=(17, 9.5))
# plt.suptitle(sheet_title)
# plt.subplot(211)
# plt.bar(first_x_aa, first_x_incidences, color='maroon', width=0.3)
# plt.ylabel("Incidence [Index -3]")
#
# # sec_col_res_list = sec_x_population.items()
# # sec_col_res_list = sorted(sec_col_res_list)
# # amino_acid, val = zip(*sec_col_res_list)
# #
# # plt.subplot(212)
# # plt.bar(amino_acid, val, color='maroon', width=0.3)
# # plt.xlabel("Amino Acid")
# # plt.ylabel("Incidence [Index -2]")
# #
# # pairs_temp = plt.figure(figsize=(17, 9.5))
# # plt.suptitle(sheet_title)
# # pairs_res_list = pairs.items()
# # pairs_res_list = sorted(pairs_res_list)
# # pair, val = zip(*pairs_res_list)
# #
# # plt.bar(pair, val, color='maroon', width=0.3)
# # plt.xlabel("Pair [Amino acids]")
# # plt.ylabel("Incidence")

plt.show()
