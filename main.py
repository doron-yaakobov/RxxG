import matplotlib.pylab as plt
import pandas as pd
from attrdict import AttrDict

# user input:
CSV_PATH = r"~/private/Din's_RXXG/candidates.xlsx"

csv = AttrDict(pd.read_excel(CSV_PATH, None))

for sheet_title in csv.keys():
    df = csv[sheet_title]
    # first -> -3 ; second -> -2 ; G -> -1
    first_col_res, sec_col_res, freq_pairs = {}, {}, {}

    for idx, row in df.iterrows():
        for inner_idx, letter in enumerate(row, start=0):
            is_relevant_amino_acids_seq = letter in ["R", "r"] \
                                          and (inner_idx + 3) < len(row) and row[inner_idx + 3] in ["G", "g"]
            if is_relevant_amino_acids_seq:
                # get amino acids value:
                first_ch = row[inner_idx + 1]
                sec_ch = row[inner_idx + 2]
                pair = f"{first_ch}{sec_ch}"

                # update counters:
                first_col_res[first_ch] = first_col_res.get(first_ch, 0) + 1
                sec_col_res[sec_ch] = sec_col_res.get(sec_ch, 0) + 1
                freq_pairs[pair] = freq_pairs.get(pair, 0) + 1

    first_col_res_list = first_col_res.items()
    first_col_res_list = sorted(first_col_res_list)
    amino_acid, val = zip(*first_col_res_list)

    col_vals = plt.figure(figsize=(17, 9.5))
    plt.suptitle(sheet_title)
    plt.subplot(211)
    plt.bar(amino_acid, val, color='maroon', width=0.3)
    plt.ylabel("Incidence [Index -3]")

    sec_col_res_list = sec_col_res.items()
    sec_col_res_list = sorted(sec_col_res_list)
    amino_acid, val = zip(*sec_col_res_list)

    plt.subplot(212)
    plt.bar(amino_acid, val, color='maroon', width=0.3)
    plt.xlabel("Amino Acid")
    plt.ylabel("Incidence [Index -2]")

    pairs = plt.figure(figsize=(17, 9.5))
    plt.suptitle(sheet_title)
    pairs_res_list = freq_pairs.items()
    pairs_res_list = sorted(pairs_res_list)
    pair, val = zip(*pairs_res_list)

    plt.bar(pair, val, color='maroon', width=0.3)
    plt.xlabel("Pair [Amino acids]")
    plt.ylabel("Incidence")

plt.show()
