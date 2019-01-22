#!/usr/bin/env python3

import os
import sys
import argparse
import random
import numpy as np
import scipy
import pandas as pd
import statsmodels.api as sm


# prepare column names and indexes
def prepare_cols_indx(alignment_file):
    alignment = open(alignment_file, "r").readlines()
    ali_len = len(alignment[1])
    column_names = []
    indexes = ['start1', 'end1']
    for i in range(0, len(alignment)-1, 2):
        protein_name = alignment[i].split(" ")[0].lstrip(">")
        if protein_name.startswith("ENS"):
            protein_name = protein_name.split(".")[0]
        start = alignment[i].split(";")[-2].strip()
        column_names.append(protein_name+"+"+start)
    for i in list(range(0, ali_len-1)):
        indexes.append(i)
    return(column_names, indexes)


# dataframe with alignment as letters
def letter_ali_dataframe(alignment_file, column_names, indexes):
    alignment = open(alignment_file, "r").readlines()
    ali_len = len(alignment[1])
    alignment_list = []
    start_row = []
    end_row = []
    for i in range(0, len(alignment)-1, 2):
        start = int(alignment[i].split(";")[1].strip())
        end = int(alignment[i].split(";")[2].strip())
        start_row.append(start)
        end_row.append(end)
    alignment_list.append(start_row)
    alignment_list.append(end_row)
    for j in range(0, ali_len-1):
        letter_row = []
        for i in range(1, len(alignment), 2):
            letter_row.append(alignment[i][j])
        alignment_list.append(letter_row)
    letter_alignment = pd.DataFrame(
        alignment_list,
        columns=column_names,
        index=indexes,
    )
    return(letter_alignment)


# dataframe with alignment as positions
def pos_dataframe(letter_alignment, column_names, indexes):
    position_alignment = pd.DataFrame(index=indexes[2:])
    for protein in column_names:
        pos_seq = []
        start = letter_alignment.loc['start1', protein]
        seq = letter_alignment[protein].tolist()[2:]
        counter = 0
        for aa in seq:
            if aa == "-":
                pos_seq.append(0)
            else:
                pos_seq.append(start+counter)
                counter += 1
        position_alignment[protein] = pos_seq
    return(position_alignment)


# 0/1 dataframe with mapped phosphorylations onto alignment
def phos_dataframe(phosp_info, letter_alignment, position_alignment,
                   column_names, indexes):
    phosp_alignment = pd.DataFrame(0, index=indexes[2:], columns=column_names)
    for column in column_names:
        protein_id = column.split("+")[0].strip()
        start_no = int(column.split("+")[1].strip())
        for i in range(0, len(phosp_info)):
            line = phosp_info[i].split(",")
            protein_id2 = line[1].strip()
            phosp_pos = int(line[2].strip())
            phosp_aa = line[3].strip()
            if protein_id == protein_id2 and start_no <= phosp_pos:
                row_with_pos = position_alignment.index[
                    position_alignment[column] == phosp_pos,
                ].tolist()
                if len(row_with_pos) == 1:
                    aa_in_pos = letter_alignment.loc[row_with_pos[0], column]
                    if phosp_aa == aa_in_pos:
                        phosp_alignment.at[row_with_pos[0], column] = 1
    return(phosp_alignment)


# rolling window for a list (returns list [2:-2] of original length)
def count_window_for_list(list):
    window_values = []
    for i in range(2, len(list)-2):
        a = int(list[i-2])
        b = int(list[i-1])
        c = int(list[i])
        d = int(list[i+1])
        e = int(list[i+2])
        how_many = a+b+c+d+e
        bg = float(how_many/5.)
        window_values.append(bg)
    return(window_values)


# background construction
def make_random_histogram(elements, size):
    my_list = [1]*elements + [0]*(size-elements)
    random.shuffle(my_list)
    return(my_list)


# returns one permutated column of rolled window on background
def permutated_dataframe(letter_alignment, phosp_alignment, column_names,
                         indexes):
    permutated_alignment = pd.DataFrame(
        0,
        index=indexes[2:],
        columns=column_names,
    )
    for protein_id in column_names:
        rows_with_S = letter_alignment.index[
            letter_alignment[protein_id] == "S",
        ].tolist()
        rows_with_T = letter_alignment.index[
            letter_alignment[protein_id] == "T",
        ].tolist()
        rows_with_Y = letter_alignment.index[
            letter_alignment[protein_id] == "Y",
        ].tolist()
        phosp_rows = phosp_alignment.index[
            phosp_alignment[protein_id] >= 1,
        ].tolist()
        pS = pT = pY = 0
        for i in phosp_rows:
            if i in rows_with_S:
                pS += 1
            if i in rows_with_T:
                pT += 1
            if i in rows_with_Y:
                pY += 1
        random_S = make_random_histogram(pS, len(rows_with_S))
        random_T = make_random_histogram(pT, len(rows_with_T))
        random_Y = make_random_histogram(pY, len(rows_with_Y))
        for i, j in zip(random_S, rows_with_S):
            permutated_alignment.at[j, protein_id] = i
        for i, j in zip(random_T, rows_with_T):
            permutated_alignment.at[j, protein_id] = i
        for i, j in zip(random_Y, rows_with_Y):
            permutated_alignment.at[j, protein_id] = i
    sum_of_phosps = permutated_alignment.loc[0:].sum(axis=1).tolist()
    one_of_bg = count_window_for_list(sum_of_phosps)
    return(one_of_bg)


def count_zscore(fg, mean, stdev):
    if stdev == 0:
        zscore = 0
    else:
        zscore = (fg-mean)/stdev
    return(zscore)


def count_pval(fg, mean, stdev):
    fg = float(fg)
    mean = float(mean)
    stdev = float(stdev)
    if fg > mean:
        z = count_zscore(fg, mean, stdev)
        p_val = scipy.stats.norm.sf(abs(z))
    if fg <= mean:
        z = 1.
        p_val = 1.
    return(p_val)


def getHotspotSitesInFile(alignment_path, ptms, repeat):
    # Starting dataframes
    (column_names, indexes) = prepare_cols_indx(alignment_path)
    letter_alignment = letter_ali_dataframe(alignment_path, column_names,
                                            indexes)
    position_alignment = pos_dataframe(letter_alignment, column_names, indexes)
    phosp_alignment = phos_dataframe(ptms, letter_alignment,
                                     position_alignment, column_names, indexes)
    # Foreground dataframe
    sum_of_phosps = phosp_alignment.loc[0:].sum(axis=1).tolist()
    foreground = count_window_for_list(sum_of_phosps)
    # dataframe (and permutations)
    bg_df = pd.DataFrame(index=indexes[4:-2], columns=list(range(0, repeat)))
    for k in range(0, repeat):
        bg_df[k] = permutated_dataframe(letter_alignment, phosp_alignment,
                                        column_names, indexes)
    bg_df['domain'] = os.path.splitext(os.path.basename(alignment_path))[0]
    bg_df['position_aln'] = indexes[4:-2]
    bg_df['bg_medians'] = bg_df.iloc[:, 0:(repeat-1)].median(axis=1).tolist()
    bg_df['bg_stdev'] = bg_df.iloc[:, 0:(repeat-1)].std(axis=1).tolist()
    bg_df['bg_means'] = bg_df.iloc[:, 0:(repeat-1)].mean(axis=1).tolist()
    bg_df['foreground'] = foreground
    # this is to count pvals and add them to the column in background_dataframe
    all_pvals = []
    for row in bg_df.index:
        fg = bg_df.at[row, 'foreground']
        mean = bg_df.at[row, 'bg_means']
        st_dev = bg_df.at[row, 'bg_stdev']
        p_val = count_pval(fg, mean, st_dev)
        all_pvals.append(p_val)
    bg_df['pvals'] = all_pvals
    bg_df['pvals'] = bg_df.loc[:, ('pvals')] + 1e-150  # pseudocount for 0s
    bg_dataframe = bg_df.loc[:, ("domain", "position_aln", "foreground",
                                 "pvals")].copy()
    # adding protein/position level information
    let_copy = letter_alignment[2:].copy()
    let_copy["position_aln"] = let_copy.index
    position_alignment["position_aln"] = position_alignment.index
    tomerge = pd.merge(let_copy.infer_objects().melt(id_vars=["position_aln"],
                                                     var_name="protein",
                                                     value_name="residue"),
                       position_alignment.melt(id_vars=["position_aln"],
                                               var_name="protein",
                                               value_name="position"),
                       on=["position_aln", "protein"])
    tomerge["protein"] = tomerge["protein"].str.split("+", expand=True)[0]
    expanded_df = bg_dataframe.merge(tomerge, on="position_aln")
    # ## fix position in alignment
    expanded_df["position_aln"] = expanded_df[["position_aln"]] + 1
    return(expanded_df)


def multipletesting(h_sites):
    tmp = h_sites.loc[:, ("domain", "position_aln", "pvals")] \
                 .drop_duplicates()
    tmp["p_adjust"] = sm.stats.multipletests(tmp.pvals.tolist(),
                                             method='b')[1].tolist()
    h_sites = h_sites.merge(tmp,
                            on=["domain", "position_aln", "pvals"])
    return(tmp)


def get_hotspot_regions(h_sites):
    # Creting hotspot regions (+/-2)
    #  Positions + neighboring positions and whether they are hotspots or not
    df = h_sites.loc[:, ("domain", "protein", "position", "position_aln",)] \
                .drop_duplicates().copy()
    df["minus2"] = df.position_aln - 2
    df["minus1"] = df.position_aln - 1
    df["exact"] = df.position_aln
    df["plus1"] = df.position_aln + 1
    df["plus2"] = df.position_aln + 2
    df = pd.melt(df,
                 id_vars=["domain", "protein", "position", "position_aln"],
                 var_name="rel_position",
                 value_name="neigh_position")
    df = df.drop(columns=["rel_position"])
    # hotspot-domain definitions
    h_dom = pd.merge(df.loc[:, ("domain", "position_aln", "neigh_position")]
                     .copy().drop_duplicates(),
                     h_sites.loc[:, ("domain", "position_aln", "hotspot")]
                     .drop_duplicates()
                     .rename(columns={"position_aln": "neigh_position"}),
                     on=["domain", "neigh_position"])
    h_dom = h_dom.groupby(["domain", "position_aln"])[["hotspot"]] \
                 .any() \
                 .reset_index()
    h_dom = h_dom.rename(columns={"hotspot": "in_hotspot_region"})
    h_dom["in_hotspot_region_inv"] = np.invert(h_dom.in_hotspot_region)
    h_dom = h_dom.sort_values(by=["domain", "position_aln"])
    h_dom["in_hotspot_region_cumsum"] = h_dom \
        .groupby(["domain"])[["in_hotspot_region_inv"]] \
        .cumsum()
    h_def = h_dom\
        .groupby(["domain"])["in_hotspot_region_cumsum"] \
        .value_counts() \
        .rename("counts") \
        .reset_index()
    h_def = h_def[h_def.counts > 1].drop(columns="counts")
    h_def = pd.merge(h_dom,
                     h_def,
                     on=["domain", "in_hotspot_region_cumsum"])
    h_def = h_def[h_def.in_hotspot_region]
    h_def["start_aln"] = h_def \
        .groupby(["domain", "in_hotspot_region_cumsum"])["position_aln"] \
        .transform("min")
    h_def["end_aln"] = h_def \
        .groupby(["domain", "in_hotspot_region_cumsum"])["position_aln"] \
        .transform("max")
    h_def["hotspot_id"] = h_def["domain"]\
        .str.cat(h_def.start_aln.astype(str), sep="_")
    h_def = h_def.loc[:, ("domain", "position_aln", "start_aln", "end_aln",
                          "hotspot_id")]
    # add minimum p_adjust
    h_def = pd.merge(h_def,
                     h_sites.loc[:, ("domain", "position_aln", "p_adjust")]
                     .drop_duplicates(),
                     on=["domain", "position_aln"])
    h_def['min_adj_pval'] = h_def \
        .groupby(["domain", "hotspot_id"])["p_adjust"] \
        .transform("min")
    h_def = h_def.drop(columns=["p_adjust"])
    return(h_def)


def f(x):
    d = {}
    d['sequence'] = x["residue"].str.cat()
    d['start'] = x[x.residue != "-"]["position"].min()
    d['end'] = x[x.residue != "-"]["position"].max()
    d['foreground_max'] = x['foreground'].max()
    return pd.Series(d, index=["start", "end", "sequence", "foreground_max"])


def find_hotspot_instances(hotspot_sites, hotspot_definitions):
    if(len(hotspot_definitions) == 0):
        return()
    # instances of the hotspots in the proteins
    h_inst = pd.merge(hotspot_sites,
                      hotspot_definitions,
                      on=["domain", "position_aln"],
                      how="left")
    h_inst["hotspot_region"] = np.invert(h_inst[["start_aln"]].isna())
    h_inst["hotspot_region_inv"] = h_inst[["start_aln"]].isna()
    h_inst = h_inst.loc[h_inst['residue'] != "-"]
    h_inst = h_inst \
        .sort_values(by=["domain", "hotspot_id", "protein", "position"])
    out = h_inst \
        .groupby(["domain", "protein", "hotspot_id", "start_aln", "end_aln",
                  "min_adj_pval"]) \
        .apply(f).reset_index()
    out["start"] = out["start"].astype("int")
    out["end"] = out["end"].astype("int")
    out["start_aln"] = out["start_aln"].astype("int")
    out["end_aln"] = out["end_aln"].astype("int")
    out["foreground_max"] = out["foreground_max"].astype("int")
    return(out)


###########
# MAIN
###########
if __name__ == '__main__':
    # Reading arguments
    parser = argparse \
        .ArgumentParser(
            description='Estimate PTM hotspots in sequence alignments')
    parser.add_argument('--dir',
                        nargs='?',
                        default="db/alignments",
                        action="store",
                        metavar="PATH",
                        dest="ALIGNMENTS_DIR",
                        help='fasta alignments dir (default: db/alignments)')
    parser.add_argument('--ptmfile',
                        nargs='?',
                        default="db/all_phosps",
                        action="store",
                        metavar="PATH",
                        dest="PTM_FILE",
                        help='file containing PTMs (default: db/all_phosps)')
    parser.add_argument('-d',
                        '--domain',
                        nargs='?',
                        action="store",
                        metavar="PFXXXXX",
                        dest="domain",
                        help='query single domain (i.e. PF00069)')
    parser.add_argument('--iter',
                        nargs='?',
                        default=100,
                        action="store",
                        dest="ITER",
                        metavar="INTEGER",
                        type=int,
                        help='number of permutations (default: 100)')
    parser.add_argument('--threshold',
                        nargs='?',
                        default=0.01,
                        action="store",
                        metavar="FLOAT",
                        dest="THRESHOLD",
                        type=float,
                        help='Corrected p-value threshold (default: 0.01)')
    parser.add_argument('--foreground',
                        nargs='?',
                        default=2,
                        action="store",
                        metavar="FLOAT",
                        dest="FORE_VAL",
                        type=float,
                        help='effect-size foreground cutoff (default: 2)')
    parser.add_argument('-o',
                        '--out',
                        required=True,
                        action="store",
                        metavar="PATH",
                        dest="OUTPUTFILE",
                        help='output csv file')
    parser.add_argument('--printSites',
                        dest="PRINTSITES",
                        action="store_true",
                        help='print all residue predictions',
                        default=False)

    results = parser.parse_args()
    ALIGNMENTS_DIR = results.ALIGNMENTS_DIR
    PTM_FILE = results.PTM_FILE
    ITER = results.ITER
    THRESHOLD = results.THRESHOLD
    FORE_VAL = results.FORE_VAL
    OUTPUTFILE = results.OUTPUTFILE
    PRINTSITES = results.PRINTSITES

    # reading all PTMs
    ptms = open(PTM_FILE, "r").readlines()
    # check if domain is specified. read all otherwise
    if results.domain:
        alignmentFiles = [results.domain + ".fasta"]
        if not os.path.isfile(os.path.join(ALIGNMENTS_DIR, alignmentFiles[0])):
            sys.stderr.write("Domain file does not exist!")
            sys.exit(1)
            # read all
    else:
        alignmentFiles = os.listdir(ALIGNMENTS_DIR)

    # estimating all hotspot sites
    allHotspots = []
    for filename in alignmentFiles:
        print("* " + filename)
        alignment_path = os.path.join(ALIGNMENTS_DIR, filename)
        expanded_dataframe = getHotspotSitesInFile(alignment_path, ptms, ITER)
        allHotspots.append(expanded_dataframe)
    h_sites = pd.concat(allHotspots)
    h_sites = h_sites.loc[:, ("domain", "protein", "position", "residue",
                              "position_aln", "foreground", "pvals")] \
                     .copy()
    h_sites = multipletesting(h_sites)
    h_sites["hotspot"] = (h_sites.p_adjust <= THRESHOLD) & \
        (h_sites.foreground >= FORE_VAL)

    # estimating hotspot ranges
    # All hotspot range-definitions
    hotspot_definitions = get_hotspot_regions(h_sites)
    # match of all hotspots into
    hotspot_regions = find_hotspot_instances(h_sites, hotspot_definitions)

    print("Done!")

    if PRINTSITES:
        # h_sites[h_sites.hotspot].to_csv(OUTPUTFILE, index = False)
        h_sites.to_csv(OUTPUTFILE, index=False)
    else:
        if(len(hotspot_regions) == 0):
            thisfile = open(OUTPUTFILE, 'w')
            thisfile.write('0 hotspots found\n')
            thisfile.close()
        else:
            hotspot_regions.to_csv(OUTPUTFILE, index=False)
