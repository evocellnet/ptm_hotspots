#!/usr/bin/env python3

from sys import argv
import os
import sys
import argparse
import random
import numpy as np
import scipy
from scipy.stats import norm
import pandas as pd
from itertools import chain 
import statsmodels.api as sm

## prepare column names and indexes
def prepare_cols_indx(alignment_file):
    alignment=open(alignment_file,"r").readlines()
    ali_len=len(alignment[1])
    column_names=[]
    indexes=['start1','end1']
    for i in range(0,len(alignment)-1,2):
        protein_name=alignment[i].split(" ")[0].lstrip(">")
        start=alignment[i].split(";")[-2].strip()
        column_names.append(protein_name+"+"+start)
    for i in list(range(0,ali_len-1)):
        indexes.append(i)
    return(column_names, indexes)       

## dataframe with alignment as letters
def letter_ali_dataframe(alignment_file, column_names, indexes):
    alignment=open(alignment_file,"r").readlines()
    ali_len=len(alignment[1])
    alignment_list=[]
    start_row=[]
    end_row=[]
    for i in range(0,len(alignment)-1,2):
        start=int(alignment[i].split(";")[1].strip())
        end=int(alignment[i].split(";")[2].strip())
        start_row.append(start)
        end_row.append(end)
    alignment_list.append(start_row)
    alignment_list.append(end_row)
    for j in range(0,ali_len-1):
        letter_row=[]
        for i in range(1,len(alignment),2):
            letter_row.append(alignment[i][j])
        alignment_list.append(letter_row)
    letter_alignment=pd.DataFrame(alignment_list, columns=column_names, index=indexes)
    return(letter_alignment)

## dataframe with alignment as positions
def pos_dataframe(letter_alignment, column_names, indexes):
    position_alignment=pd.DataFrame(index=indexes[2:])
    for protein in column_names:
        pos_seq=[]
        start=letter_alignment.loc['start1',protein]
        seq=letter_alignment[protein].tolist()[2:]
        counter=0
        for aa in seq:
            if aa=="-":
                pos_seq.append(0)
            else:
                pos_seq.append(start+counter)
                counter+=1
        position_alignment[protein]=pos_seq
    return(position_alignment)

###### 0/1 dataframe with mapped phosphorylations onto alignment
def phos_dataframe(phosp_info, letter_alignment, position_alignment, column_names, indexes):
    phosp_alignment=pd.DataFrame(0,index=indexes[2:], columns=column_names)
    for column in column_names:
        protein_id=column.split("+")[0].strip()
        start_no=int(column.split("+")[1].strip())
        for i in range(0,len(phosp_info)):
            l=phosp_info[i].split(",")
            protein_id2=l[1].strip()
            phosp_pos=int(l[2].strip())
            phosp_aa=l[3].strip()
            if protein_id==protein_id2 and start_no<=phosp_pos:
                row_with_pos=position_alignment.index[position_alignment[column]==phosp_pos].tolist()
                if len(row_with_pos)==1:
                    phosp_alignment.at[row_with_pos[0],column]+=1
    return(phosp_alignment)

#######rolling window for a list (returns list [2:-2] of original length)
def count_window_for_list(list):                                                 
    window_values=[]                                                             
    for i in range(2,len(list)-2):                                               
        a=int(list[i-2])                                                         
        b=int(list[i-1])                                                         
        c=int(list[i])                                                           
        d=int(list[i+1])                                                         
        e=int(list[i+2])                                                         
        how_many=a+b+c+d+e                                                       
        bg=float(how_many/5.)                                                    
        window_values.append(bg)                                                 
    return(window_values)                                                        

## background construction
## make_random_histogram(phosp,total_S) of phosphosites across all available S (or T or Y or whatever)
def make_random_histogram(elements, size):
    rlist = [0 for i in range(size)]
    for i in [random.randint(0, size-1) for i in range(elements)]:
        rlist[i] += 1
    return(rlist)

## returns one permutated column of rolled window on background
def permutated_dataframe(letter_alignment,phosp_alignment, column_names, indexes):
    permutated_alignment=pd.DataFrame(0, index=indexes[2:], columns=column_names)
    for protein_id in column_names:
        rows_with_S=letter_alignment.index[letter_alignment[protein_id]=="S"].tolist()
        rows_with_T=letter_alignment.index[letter_alignment[protein_id]=="T"].tolist()
        rows_with_Y=letter_alignment.index[letter_alignment[protein_id]=="Y"].tolist()
        phosp_rows=phosp_alignment.index[phosp_alignment[protein_id]>=1].tolist()
        pS=pT=pY=0
        for i in phosp_rows:
                if i in rows_with_S:
                    pS+=1
                if i in rows_with_T:
                    pT+=1
                if i in rows_with_Y:
                    pY+=1
        random_S=make_random_histogram(pS, len(rows_with_S))
        random_T=make_random_histogram(pT, len(rows_with_T))
        random_Y=make_random_histogram(pY, len(rows_with_Y))
        for i,j in zip(random_S, rows_with_S):
            permutated_alignment.at[j,protein_id]=i
        for i,j in zip(random_T, rows_with_T):
            permutated_alignment.at[j,protein_id]=i
        for i,j in zip(random_Y, rows_with_Y):
            permutated_alignment.at[j,protein_id]=i
    sum_of_phosps=permutated_alignment.loc[0:].sum(axis=1).tolist()
    one_of_bg=count_window_for_list(sum_of_phosps)
    return(one_of_bg)

#permutated_alignment=permutated_dataframe(letter_alignment,phosp_alignment, column_names, indexes)
def count_zscore(fg,mean,stdev):
    if stdev==0:
        zscore=0
    else:
        zscore=(fg-mean)/stdev
    return(zscore)

def count_pval(fg, mean, stdev):
    fg=float(fg)
    mean=float(mean)
    stdev=float(stdev)
    if fg>mean:
        z=count_zscore(fg, mean, stdev)
        p_val = scipy.stats.norm.sf(abs(z))
    if fg<=mean:
        z=1.
        p_val=1.
    return(p_val)

def getHotspotSitesInFile(alignment_path, ptms, how_many_permuts):
    ## Starting dataframes
    (column_names, indexes)=prepare_cols_indx(alignment_path)
    letter_alignment=letter_ali_dataframe(alignment_path, column_names, indexes)
    position_alignment=pos_dataframe(letter_alignment, column_names, indexes)
    phosp_alignment=phos_dataframe(ptms, letter_alignment, position_alignment, column_names, indexes)
    ## Foreground dataframe
    sum_of_phosps=phosp_alignment.loc[0:].sum(axis=1).tolist()
    foreground=count_window_for_list(sum_of_phosps)
    ## dataframe (and permutations)
    background_dataframe=pd.DataFrame(index=indexes[4:-2],columns=list(range(0,how_many_permuts)))
    for k in range(0,how_many_permuts):
        background_dataframe[k]=permutated_dataframe(letter_alignment, phosp_alignment, column_names, indexes)
        background_dataframe['domain'] = os.path.splitext(os.path.basename(alignment_path))[0]
        background_dataframe['position_aln'] = indexes[4:-2]
        background_dataframe['bg_medians'] = background_dataframe.iloc[:,0:(how_many_permuts-1)].median(axis=1).tolist()
        background_dataframe['bg_stdev'] = background_dataframe.iloc[:,0:(how_many_permuts-1)].std(axis=1).tolist()
        background_dataframe['bg_means'] = background_dataframe.iloc[:,0:(how_many_permuts-1)].mean(axis=1).tolist()
        background_dataframe['foreground']=foreground
        ## this is to count pvals and add them to the column in background_dataframe
    all_pvals=[]
    for row in background_dataframe.index:
        fg=background_dataframe.at[row,'foreground']
        mean=background_dataframe.at[row,'bg_means']
        st_dev=background_dataframe.at[row,'bg_stdev']
        p_val=count_pval(fg, mean, st_dev)
        all_pvals.append(p_val)
        background_dataframe['pvals']=all_pvals
        background_dataframe['pvals']=background_dataframe.loc[:,('pvals')]+0.00000000000000001 # pseudocount for 0 values
        # background_dataframe = background_dataframe[["domain", "position_aln", "foreground", "pvals"]]
    background_dataframe = background_dataframe.loc[:,("domain", "position_aln", "foreground", "pvals")].copy()
    # adding protein/position level information
    letter_alignment_tomerge = letter_alignment[2:].copy()
    letter_alignment_tomerge["position_aln"] = letter_alignment_tomerge.index
    position_alignment["position_aln"] = position_alignment.index
    tomerge = pd.merge(letter_alignment_tomerge.infer_objects().melt(id_vars=["position_aln"],
                                                                     var_name = "protein",
                                                                     value_name = "residue"),
                       position_alignment.melt(id_vars=["position_aln"],
                                               var_name = "protein",
                                               value_name = "position"),
                       on = ["position_aln", "protein"])
    tomerge["protein"] = tomerge["protein"].str.split("+", expand = True)[0]
    tomerge["protein"] = tomerge["protein"].str.replace("#human", "")
    expanded_dataframe = background_dataframe.merge(tomerge, on = "position_aln")
    # ## fix position in alignment
    expanded_dataframe["position_aln"] = expanded_dataframe[["position_aln"]] + 1
    return(expanded_dataframe)

def get_hotspot_regions(hotspot_sites):
    ##Creting hotspot regions (+/-2)
    ## Positions + neighboring positions and whether they are hotspots or not
    df = hotspot_sites.loc[:,("domain", "protein","position", "position_aln")].drop_duplicates().copy()
    df["minus2"] = df.position_aln - 2
    df["minus1"] = df.position_aln - 1
    df["exact"] = df.position_aln
    df["plus1"] = df.position_aln + 1
    df["plus2"] = df.position_aln + 2
    df = pd.melt(df,
                 id_vars = ["domain", "protein", "position", "position_aln"],
                 var_name = "rel_position",
                 value_name = "neigh_position")
    df = df.drop(columns = ["rel_position"])
    ##hotspot-domain definitions
    hotspot_domains = pd.merge(df.loc[:,("domain", "position_aln", "neigh_position")].copy().drop_duplicates(),
                               hotspot_sites.loc[:,("domain", "position_aln", "hotspot")].drop_duplicates().rename(columns = {"position_aln":"neigh_position"}),
                               on = ["domain", "neigh_position"])
    hotspot_domains = hotspot_domains.groupby(["domain", "position_aln"])[["hotspot"]].any().reset_index()
    hotspot_domains = hotspot_domains.rename(columns = {"hotspot":"in_hotspot_region"})
    hotspot_domains["in_hotspot_region_inv"] = np.invert(hotspot_domains.in_hotspot_region)
    hotspot_domains = hotspot_domains.sort_values(by = ["domain", "position_aln"])
    hotspot_domains["in_hotspot_region_cumsum"] = hotspot_domains.groupby(["domain"])[["in_hotspot_region_inv"]].cumsum()
    hotspot_definitions = hotspot_domains.groupby(["domain"])["in_hotspot_region_cumsum"].value_counts().rename("counts").reset_index()
    hotspot_definitions = hotspot_definitions [hotspot_definitions.counts > 1].drop(columns = "counts")
    hotspot_definitions = pd.merge(hotspot_domains,
                                   hotspot_definitions,
                                   on = ["domain", "in_hotspot_region_cumsum"])
    hotspot_definitions = hotspot_definitions[hotspot_definitions.in_hotspot_region]
    hotspot_definitions["start_aln"] = hotspot_definitions.groupby(["domain", "in_hotspot_region_cumsum"])["position_aln"].transform("min")
    hotspot_definitions["end_aln"] = hotspot_definitions.groupby(["domain", "in_hotspot_region_cumsum"])["position_aln"].transform("max")
    hotspot_definitions["hotspot_id"] = hotspot_definitions["domain"].str.cat(hotspot_definitions.start_aln.astype(str), sep = "_")
    hotspot_definitions = hotspot_definitions.loc[:,("domain", "position_aln", "start_aln", "end_aln", "hotspot_id")]
    return(hotspot_definitions)

def f(x):
    d = {}
    d['sequence'] = x["residue"].str.cat()
    d['start'] = x[x.residue != "-"]["position"].min()
    d['end'] = x[x.residue != "-"]["position"].max()
    d['min_adj_pval'] = x['p_adjust'].min()
    d['foreground_max'] = x['foreground'].max()
    return pd.Series(d, index=["start", "end", "sequence", "foreground_max", 'min_adj_pval'])

def find_hotspot_instances(hospot_sites, hotspot_definitions):
    ## instances of the hotspots in the proteins
    hotspot_instances = pd.merge(hotspot_sites,
                                 hotspot_definitions,
                                 on = ["domain", "position_aln"],
                                 how = "left")
    hotspot_instances["hotspot_region"] = np.invert(hotspot_instances[["start_aln"]].isna())
    hotspot_instances["hotspot_region_inv"] = hotspot_instances[["start_aln"]].isna()
    hotspot_instances = hotspot_instances.loc[hotspot_instances['residue'] != "-"]
    hotspot_instances = hotspot_instances.sort_values(by = ["domain", "hotspot_id","protein", "position"])
    final_output = hotspot_instances.groupby(["domain", "protein", "hotspot_id", "start_aln", "end_aln"]).apply(f).reset_index()
    final_output["start"] = final_output["start"].astype("int")
    final_output["end"] = final_output["end"].astype("int")
    final_output["start_aln"] = final_output["start_aln"].astype("int")
    final_output["end_aln"] = final_output["end_aln"].astype("int")
    final_output["foreground_max"] = final_output["foreground_max"].astype("int")
    return(final_output)

###########
## MAIN
###########
if __name__ == '__main__':
    ## Reading arguments
    parser = argparse.ArgumentParser(description='Estimate PTM hotspots in sequence alignments')
    parser.add_argument('--dir',
                        nargs = '?',
                        default = "db/alignments",
                        action="store",
                        metavar="PATH",
                        dest = "alignments_dir",
                        help='directory containing fasta alignments (default: db/alignments)')
    parser.add_argument('--ptmfile',
                        nargs = '?',
                        default = "db/all_phosps",
                        action="store",
                        metavar="PATH",
                        dest = "ptm_file",
                        help='file containing PTMs (default: db/all_phosps)')
    parser.add_argument('-d',
                        '--domain',
                        nargs = '?',
                        action="store",
                        metavar="PFXXXXX",
                        dest = "domain",
                        help='predictions in domain (i.e. Protein kinase domain PF00069)')
    parser.add_argument('--iter',
                        nargs = '?',
                        default = 100,
                        action="store",
                        dest="how_many_permuts",
                        metavar="INTEGER",
                        type=int,
                        help='number of permutations (default: 100)')
    parser.add_argument('--threshold',
                        nargs = '?',
                        default = 0.05,
                        action="store",
                        metavar = "FLOAT",
                        dest="threshold",
                        type=float,
                        help='Bonferroni-corrected p-value threshold for calling hotspots (default: 0.05)')
    parser.add_argument('--foreground',
                        nargs = '?',
                        default = 2,
                        action="store",
                        metavar = "FLOAT",
                        dest="fore_val",
                        type=float,
                        help='effect-size cutoff for calling hotspots (default: 2)')
    parser.add_argument('-o',
                        '--out',
                        required=True,
                        action="store",
                        metavar = "PATH",
                        dest="outputFile",
                        help='output csv file')
    parser.add_argument('--printSites',
                        dest = "printSites",
                        action="store_true",
                        help='print all site predictions instead of hotspot regions',
                        default=False)

    results = parser.parse_args()
    alignments_dir = results.alignments_dir
    ptm_file = results.ptm_file
    how_many_permuts = results.how_many_permuts
    threshold = results.threshold
    fore_val = results.fore_val
    outputFile = results.outputFile
    printSites = results.printSites

    ## reading all PTMs
    ptms=open(ptm_file,"r").readlines()
    ## check if domain is specified. read all otherwise
    if results.domain:
        alignmentFiles = [results.domain + ".fasta"]
        if not os.path.isfile(os.path.join(alignments_dir, alignmentFiles[0])):
            sys.stderr.write("Domain file does not exist!")
            sys.exit(1)
            ## read all
    else:
        alignmentFiles = os.listdir(alignments_dir)

    ## estimating all hotspot sites
    allHotspots = []
    for filename in alignmentFiles:
        print("* " + filename)
        alignment_path = os.path.join(alignments_dir, filename)
        getHotspotSitesInFile(alignment_path, ptms, how_many_permuts)
        allHotspots.append(expanded_dataframe)
    hotspot_sites = pd.concat(allHotspots)
    hotspot_sites = hotspot_sites.loc[:,("domain", "protein", "position", "residue", "position_aln", "foreground", "pvals")].copy()
    hotspot_sites["p_adjust"] = sm.stats.multipletests(hotspot_sites.pvals.tolist(), method='fdr_bh')[1]
    hotspot_sites["hotspot"] = (hotspot_sites.p_adjust <= threshold) & (hotspot_sites.foreground >= fore_val)

    ## estimating hotspot ranges
    #All hotspot range-definitions 
    hotspot_definitions = get_hotspot_regions(hotspot_sites)
    #Match of all hotspots into   
    hotspot_regions = find_hotspot_instances(hotspot_sites, hotspot_definitions)

    print("Done!")

    if printSites:
        hotspot_sites[hotspot_sites.hotspot].to_csv(outputFile, index = False)
    else:
        hotspot_regions.to_csv(outputFile, index = False)
