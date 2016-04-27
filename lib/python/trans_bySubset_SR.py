#!/usr/bin/env python

from __future__ import division
#from csv   import reader
import pandas as pd
import numpy as np
from collections import defaultdict
from math import log
#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from sys import stderr, argv, exit
from json import load as json_load
from argparse import ArgumentParser
# import matplotlib.pyplot as plt
# import numpy as npe

#############################
## Functions
##########    
# def get_transitions(in_df, dict_all_trans, pseudo=True, mode="dosage"):
def get_transitions(in_df, set_all_drugs, pseudo=True, mode="dosage"):
    '''
    Get transitions by patient inside a pandas dataframe
    '''
    
    ## Variables
    # keep adding dosage into a dictionary until I reach t2
    # dict_trans = defaultdict(int)
    index_drug_set = set()
    dict_trans = defaultdict(lambda: defaultdict(int))
    
    # dict_trans_counts = defaultdict(lambda: defaultdict(int))
    
    # Index to get all dosage of a drug for calculating the odds 
    # of index drug not having transition to marker drug
    dict_index_dosage = defaultdict(int)
    
    for name, group in in_df.groupby('id'):
        
        index_drug_set_by_id = set()
    
        ## We reset the index of the rows otherwise the original is conserved
        ## Everything at the same step below
        # group = group.reset_index(drop=True)
        ## We sort by month important for reversed table and reset indexes
        group = group.sort_values(ascending=True, by = ['month']).reset_index(drop=True)

        # print ("table for group is:", group) #del
    
        for i, row in group.iterrows():
            last_row_index = group.shape[0]
        
            index_drug = row['drug']
        
            # We just take into account the first appearance of index drug
            if index_drug in index_drug_set_by_id:
                # print ("I am coming here after first appearance", index_drug)
                continue
        
            # Only first time that index drug appears
            # Get all dosage of a index drug to use in cases where there is any transitions
            # to a given transition drug
            if mode == "dosage":
                dosage_index_drug_all = group.loc[group['drug'] == index_drug, 'dosage'].astype(int).sum()
                # if index_drug == 'A02BC01': print "====================", dosage_index_drug_all
                dict_index_dosage[index_drug] += dosage_index_drug_all
            elif mode == "counts":
                counts_index_drug_all = len(group.loc[group['drug'] == index_drug, 'dosage'])
                dict_index_dosage[index_drug] += counts_index_drug_all
        
            t1 = row['month']
            
            index_drug_set.add(index_drug)
            index_drug_set_by_id.add(index_drug)
            
            month_index_drug = row['month']
        
            group_cp = group.copy(deep=True)
            
            # This dictionary keeps the appearance of marker drugs so that if it appears
            # again it is not considered twice
            i_m_drug_trans_set_by_id = set()
            
            for j in range(i+1, last_row_index):
            
                if month_index_drug == group_cp.iloc[j]['month']:
                    continue
                
                marker_drug = group_cp.iloc[j]['drug']
                trans_i_m = (index_drug, marker_drug)
                if index_drug == marker_drug or trans_i_m in i_m_drug_trans_set_by_id:
                    continue
                
                i_m_drug_trans_set_by_id.add (trans_i_m)

                # Now i have to find last appearence of marker drug => t2
                # index of all appearance of marker drug
                idx = group_cp[group_cp['drug']==marker_drug].index.tolist()
            
                t2 = group_cp.iloc[idx[-1]]['month']

                # print ("index drug, time first index, time last marker", index_drug, marker_drug, t1, t2)
                mask_index_drug = (group['month'] >= t1) & (group['month'] <= t2) & (group['drug'] == index_drug)
                mask_marker_drug = (group['month'] > t1) & (group['month'] <= t2) & (group['drug'] == marker_drug)
                
                if mode == "dosage":    
                    dosage_index_drug = group.loc[mask_index_drug, 'dosage'].astype(int).sum()
                    dosage_marker_drug = group.loc[mask_marker_drug, 'dosage'].astype(int).sum()
                    
                    dict_trans[(index_drug, marker_drug)][index_drug] += dosage_index_drug 
                    dict_trans[(index_drug, marker_drug)][marker_drug] += dosage_marker_drug
                    # print "dosage sum index marker", index_drug, dosage_index_drug, marker_drug, dosage_marker_drug
                    # print "dict_trans", dict_trans
                    # print "dosage sum index", index_drug, dosage_index_drug
                elif mode == "counts":
                    counts_index_drug = len(group.loc[mask_index_drug, 'dosage'])
                    counts_marker_drug = len(group.loc[mask_marker_drug, 'dosage'])
                    # Naming differently dict_trans and dict_trans_counts could lead to problems
                    # dict_trans_counts[(index_drug, marker_drug)][index_drug] += counts_index_drug 
                    # dict_trans_counts[(index_drug, marker_drug)][marker_drug] += counts_marker_drug
                    dict_trans[(index_drug, marker_drug)][index_drug] += counts_index_drug 
                    dict_trans[(index_drug, marker_drug)][marker_drug] += counts_marker_drug
                # REVIEW
                # Not compatible bw them, I have to declare the dictionary again
                # tengo que sumar los divisiones o sumar las dosis primero y luego dividirlo todo
                # dict_trans[(index_drug, marker_drug)] += dosage_index_drug/dosage_marker_drug
            
                # The other way is how it is done

    
    if pseudo:
        dict_trans = add_pseudo_acc_dict (set_all_drugs = set_all_drugs, index_drug_set = index_drug_set, d_trans_subset = dict_trans, d_index_dosage = dict_index_dosage)
    
    print ("End of function")
    return(dict_trans, dict_index_dosage, index_drug_set)

####################
# dict_all_trans contains in principle all the drugs
# becuase it comes from the orignal table
# def add_pseudo_acc_dict (dict_all_trans, index_drug_set, d_trans_subset, d_index_dosage):
def add_pseudo_acc_dict (set_all_drugs, index_drug_set, d_trans_subset, d_index_dosage):
    # for i in sorted(dict_all_trans.keys()):
    for i in sorted(set_all_drugs):    
        #for j in list_all_drugs:
        # for j in sorted(dict_all_trans[i].keys()):
        # for j in sorted(dict_all_trans.keys()):
        for j in sorted(set_all_drugs):    
            # We need pseudo for transition between same 
            # if i == j:
            #     continue
            ## If we do not see the index drug then the transition is never found
            if i not in index_drug_set:
                d_trans_subset[(i, j)][i] += 1
                d_trans_subset[(i, j)][j] += 1
                continue
            
            # If the index is present but not the transition to the marker
            elif not d_trans_subset.get((i,j), False):
                # Changed because otherwise this guys are always showing up
#                 d_trans_subset[(i, j)][i] += d_index_dosage[i] + 1
                d_trans_subset[(i, j)][i] += 1
                d_trans_subset[(i, j)][j] += 1
            
            else:
                d_trans_subset[(i, j)][i] += 1
                d_trans_subset[(i, j)][j] += 1
                
    return (d_trans_subset)   

##########
# def SR_dosage_dict (dict_all_trans, dict_trans_drug, dict_trans_drug_rev, mode="ratio_avg"):
def SR_dosage_dict (set_all_drugs, dict_trans_drug, dict_trans_drug_rev, mode="ratio_avg"):    
    d_trans_for = dict_trans_drug   
    d_trans_rev = dict_trans_drug_rev
    d_sr = defaultdict(lambda: defaultdict(int))
    #drug_set = list_drugs
    #dict_all_trans
    
    # # Comparing set of drugs between forward and reverse
    # if (drug_set != drug_set_rev):
    #     raise ValueError("Drug sets contain different elements!!!")
    
    sr = list()
    sr_trans = list()
    
    # for i in sorted(dict_all_trans.keys()):
    for i in sorted(set_all_drugs):
        # for j in sorted (drug_set):
        # for j in sorted(dict_all_trans[i].keys()):
        for j in sorted(set_all_drugs):
            # print( "------", i,j)
            if i == j:
                sr.append(0)
                d_sr[i][j] = 0
                sr_trans.append((i,j))
            else:
                #dict_tr[(i, j)][i] / dict_tr[(i, j)][j]
                # sr.append (d_trans_for[(i,j)] / d_trans_rev[(i,j)])
                # is faster to append to a list and then transform to numpy 
                # than append in numpy array
                # sr_val = log(dict_tr[(i, j)][i] /  d_trans_rev[(i, j)][j])
                #review calculation
                marker_for = dict_trans_drug[(i, j)][j]
                index_for = dict_trans_drug[(i, j)][i]
                marker_rev = dict_trans_drug_rev[(i, j)][j]
                index_rev = dict_trans_drug_rev[(i, j)][i]
                # sr_val = log((dict_tr[(i, j)][j]/dict_tr[(i, j)][i]) / (d_trans_rev[(i, j)][j]/d_trans_rev[(i, j)][i]))
                # if (marker_rev/index_rev) == 0 : #del
                #     print >> stderr, "index ", i
                #     print >> stderr, "marker ", j
                sr_val = log((marker_for/index_for) / (marker_rev/index_rev))
                # sr.append (log(dict_tr[(i, j)][i] /  d_trans_rev[(i, j)][j]))
                sr.append(sr_val)
                d_sr[i][j] = sr_val
                sr_trans.append((i,j))
            
    # return (sr, sr_trans, dict_all_trans)
    return (sr, sr_trans, d_sr)

def plot_heatmap (ary_data, lab_cols="", lab_rows="", colors=plt.cm.jet, path_fig='/Users/jespinosa/2015_viscMes/foo.png'):
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(ary_data, cmap=colors)
    
    #legend
    cbar = plt.colorbar(heatmap)
    #cbar.ax.set_yticklabels(['0','1','2','>3'])
    cbar.set_label('SR', rotation=270)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(ary_data.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(ary_data.shape[0]) + 0.5, minor=False)
    ax.invert_yaxis()
    # plt.xticks(xdata,catnames,rotation=90)
    plt.xticks(rotation=90)
    #labels
    column_labels = lab_cols
    row_labels = lab_rows
    ax.set_xticklabels(lab_cols, minor=False, size=0.2)
    ax.set_yticklabels(lab_rows, minor=False, size=0.2)

    #plt.show()
    #plt.figure(figsize=(15,15), dpi=300)
    plt.savefig(path_fig, figsize=(15,15), dpi=300)    

###################################
### code
###################################
##$HOME/git/wasser/lib/python/trans_bySubset_SR.py -i $HOME/2015_viscMes/data/SR_analysis/traj_annotated.csv -t $HOME/2015_viscMes/data/SR_analysis/dict_traj_found.json -s 'm' -a '30 35' -m 'dosage'
####################################

parser = ArgumentParser(description='File input arguments')
parser.add_argument('-i','--input', help='Input Visc + file name', required=True)
parser.add_argument('-t','--transitions', help='Transitions dictionary file', required=True)
parser.add_argument('-s','--sex', help='Sex to filter from table', required=True)
parser.add_argument('-a','--age', help='Age to filter from data', required=True)
parser.add_argument('-m','--mode', help='Mode to calculate SR, counts or dosage', required=True)

args = parser.parse_args()

print ("Input file: %s" % args.input)
print ("Transition file: %s" % args.transitions)
print ("Sex option: %s" % args.sex)
print ("Age option: %s" % args.age)
print ("Age option: %s" % args.mode)

# Input files
dosage_traj_file = args.input
transitions_file = args.transitions

##########################
## Options to subset data
sex_opt = list(args.sex.split(" "))
age_opt_str = list(args.age.split(" "))
age_opt = map(int, age_opt_str)
mode_opt = args.mode

print ("sex option is:", sex_opt)
# print ("age option is:", age_opt)
print ("type of sex option is:", type(sex_opt))
print ("mode option is:", mode_opt)

# print ("type of age option is:", type(age_opt))
print >> stderr, "File trajectories annotation ", dosage_traj_file
print >> stderr, "File trajectories annotation",  type(dosage_traj_file)
print >> stderr, "File trajectories annotation",  dosage_traj_file
print >> stderr, "File transition dictionary", transitions_file

#/Library/Frameworks/Python.framework/Versions/2.7/envs/pytables/bin/python
# df_dosage_traj_groups = pd.read_csv('/Users/jespinosa/2015_viscMes/data/SR_analysis/traj_annotated.csv', index_col=0)#comment
df_dosage_traj_groups = pd.read_csv(dosage_traj_file, index_col=0)

## reading dictionary with all the transitions found in the original data
# It is better just to recover a set with all possible drugs 
set_drugs = set(df_dosage_traj_groups['drug'])



# with open('/Users/jespinosa/2015_viscMes/data/SR_analysis/dict_traj_found.json', 'r') as f: #comment
    # trans_dict = json_load(f)
# with open(transitions_file, 'r') as f:
#         trans_dict = json_load(f)
        
# print >> stderr, ("Transitions dictionary contains %i keys" %  len(trans_dict.keys()))

df_subset_for = df_dosage_traj_groups.loc[df_dosage_traj_groups['sex'].isin(sex_opt) & df_dosage_traj_groups['ages'].isin(age_opt)]

### First reverse the data
# Reverse dataframe
df_subset_rev = df_subset_for.copy()
df_subset_rev = df_subset_rev.iloc[::-1]
max_t = int(max(df_subset_rev["month"]))
df_subset_rev["month"] = (df_subset_rev["month"].astype(int) - max_t).abs()

# This is insde the function because I just want to change the order for each id (pacient)
# df_subset_rev = df_subset_rev.sort_values(ascending=True, by = ['month']).reset_index(drop=True)

(dict_tr, dict_index_dosage, all_index_drug) = get_transitions(in_df=df_subset_for,  set_all_drugs=set_drugs, pseudo=True, mode=mode_opt)
(dict_tr_rev, dict_index_dosage_rev, all_index_drug_rev) = get_transitions(in_df=df_subset_rev,  set_all_drugs=set_drugs, pseudo=True, mode=mode_opt)

#(sr_log, list_trans, list_drug_returned) = SR_dosage(list_drugs=list_all_drug, dict_trans_drug = dict_tr, dict_trans_drug_rev= dict_tr_rev)
# (sr_log, list_trans, dict_sr_log) = SR_dosage_dict(dict_all_trans = trans_dict, dict_trans_drug = dict_tr, dict_trans_drug_rev= dict_tr_rev, mode="ratio_avg")
(sr_log, list_trans, dict_sr_log) = SR_dosage_dict(set_all_drugs=set_drugs, dict_trans_drug = dict_tr, dict_trans_drug_rev= dict_tr_rev, mode="ratio_avg")


ext_file = ".csv"
sr_file_name = "sr_" + "_".join(age_opt_str) + "_" + "_".join(sex_opt) + "_" + mode_opt
# print >> stderr, (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", sr_file_name)
# print >> stderr, (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", sr_file_name + ext_file)
print >> stderr, (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", sr_file_name + "_tbl" + ext_file)
print >> stderr, (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", sr_file_name + "_mat" + ext_file)

set_index, set_marker = map(set,zip(*list_trans))

path_tbl = sr_file_name + "_tbl" + ext_file
path_mat = sr_file_name + "_mat" + ext_file

sr_file_mat = open(path_mat, "w")
sr_file_tbl = open(path_tbl, "w")

sr_file_mat.write("index_drug\t") # let column with index below this label
# sr_file_mat.write('\t'.join(item for item in sorted (drug_set)))
sr_file_mat.write('\t'.join(item for item in sorted (set_marker)))
sr_file_mat.write("\n")

# All missing transitions can be inputed as 0 because the log will be 0, 0

sr_list_all_trans = list()

# ('N06AB03', 'N06AB03')
# for n, i in enumerate(sorted(drugs_set_for)):
# for index_dr in sorted(set_index):
for index_dr in sorted(set_drugs):
    sr_file_mat.write('%s' % (index_dr))

    # for j in sorted (drugs_set_for):
    # for marker_dr in sorted(set_marker):
    for marker_dr in sorted(set_drugs):
        #sr_file_mat.write('%s' % (marker_dr))
        sr_file_mat.write('\t%f' % (dict_sr_log[index_dr][marker_dr]))
        # if dict_sr_log[index_dr][marker_dr] != 0 : 
            # print >> stderr, (">>>>>>>>>>>>>>>>>>>>>dict>>>>>>>>>>>>>>>", index_dr, marker_dr, dict_sr_log[index_dr][marker_dr])
        # sr_file_mat.write('\t%f' % (sr_v[n]))
        sr_list_all_trans.append(dict_sr_log[index_dr][marker_dr])
        sr_file_tbl.write('%s\t%s\t%f\n' % (index_dr, marker_dr, dict_sr_log[index_dr][marker_dr]))

    sr_file_mat.write('\n')    

sr_file_tbl.close()
sr_file_mat.close()

n = len(set_index)
m = len(set_marker)
len_sr_log = len(sr_log)
len_n_m = n*m
len_all_trans = len(sr_list_all_trans)

#print >> stderr, (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", len_sr_log, len_n_m, len_all_trans)

#dict_tr_test = defaultdict(lambda: defaultdict(int)) #del test
#print dict_tr_test['a']['b'] #del test
# exit("This is execution has been fucking good")

sr_ary = np.array(sr_list_all_trans).reshape(n, m)
labels_heatmap_col = list(sorted(set_index))
labels_heatmap_row = list(sorted(set_marker))
#name_fig = k_gr + type_img
type_img = '.pdf'
name_fig = sr_file_name + type_img
#sr_file_name = "sr_" + "_".join(age_opt_str) + "_" + "_".join(sex_opt) + "_" + mode_opt

#plot_heatmap (sr_ary, labels_heatmap, colors=plt.cm.Blues, path_fig=path_fig)
plot_heatmap (sr_ary, lab_cols=labels_heatmap_col, lab_rows=labels_heatmap_row, colors=plt.cm.jet, path_fig=name_fig)

# sr_file_mat.write('\t'.join(item for item in sorted (drug_set)))
# sr_file_mat.write("\n")
print >> stderr, (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", len_sr_log, len_n_m, len_all_trans)
