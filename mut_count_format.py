#this code is designed to take a mutation-counts file from shapemapper for both the untreated and treated samples to:
##1 calculate mutation rates for each mutation type
##2 calculate raw reactivities for each mutation type (modified mutation rate - untreated mutation rate)
##3 output one file



import os
import sys
import argparse
import pandas as pd
import numpy as np

#set arguments
parser = argparse.ArgumentParser(description="This script takes mutation counts fiiles from shapemapper and outputs mutation rates and raw reactivities")
parser.add_argument("-ut","--ut", help="The unreated input file")
parser.add_argument("-t","--t", help="The treated input file")
parser.add_argument("output", help="The output file identifier")
parser.add_argument("-muttype", help="The mutation type csv file")
parser.add_argument("-primer", help="length of RT primer, default is 10", default=10, type=int)
parser.add_argument("--raw-react", help="output raw reactivity file", action="store_true")
parser.add_argument("--raw-react-filt", help="output raw reactivity file with filtering", action = "store_true")
parser.add_argument("--nind", help="normalize reactivity individually", action="store_true")
args = parser.parse_args()


############## step 1 ##############
#read in the mutation counts
def read_input(input_file):
    df = pd.read_csv(input_file, sep="\t", header=0)
    return df

df_ut = read_input(args.ut)
df_t = read_input(args.t)


############## step 2 ##############

#combine columns (ie A-, T-, G-, and C- to 1nt_del)
def combine_cols(df):
    df1 = df.copy() #make a copy to avoid modifying the original DataFrame
    #combine columns
    df1["1nt_del"] = df1["A-"] + df1["T-"] + df1["G-"] + df1["C-"]
    df1["1nt_ins"] = df1["-A"] + df1["-T"] + df1["-G"] + df1["-C"] + df1["-N"]
    df1["multinuc_mismatch"]= df1["A_multinuc_mismatch"] + df1["T_multinuc_mismatch"] + df1["G_multinuc_mismatch"] + df1["C_multinuc_mismatch"] + df1["N_multinuc_mismatch"]
    #drop the columns that were combined
    df1 = df1.drop(columns=["A-","T-","G-","C-","-A","-T","-G","-C","-N","A_multinuc_mismatch","T_multinuc_mismatch","G_multinuc_mismatch","C_multinuc_mismatch","N_multinuc_mismatch"])
    #drop extra depth columns
    df1 = df1.drop(columns=["read_depth","off_target_mapped_depth","low_mapq_mapped_depth","mapped_depth"])
    return df1

df_ut = combine_cols(df_ut)
df_t = combine_cols(df_t)

############## step 3 ##############
def calc_mutation_rate(df,sample):
    df1 = df.copy()  # Create a copy to avoid modifying the original DataFrame
    #calculate mutation rates
    df1[sample+"_1nt_del_rate"] = df1["1nt_del"] / df1["effective_depth"]
    df1[sample+"_1nt_ins_rate"] = df1["1nt_ins"] / df1["effective_depth"]
    df1[sample+"_AT_rate"] = df1["AT"] / df1["effective_depth"]
    df1[sample+"_AG_rate"] = df1["AG"] / df1["effective_depth"]
    df1[sample+"_AC_rate"] = df1["AC"] / df1["effective_depth"]
    df1[sample+"_TA_rate"] = df1["TA"] / df1["effective_depth"]
    df1[sample+"_TG_rate"] = df1["TG"] / df1["effective_depth"]
    df1[sample+"_TC_rate"] = df1["TC"] / df1["effective_depth"]
    df1[sample+"_GA_rate"] = df1["GA"] / df1["effective_depth"]
    df1[sample+"_GT_rate"] = df1["GT"] / df1["effective_depth"]
    df1[sample+"_GC_rate"] = df1["GC"] / df1["effective_depth"]
    df1[sample+"_CA_rate"] = df1["CA"] / df1["effective_depth"]
    df1[sample+"_CT_rate"] = df1["CT"] / df1["effective_depth"]
    df1[sample+"_CG_rate"] = df1["CG"] / df1["effective_depth"]
    df1[sample+"_multint_mismatch_rate"] = df1["multinuc_mismatch"] / df1["effective_depth"]
    df1[sample+"_multint_del_rate"] = df1["multinuc_deletion"] / df1["effective_depth"]
    df1[sample+"_multint_ins_rate"] = df1["multinuc_insertion"] / df1["effective_depth"]
    df1[sample+"_complex_del_rate"] = df1["complex_deletion"] / df1["effective_depth"]
    df1[sample+"_complex_ins_rate"] = df1["complex_insertion"] / df1["effective_depth"]
    df1[sample+"_effective_depth"] = df1["effective_depth"]
    #drop raw counts columns
    df1 = df1.drop(columns=["1nt_del","1nt_ins","AT","AG","AC","TA","TG","TC","GA","GT","GC","CA","CT","CG","multinuc_mismatch","multinuc_deletion","multinuc_insertion","complex_deletion","complex_insertion"])
    return df1

df_ut = calc_mutation_rate(df_ut,"UT")
df_t = calc_mutation_rate(df_t,"T")





############## step 4 ##############

#combine dataframes and calculate raw reactivity
def calc_raw_reactivity(ut_df,t_df):
    #combine dataframes
    df = pd.concat([ut_df, t_df], axis=1)
    #calculate raw reactivity for each mutation type, save as NA if UT_rate < 0.05
    df["1nt_del_raw_react"] = np.where(df["UT_1nt_del_rate"] < 0.05, df["T_1nt_del_rate"] - df["UT_1nt_del_rate"], np.nan)
    df["1nt_ins_raw_react"] = np.where(df["UT_1nt_ins_rate"] < 0.05, df["T_1nt_ins_rate"] - df["UT_1nt_ins_rate"], np.nan)
    df["AT_raw_react"] = np.where(df["UT_AT_rate"] < 0.05, df["T_AT_rate"] - df["UT_AT_rate"], np.nan)
    df["AG_raw_react"] = np.where(df["UT_AG_rate"] < 0.05, df["T_AG_rate"] - df["UT_AG_rate"], np.nan)
    df["AC_raw_react"] = np.where(df["UT_AC_rate"] < 0.05, df["T_AC_rate"] - df["UT_AC_rate"], np.nan)
    df["TA_raw_react"] = np.where(df["UT_TA_rate"] < 0.05, df["T_TA_rate"] - df["UT_TA_rate"], np.nan)
    df["TG_raw_react"] = np.where(df["UT_TG_rate"] < 0.05, df["T_TG_rate"] - df["UT_TG_rate"], np.nan)
    df["TC_raw_react"] = np.where(df["UT_TC_rate"] < 0.05, df["T_TC_rate"] - df["UT_TC_rate"], np.nan)
    df["GA_raw_react"] = np.where(df["UT_GA_rate"] < 0.05, df["T_GA_rate"] - df["UT_GA_rate"], np.nan)
    df["GT_raw_react"] = np.where(df["UT_GT_rate"] < 0.05, df["T_GT_rate"] - df["UT_GT_rate"], np.nan)
    df["GC_raw_react"] = np.where(df["UT_GC_rate"] < 0.05, df["T_GC_rate"] - df["UT_GC_rate"], np.nan)
    df["CA_raw_react"] = np.where(df["UT_CA_rate"] < 0.05, df["T_CA_rate"] - df["UT_CA_rate"], np.nan)
    df["CT_raw_react"] = np.where(df["UT_CT_rate"] < 0.05, df["T_CT_rate"] - df["UT_CT_rate"], np.nan)
    df["CG_raw_react"] = np.where(df["UT_CG_rate"] < 0.05, df["T_CG_rate"] - df["UT_CG_rate"], np.nan)
    df["multint_mismatch_raw_react"] = np.where(df["UT_multint_mismatch_rate"] < 0.05, df["T_multint_mismatch_rate"] - df["UT_multint_mismatch_rate"], np.nan)
    df["multint_del_raw_react"] = np.where(df["UT_multint_del_rate"] < 0.05, df["T_multint_del_rate"] - df["UT_multint_del_rate"], np.nan)
    df["multint_ins_raw_react"] = np.where(df["UT_multint_ins_rate"] < 0.05, df["T_multint_ins_rate"] - df["UT_multint_ins_rate"], np.nan)
    df["complex_del_raw_react"] = np.where(df["UT_complex_del_rate"] < 0.05, df["T_complex_del_rate"] - df["UT_complex_del_rate"], np.nan)
    df["complex_ins_raw_react"] = np.where(df["UT_complex_ins_rate"] < 0.05, df["T_complex_ins_rate"] - df["UT_complex_ins_rate"], np.nan)
    #return the dataframe
    return df

raw_react_df = calc_raw_reactivity(df_ut,df_t)







############## step 5 ##############

#add in sequence column
def add_seq(df):
    #add sequence column and set every row equal to N
    df["sequence"] = "N"
    #if AG_raw_react and CT_raw_react and GA_raw_react are equal to 0 and TC_raw_react is not zero, change sequence column to T
    df.loc[(df["AG_raw_react"] == 0) & (df["CT_raw_react"] == 0) & (df["GA_raw_react"] == 0) & (df["TC_raw_react"] != 0), "sequence"] = "T"
    #if AG_raw_react and CT_raw_react and TC_raw_react are equal to 0 and GA_raw_react is not zero, change sequence column to G
    df.loc[(df["AG_raw_react"] == 0) & (df["CT_raw_react"] == 0) & (df["TC_raw_react"] == 0) & (df["GA_raw_react"] != 0), "sequence"] = "G"
    #if AG_raw_react and TC_raw_react and GA_raw_react are equal to 0 and CT_raw_react is not zero, change sequence column to C
    df.loc[(df["AG_raw_react"] == 0) & (df["TC_raw_react"] == 0) & (df["GA_raw_react"] == 0) & (df["CT_raw_react"] != 0), "sequence"] = "C"
    #if TC_raw_react and CT_raw_react and TC_raw_react and GA_raw_react are equal to 0 and AT_raw_react is not zero, change sequence column to A
    df.loc[(df["TC_raw_react"] == 0) & (df["CT_raw_react"] == 0) & (df["GA_raw_react"] == 0) & (df["AT_raw_react"] != 0), "sequence"] = "A"
    #return the dataframe
    return df

raw_react_df = add_seq(raw_react_df)

#if raw-react flag is set, write the raw reactivity file
if args.raw_react:
    #get name of output file from the input file, save the text before the second underscore
    output_file = args.output + "_muttype_raw_react.txt"
    #output the profile file to the current directory
    raw_react_df.to_csv(output_file, sep="\t", index=False)
#write to csv






############## step 6 ##############

#remove rates that are filtered out
def remove_filtered(df,mutttype):
    #read in muttype csv as a dataframe
    muttype_df = pd.read_csv(mutttype, sep=",", header=0)
    #add "_raw_react" to the values in the Mutation column
    muttype_df["Mutation"] = muttype_df["Mutation"] + "_raw_react"
    #make a list of "mutations" that equal X in column A
    A_mut_list = muttype_df[muttype_df["A"] == "X"]["Mutation"].tolist()
    C_mut_list = muttype_df[muttype_df["C"] == "X"]["Mutation"].tolist()
    G_mut_list = muttype_df[muttype_df["G"] == "X"]["Mutation"].tolist()
    T_mut_list = muttype_df[muttype_df["U"] == "X"]["Mutation"].tolist()
    #if "sequence" in df == A, set values in columns in A_mut_list to NA
    df.loc[df["sequence"] == "A", A_mut_list] = np.nan
    df.loc[df["sequence"] == "C", C_mut_list] = np.nan
    df.loc[df["sequence"] == "G", G_mut_list] = np.nan
    df.loc[df["sequence"] == "T", T_mut_list] = np.nan
    return df

#if args.muttype is not empty, remove filtered rates
if args.muttype:
    raw_react_df = remove_filtered(raw_react_df,args.muttype)

#if raw-react flag is set, write the raw reactivity file
if args.raw_react_filt:
    #get name of output file from the input file, save the text before the second underscore
    output_file = args.output + "_muttype_raw_react_filt.txt"
    #output the profile file to the current directory
    raw_react_df.to_csv(output_file, sep="\t", index=False)
#write to csv

############## step 7 ##############

#normalize for reactivity
def norm_individual(df,primer):
    #add mutation rates for the untreated to a new column "UT_all"
    df["UT_all"] = df["UT_1nt_del_rate"] + df["UT_1nt_ins_rate"] + df["UT_AT_rate"] + df["UT_AG_rate"] + df["UT_AC_rate"] + df["UT_TA_rate"] + df["UT_TG_rate"] + df["UT_TC_rate"] + df["UT_GA_rate"] + df["UT_GT_rate"] + df["UT_GC_rate"] + df["UT_CA_rate"] + df["UT_CT_rate"] + df["UT_CG_rate"] + df["UT_multint_mismatch_rate"] + df["UT_multint_del_rate"] + df["UT_multint_ins_rate"] + df["UT_complex_del_rate"] + df["UT_complex_ins_rate"]
    #if UT_all > 0.05, set reactivity to NA
    df.loc[df["UT_all"] > 0.05, "reactivity"] = np.nan
    #set last #rows-primer : last row of reactivity to NA
    df.iloc[-primer:,-1] = np.nan
    #add column "nt" that counts from 1 to end of dataframe
    df["nt"] = range(1,len(df)+1)
    #add all raw_reactivities to a new column "all_raw_react"
    df["all_raw_react"] = df[["1nt_del_raw_react", "1nt_ins_raw_react", "AT_raw_react", "AG_raw_react", "AC_raw_react", "TA_raw_react", "TG_raw_react", "TC_raw_react", "GA_raw_react", "GT_raw_react", "GC_raw_react", "CA_raw_react", "CT_raw_react", "CG_raw_react", "multint_mismatch_raw_react", "multint_del_raw_react", "multint_ins_raw_react", "complex_del_raw_react", "complex_ins_raw_react"]].sum(axis=1, skipna=True)
    #save the 90th percentile of all_raw_react where sequence == A
    A_90, A_98 = df.loc[df["sequence"] == "A", "all_raw_react"].quantile([0.9, 0.98])
    #A div is the average of all values between the 90th and 98th percentile
    A_div = df.loc[(df["sequence"] == "A") & (df["all_raw_react"] >= A_90) & (df["all_raw_react"] <= A_98), "all_raw_react"].mean()
    #save the same information for other nucleotides 
    #save 90th percentile for all_raw_react where sequence == C
    C_90, C_98 = df.loc[df["sequence"] == "C", "all_raw_react"].quantile([0.9, 0.98])
    C_div = df.loc[(df["sequence"] == "C") & (df["all_raw_react"] >= C_90) & (df["all_raw_react"] <= C_98), "all_raw_react"].mean()
    G_90, G_98 = df.loc[df["sequence"] == "G", "all_raw_react"].quantile([0.9, 0.98])
    G_div = df.loc[(df["sequence"] == "G") & (df["all_raw_react"] >= G_90) & (df["all_raw_react"] <= G_98), "all_raw_react"].mean()
    T_90, T_98 = df.loc[df["sequence"] == "T", "all_raw_react"].quantile([0.9, 0.98])
    T_div = df.loc[(df["sequence"] == "T") & (df["all_raw_react"] >= T_90) & (df["all_raw_react"] <= T_98), "all_raw_react"].mean()
    #normalize the reactivity values to new column, "reactivity"
    df["reactivity"] = np.nan
    df.loc[df["sequence"] == "A", "reactivity"] = df["all_raw_react"] / A_div
    df.loc[df["sequence"] == "C", "reactivity"] = df["all_raw_react"] / C_div
    df.loc[df["sequence"] == "G", "reactivity"] = df["all_raw_react"] / G_div
    df.loc[df["sequence"] == "T", "reactivity"] = df["all_raw_react"] / T_div
    #make a final dataframe with only nt# counting from 1+, sequence, and reactivity
    final_df = df[["nt","sequence","reactivity"]]
    #replace NA with "NA"
    final_df = final_df.fillna("NA")
    #write final df to csv
    final_df.to_csv(args.output + "_reactivity.csv", sep=",", index=False)


def norm_together(df,primer):
    #add mutation rates for the untreated to a new column "UT_all"
    df["UT_all"] = df["UT_1nt_del_rate"] + df["UT_1nt_ins_rate"] + df["UT_AT_rate"] + df["UT_AG_rate"] + df["UT_AC_rate"] + df["UT_TA_rate"] + df["UT_TG_rate"] + df["UT_TC_rate"] + df["UT_GA_rate"] + df["UT_GT_rate"] + df["UT_GC_rate"] + df["UT_CA_rate"] + df["UT_CT_rate"] + df["UT_CG_rate"] + df["UT_multint_mismatch_rate"] + df["UT_multint_del_rate"] + df["UT_multint_ins_rate"] + df["UT_complex_del_rate"] + df["UT_complex_ins_rate"]
    #if UT_all > 0.05, set reactivity to NA
    df.loc[df["UT_all"] > 0.05, "reactivity"] = np.nan
    #set last #rows-primer : last row of reactivity to NA
    df.iloc[-primer:,-1] = np.nan
    #add column "nt" that counts from 1 to end of dataframe
    df["nt"] = range(1,len(df)+1)
    #add all raw_reactivities to a new column "all_raw_react"
    df["all_raw_react"] = df[["1nt_del_raw_react", "1nt_ins_raw_react", "AT_raw_react", "AG_raw_react", "AC_raw_react", "TA_raw_react", "TG_raw_react", "TC_raw_react", "GA_raw_react", "GT_raw_react", "GC_raw_react", "CA_raw_react", "CT_raw_react", "CG_raw_react", "multint_mismatch_raw_react", "multint_del_raw_react", "multint_ins_raw_react", "complex_del_raw_react", "complex_ins_raw_react"]].sum(axis=1, skipna=True)
    #save the 90th percentile of all_raw_react
    q_90, q_98 = df["all_raw_react"].quantile([0.9, 0.98])
    div = df.loc[(df["all_raw_react"] >= q_90) & (df["all_raw_react"] <= q_98), "all_raw_react"].mean()
    #normalize the reactivity values to new column, "reactivity" by dividing each "all_raw_react" by div
    df["reactivity"] = df["all_raw_react"] / div
    #make a final dataframe with only nt# counting from 1+, sequence, and reactivity
    final_df = df[["nt","sequence","reactivity"]]
    #replace NA with "NA"
    final_df = final_df.fillna("NA")
    #write final df to csv
    final_df.to_csv(args.output + "_reactivity.csv", sep=",", index=False)

#if nind flag is set, normalize reactivity individually
if args.nind:
    norm_individual(raw_react_df,args.primer)
else:
    norm_together(raw_react_df,args.primer)

    
