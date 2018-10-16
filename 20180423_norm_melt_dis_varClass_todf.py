
# coding: utf-8

from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

bd=sys.argv[1]

# this example file is necessary. anchor all the files
infile1=sys.argv[2]
# In[305]:

# these variants will affect all later blocks. must define well
# bd='/Users/jia/Desktop/Kitzman_lab/data/20180330_tn5seq/varcall_w/'
# bd='/Users/Isaac/box/Kitzman_lab/data/20180419_tn5seq/realign/'
# infile1= 'XJ_XJ_1.csv'


def normalize(infile):
    df = pd.read_csv(bd+infile, sep='\t') #need bd and infile well defined. 
    
    df = df.rename(columns = {'Unnamed: 0':'codon_pos'})
    df = df.set_index('codon_pos')
    df["sum"] = df.sum(axis=1)
    df["sum"] =df["sum"].astype(float)
    df_norm = df.loc[:,"ATA":"TGG"].div(df["sum"], axis=0)
    return df_norm


def melt(infile): #infile will be used twice in this function
    df = normalize(infile)
    df = df.reset_index('codon_pos')
    AF=infile.rstrip('.csv')
    df_melted = pd.melt(df, id_vars=["codon_pos"], var_name="codon_var", value_name=AF)
    return df_melted


df=melt(infile1)

df_all=df.copy()


# loop through all variant call .csv files and add to the melt dataframe. 
for infile in os.listdir(bd):
#     if infile.startswith('XJ') and '.csv' in infile:
# make sure all the files are .csv before proceeding. 
    if '.csv' in infile:
        print infile
        AF=infile.rstrip('.csv')
#         print AF
        df_all[AF]=melt(infile)[AF]


# In[313]:

# Add a real_index column to the dataframe, to be used for adding the distance column 
# df_all['real_index'] = df_all["codon_pos"].astype(str)+'_' +df_all["codon_var"]
df_all['real_index'] = df_all["codon_pos"].astype(str)+'_' +df_all["codon_var"]
df_all = df_all.set_index('real_index')
df_all.head()





MSH2_cds='ATGGCGGTGCAGCCGAAGGAGACGCTGCAGTTGGAGAGCGCGGCCGAGGTCGGCTTCGTGCGCTTCTTTCAGGGCATGCCGGAGAAGCCGACCACCACAGTGCGCCTTTTCGACCGGGGCGACTTCTATACGGCGCACGGCGAGGACGCGCTGCTGGCCGCCCGGGAGGTGTTCAAGACCCAGGGGGTGATCAAGTACATGGGGCCGGCAGGAGCAAAGAATCTGCAGAGTGTTGTGCTTAGTAAAATGAATTTTGAATCTTTTGTAAAAGATCTTCTTCTGGTTCGTCAGTATAGAGTTGAAGTTTATAAGAATAGAGCTGGAAATAAGGCATCCAAGGAGAATGATTGGTATTTGGCATATAAGGCTTCTCCTGGCAATCTCTCTCAGTTTGAAGACATTCTCTTTGGTAACAATGATATGTCAGCTTCCATTGGTGTTGTGGGTGTTAAAATGTCCGCAGTTGATGGCCAGAGACAGGTTGGAGTTGGGTATGTGGATTCCATACAGAGGAAACTAGGACTGTGTGAATTCCCTGATAATGATCAGTTCTCCAATCTTGAGGCTCTCCTCATCCAGATTGGACCAAAGGAATGTGTTTTACCCGGAGGAGAGACTGCTGGAGACATGGGGAAACTGAGACAGATAATTCAAAGAGGAGGAATTCTGATCACAGAAAGAAAAAAAGCTGACTTTTCCACAAAAGACATTTATCAGGACCTCAACCGGTTGTTGAAAGGCAAAAAGGGAGAGCAGATGAATAGTGCTGTATTGCCAGAAATGGAGAATCAGGTTGCAGTTTCATCACTGTCTGCGGTAATCAAGTTTTTAGAACTCTTATCAGATGATTCCAACTTTGGACAGTTTGAACTGACTACTTTTGACTTCAGCCAGTATATGAAATTGGATATTGCAGCAGTCAGAGCCCTTAACCTTTTTCAGGGTTCTGTTGAAGATACCACTGGCTCTCAGTCTCTGGCTGCCTTGCTGAATAAGTGTAAAACCCCTCAAGGACAAAGACTTGTTAACCAGTGGATTAAGCAGCCTCTCATGGATAAGAACAGAATAGAGGAGAGATTGAATTTAGTGGAAGCTTTTGTAGAAGATGCAGAATTGAGGCAGACTTTACAAGAAGATTTACTTCGTCGATTCCCAGATCTTAACCGACTTGCCAAGAAGTTTCAAAGACAAGCAGCAAACTTACAAGATTGTTACCGACTCTATCAGGGTATAAATCAACTACCTAATGTTATACAGGCTCTGGAAAAACATGAAGGAAAACACCAGAAATTATTGTTGGCAGTTTTTGTGACTCCTCTTACTGATCTTCGTTCTGACTTCTCCAAGTTTCAGGAAATGATAGAAACAACTTTAGATATGGATCAGGTGGAAAACCATGAATTCCTTGTAAAACCTTCATTTGATCCTAATCTCAGTGAATTAAGAGAAATAATGAATGACTTGGAAAAGAAGATGCAGTCAACATTAATAAGTGCAGCCAGAGATCTTGGCTTGGACCCTGGCAAACAGATTAAACTGGATTCCAGTGCACAGTTTGGATATTACTTTCGTGTAACCTGTAAGGAAGAAAAAGTCCTTCGTAACAATAAAAACTTTAGTACTGTAGATATCCAGAAGAATGGTGTTAAATTTACCAACAGCAAATTGACTTCTTTAAATGAAGAGTATACCAAAAATAAAACAGAATATGAAGAAGCCCAGGATGCCATTGTTAAAGAAATTGTCAATATTTCTTCAGGCTATGTAGAACCAATGCAGACACTCAATGATGTGTTAGCTCAGCTAGATGCTGTTGTCAGCTTTGCTCACGTGTCAAATGGAGCACCTGTTCCATATGTACGACCAGCCATTTTGGAGAAAGGACAAGGAAGAATTATATTAAAAGCATCCAGGCATGCTTGTGTTGAAGTTCAAGATGAAATTGCATTTATTCCTAATGACGTATACTTTGAAAAAGATAAACAGATGTTCCACATCATTACTGGCCCCAATATGGGAGGTAAATCAACATATATTCGACAAACTGGGGTGATAGTACTCATGGCCCAAATTGGGTGTTTTGTGCCATGTGAGTCAGCAGAAGTGTCCATTGTGGACTGCATCTTAGCCCGAGTAGGGGCTGGTGACAGTCAATTGAAAGGAGTCTCCACGTTCATGGCTGAAATGTTGGAAACTGCTTCTATCCTCAGGTCTGCAACCAAAGATTCATTAATAATCATAGATGAATTGGGAAGAGGAACTTCTACCTACGATGGATTTGGGTTAGCATGGGCTATATCAGAATACATTGCAACAAAGATTGGTGCTTTTTGCATGTTTGCAACCCATTTTCATGAACTTACTGCCTTGGCCAATCAGATACCAACTGTTAATAATCTACATGTCACAGCACTCACCACTGAAGAGACCTTAACTATGCTTTATCAGGTGAAGAAAGGTGTCTGTGATCAAAGTTTTGGGATTCATGTTGCAGAGCTTGCTAATTTCCCTAAGCATGTAATAGAGTGTGCTAAACAGAAAGCCCTGGAACTTGAGGAGTTTCAGTATATTGGAGAATCGCAAGGATATGATATCATGGAACCAGCAGCAAAGAAGTGCTATCTGGAAAGAGAGCAAGGTGAAAAAATTATTCAGGAGTTCCTGTCCAAGGTGAAACAAATGCCCTTTACTGAAATGTCAGAAGAAAACATCACAATAAAGTTAAAACAGCTAAAAGCTGAAGTAATAGCAAAGAATAATAGCTTTGTAAATGAAATCATTTCACGAATAAAAGTTACTACGTGA'
codon64 = ['ATA','ATC','ATT','ATG','ACA','ACC','ACG','ACT','AAC','AAT','AAA','AAG','AGC','AGT','AGA','AGG','CTA','CTC','CTG','CTT','CCA','CCC','CCG','CCT','CAC','CAT','CAA','CAG','CGA','CGC','CGG','CGT','GTA','GTC','GTG','GTT','GCA','GCC','GCG','GCT','GAC','GAT','GAA','GAG','GGA','GGC','GGG','GGT','TCA','TCC','TCG','TCT','TTC','TTT','TTA','TTG','TAC','TAT','TAA','TAG','TGC','TGT','TGA','TGG']


# In[315]:

MSH2_codonWT ={}
for codon_num in range(1,936,1):

    cpWT ='codon'+str(codon_num)+'WT'
    cpWT = MSH2_cds[codon_num*3-3:codon_num*3]
    MSH2_codonWT[codon_num]=cpWT
# print MSH2_codonWT
# print MSH2_codonWT[1],MSH2_codonWT[934]


# In[316]:

# idea is loop through 64 codon as codon 1, and loop through as codon2. 
# compare each of the 3 elements in one particular codon, counter go up 1 only if they are different
# write each codon pair as key of the dictionay and number of mutations as value. 

mut_number = {}
for codon1 in codon64:
    for codon2 in codon64:
        mut = codon1 + 'to'+codon2
        c = 0
        for b in range(0,3,1):
            if codon1[b] != codon2[b]:
                c+=1
        mut_number[mut] = c
# print mut_number


# In[317]:

codon_AA = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}


# In[318]:

all_codon_dict = {'I':['ATT', 'ATC', 'ATA'],'L':['TTG','CTT', 'CTC', 'CTA', 'CTG', 'TTA',], 'V':['GTT', 'GTC', 'GTA', 'GTG'], 'F':['TTT', 'TTC'], 'M':['ATG'], 'C':['TGT', 'TGC'], 'A':['GCT', 'GCC', 'GCA', 'GCG'], 'G':['GGT', 'GGC', 'GGA', 'GGG'], 'P':['CCA', 'CCC', 'CCT', 'CCG'], 'T':['ACT', 'ACC', 'ACA', 'ACG'], 'S':['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], 'Y':['TAT', 'TAC'], 'W':['TGG'], 'Q':['CAA', 'CAG'], 'N':['AAT', 'AAC'], 'H':['CAT', 'CAC'], 'E':['GAA', 'GAG'],'D':['GAT', 'GAC'], 'K':['AAA', 'AAG'], 'R':['AGA', 'CGT', 'CGC', 'CGA', 'CGG', 'AGG'],'*':['TAA','TGA','TAG']}


# In[319]:

codon_syn_dict ={}
for codon in codon_AA.keys():
    AA = codon_AA[codon]
    AA_all_codon = all_codon_dict[AA]
    codon_syn=[]
    for i in AA_all_codon:
        if i != codon:
            codon_syn.append(i)
    codon_syn_dict[codon]=codon_syn
# print codon_syn_dict


for i in df_all.index:
    codon_pos = df_all.at[i,"codon_pos"]
    codon_var = df_all.at[i,"codon_var"]
    
    WT_codon = MSH2_codonWT[codon_pos]
    dict_key = WT_codon + 'to' + codon_var
    mut_distance = mut_number[dict_key]
    df_all.at[i,"mut_distance"] = int(mut_distance)
    
    if codon_var == WT_codon:
        var_class = 'WT'
    elif codon_AA[codon_var] == codon_AA[WT_codon]:
        var_class = 'Syns'
    elif codon_var in ["TAA","TAG","TGA"]:
        var_class = 'STOP'
    else:
        var_class = 'VUS'
        
    
    df_all.at[i,"var_class"] = var_class
    


df_all.to_csv(bd+'norm_melt_dis_varClass.csv',sep='\t', index=True) #write df to file


# In[ ]:



