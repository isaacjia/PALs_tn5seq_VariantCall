
# coding: utf-8

# In[ ]:




# In[198]:

import pysam
import sys
import os
from collections import Counter
import pandas as pd
import numpy as np


# In[ ]:
bd=sys.argv[1]
infile=sys.argv[2]


# In[199]:

# bd='/Users/Isaac/box/Kitzman_lab/data/20180114_tn5seq/'
# bd='/nfs/kitzman2/isaac/miseq/20171219_miseq/maptoMSH2_NNN501-800/20180119MSH2NNN_out/'
filename = infile.rstrip('.bam')


# In[200]:

WT='ATGGCGGTGCAGCCGAAGGAGACGCTGCAGTTGGAGAGCGCGGCCGAGGTCGGCTTCGTGCGCTTCTTTCAGGGCATGCCGGAGAAGCCGACCACCACAGTGCGCCTTTTCGACCGGGGCGACTTCTATACGGCGCACGGCGAGGACGCGCTGCTGGCCGCCCGGGAGGTGTTCAAGACCCAGGGGGTGATCAAGTACATGGGGCCGGCAGGAGCAAAGAATCTGCAGAGTGTTGTGCTTAGTAAAATGAATTTTGAATCTTTTGTAAAAGATCTTCTTCTGGTTCGTCAGTATAGAGTTGAAGTTTATAAGAATAGAGCTGGAAATAAGGCATCCAAGGAGAATGATTGGTATTTGGCATATAAGGCTTCTCCTGGCAATCTCTCTCAGTTTGAAGACATTCTCTTTGGTAACAATGATATGTCAGCTTCCATTGGTGTTGTGGGTGTTAAAATGTCCGCAGTTGATGGCCAGAGACAGGTTGGAGTTGGGTATGTGGATTCCATACAGAGGAAACTAGGACTGTGTGAATTCCCTGATAATGATCAGTTCTCCAATCTTGAGGCTCTCCTCATCCAGATTGGACCAAAGGAATGTGTTTTACCCGGAGGAGAGACTGCTGGAGACATGGGGAAACTGAGACAGATAATTCAAAGAGGAGGAATTCTGATCACAGAAAGAAAAAAAGCTGACTTTTCCACAAAAGACATTTATCAGGACCTCAACCGGTTGTTGAAAGGCAAAAAGGGAGAGCAGATGAATAGTGCTGTATTGCCAGAAATGGAGAATCAGGTTGCAGTTTCATCACTGTCTGCGGTAATCAAGTTTTTAGAACTCTTATCAGATGATTCCAACTTTGGACAGTTTGAACTGACTACTTTTGACTTCAGCCAGTATATGAAATTGGATATTGCAGCAGTCAGAGCCCTTAACCTTTTTCAGGGTTCTGTTGAAGATACCACTGGCTCTCAGTCTCTGGCTGCCTTGCTGAATAAGTGTAAAACCCCTCAAGGACAAAGACTTGTTAACCAGTGGATTAAGCAGCCTCTCATGGATAAGAACAGAATAGAGGAGAGATTGAATTTAGTGGAAGCTTTTGTAGAAGATGCAGAATTGAGGCAGACTTTACAAGAAGATTTACTTCGTCGATTCCCAGATCTTAACCGACTTGCCAAGAAGTTTCAAAGACAAGCAGCAAACTTACAAGATTGTTACCGACTCTATCAGGGTATAAATCAACTACCTAATGTTATACAGGCTCTGGAAAAACATGAAGGAAAACACCAGAAATTATTGTTGGCAGTTTTTGTGACTCCTCTTACTGATCTTCGTTCTGACTTCTCCAAGTTTCAGGAAATGATAGAAACAACTTTAGATATGGATCAGGTGGAAAACCATGAATTCCTTGTAAAACCTTCATTTGATCCTAATCTCAGTGAATTAAGAGAAATAATGAATGACTTGGAAAAGAAGATGCAGTCAACATTAATAAGTGCAGCCAGAGATCTTGGCTTGGACCCTGGCAAACAGATTAAACTGGATTCCAGTGCACAGTTTGGATATTACTTTCGTGTAACCTGTAAGGAAGAAAAAGTCCTTCGTAACAATAAAAACTTTAGTACTGTAGATATCCAGAAGAATGGTGTTAAATTTACCAACAGCAAATTGACTTCTTTAAATGAAGAGTATACCAAAAATAAAACAGAATATGAAGAAGCCCAGGATGCCATTGTTAAAGAAATTGTCAATATTTCTTCAGGCTATGTAGAACCAATGCAGACACTCAATGATGTGTTAGCTCAGCTAGATGCTGTTGTCAGCTTTGCTCACGTGTCAAATGGAGCACCTGTTCCATATGTACGACCAGCCATTTTGGAGAAAGGACAAGGAAGAATTATATTAAAAGCATCCAGGCATGCTTGTGTTGAAGTTCAAGATGAAATTGCATTTATTCCTAATGACGTATACTTTGAAAAAGATAAACAGATGTTCCACATCATTACTGGCCCCAATATGGGAGGTAAATCAACATATATTCGACAAACTGGGGTGATAGTACTCATGGCCCAAATTGGGTGTTTTGTGCCATGTGAGTCAGCAGAAGTGTCCATTGTGGACTGCATCTTAGCCCGAGTAGGGGCTGGTGACAGTCAATTGAAAGGAGTCTCCACGTTCATGGCTGAAATGTTGGAAACTGCTTCTATCCTCAGGTCTGCAACCAAAGATTCATTAATAATCATAGATGAATTGGGAAGAGGAACTTCTACCTACGATGGATTTGGGTTAGCATGGGCTATATCAGAATACATTGCAACAAAGATTGGTGCTTTTTGCATGTTTGCAACCCATTTTCATGAACTTACTGCCTTGGCCAATCAGATACCAACTGTTAATAATCTACATGTCACAGCACTCACCACTGAAGAGACCTTAACTATGCTTTATCAGGTGAAGAAAGGTGTCTGTGATCAAAGTTTTGGGATTCATGTTGCAGAGCTTGCTAATTTCCCTAAGCATGTAATAGAGTGTGCTAAACAGAAAGCCCTGGAACTTGAGGAGTTTCAGTATATTGGAGAATCGCAAGGATATGATATCATGGAACCAGCAGCAAAGAAGTGCTATCTGGAAAGAGAGCAAGGTGAAAAAATTATTCAGGAGTTCCTGTCCAAGGTGAAACAAATGCCCTTTACTGAAATGTCAGAAGAAAACATCACAATAAAGTTAAAACAGCTAAAAGCTGAAGTAATAGCAAAGAATAATAGCTTTGTAAATGAAATCATTTCACGAATAAAAGTTACTACGTGA'


# In[ ]:




# In[201]:

# for n in range(1758,2007,3):
#     print n/3


# In[202]:

# seq = []
# for n in range(1758,2007,3):
#     Mut_codon = n/3
#     Mut_codon_seq= WT[0:3*Mut_codon]+'NNN'+WT[3*Mut_codon+3:]
#     seq.append(Mut_codon_seq)


# In[203]:

# print seq[82]


# In[204]:

# with open (bd+'MSH2_exon12_NNN.fa', 'w') as NNN:
#     c=586
#     for s in seq:
#         NNN.write('>MSH2_Mut_codon'+str(c)+'\n'+s+'\n')
#         c+=1
# NNN.close()


# In[205]:

def list_to_dict(li):  
    dct = {}  
    for item in li:
        if dct.has_key(item):
            dct[item] = dct[item] + 1  
        else:
            dct[item] = 1
    return dct  


# In[206]:

# bamfile = pysam.AlignmentFile(bd+'tn5seq_XJ_16TG-6p2_dox_Blast__9.bam', "rb")
# WT_depth =0
# Var_list = []
# WT_codon='TCA'
# for read in bamfile:
#     if read.reference_start < 1756 and read.reference_end >1758:
#         left_slice = 1756-read.reference_start
# #         print left_slice
#         read_codon_start = left_slice-1
#         read_codon_end = left_slice+2

#         read_seq = read.query_alignment_sequence
#         codon_seq = read.seq[read_codon_start:read_codon_end]
#         if codon_seq == WT_codon:
#             WT_depth +=1
#         else:
#             Var_list.append(codon_seq)
# total_depth = WT_depth + len(Var_list)
# print WT_depth, len(Var_list)

# bamfile.close()


# In[ ]:




# In[207]:

# codon_num = 586
# bamfile = pysam.AlignmentFile(bd+'tn5seq_XJ_16TG-6p2_dox_Blast__9.bam', "rb")

# codon_list = []
# for read in bamfile:
#     if read.reference_start < codon_num*3 and read.reference_end >codon_num*3+2:
#         left_slice = codon_num*3-read.reference_start
# #         print left_slice
#         read_codon_start = left_slice-1
#         read_codon_end = left_slice+2

#         read_seq = read.query_alignment_sequence
#         codon_seq = read.seq[read_codon_start:read_codon_end]
#         codon_list.append(codon_seq)
# print list_to_dict(codon_list)

# bamfile.close()


# In[208]:

def codon_var_dict_gen(infile, codon_num):
    bamfile = pysam.AlignmentFile(infile, "rb")
    codon_list = []
    for read in bamfile:
        if read.reference_start < (549 + codon_num*3-10) and read.reference_end > (549+ codon_num*3+12):
            left_slice = (549 + codon_num*3)-read.reference_start
#             print left_slice
            read_codon_start = left_slice-3
            read_codon_end = left_slice
    
            read_seq = read.query_alignment_sequence
            codon_seq = read.seq[read_codon_start:read_codon_end]
            codon_list.append(codon_seq)
    return list_to_dict(codon_list)
    bamfile.close()


# In[209]:

# print codon_var_dict_gen(bd+'tn5seq_XJ_16TG-6p2_dox_Blast__9.bam',590)


# In[210]:

# for codon_num in range(580,680,1):
#     print codon_var_dict_gen(bd+'tn5seq_XJ_16TG-6p2_dox_Blast__9.bam',codon_num)
# #     'dict_codon_'+str(codon_num) = codon_var_dict_gen(bd+'tn5seq_XJ_16TG-6p2_dox_Blast__9.bam',codon_num)
    
# #     if codon_num < 590:
# #         print 
# # #     print codon_num
    


# In[211]:

codon64 = ['ATA','ATC','ATT','ATG','ACA','ACC','ACG','ACT','AAC','AAT','AAA','AAG','AGC','AGT','AGA','AGG','CTA','CTC','CTG','CTT','CCA','CCC','CCG','CCT','CAC','CAT','CAA','CAG','CGA','CGC','CGG','CGT','GTA','GTC','GTG','GTT','GCA','GCC','GCG','GCT','GAC','GAT','GAA','GAG','GGA','GGC','GGG','GGT','TCA','TCC','TCG','TCT','TTC','TTT','TTA','TTG','TAC','TAT','TAA','TAG','TGC','TGT','TGA','TGG']
codon_pos_list = []
for codon_pos in range(1,935,1):
    codon_pos_list.append(codon_pos)
df = pd.DataFrame(None, index=codon_pos_list, columns = codon64)
df.head()


# In[212]:

from collections import defaultdict
for codon_num in range(1,935,1):
    codon_dict = codon_var_dict_gen(bd+infile,codon_num)
    for tri_base in codon64:
        if tri_base in codon_dict.keys():
            df.at[codon_num, tri_base] = codon_dict[tri_base]
        else:
            df.at[codon_num, tri_base] = 0


# In[213]:

df.head()


# In[ ]:

df.to_csv(bd+filename+'.csv',sep='\t') #write df to file

