##Purpose: Add FDR corrected P-value column to each screen result file
# Directory: cd Desktop/Dropbox/TR_Cas9Paper_CodeCheck
#run: python Scripts/TR_SC345678_FDR.py

import pandas as pd
import numpy as np
import scipy
import statsmodels.stats.multitest as smm
import matplotlib.pyplot as plt
from pandas import DataFrame, Series
import csv

#####################################
#Options:
##fdr: input desired FDR cut off for TRUE/FALSE statments (not neccessary to create column)
##fdrmeth: input desired FDR method (I used Benjamini-Hochberg)
##savedirectory: can specify save directory (Best to save in results for downstream analyses)
fdr = float('0.1')
sfdr = str(fdr)
fdrmeth = 'fdr_bh'
sfdrmeth = str(fdrmeth)
#savedirectory = str('Results/')

###################################

data1 = 'Results/SC3_KpvLib.csv'
data2 = 'Results/SC3_RAvLib.csv'
data3 = 'Results/SC4_KpvLib.csv'
data4 = 'Results/SC4_RAvLib.csv'
data5 = 'Results/SC5_KpvLib.csv'
data6 = 'Results/SC5_RAvLib.csv'
data7 = 'Results/SC6_KpvLib.csv'
data8 = 'Results/SC6_RAvLib.csv'
data9 = 'Results/SC7_KpvLib.csv'
data10 = 'Results/SC7_RAvLib.csv'
data11 = 'Results/SC8_KpvLib.csv'
data12 = 'Results/SC8_RAvLib.csv'



sdata1 = str(data1[0:data1.find('v')])
sdata2 = str(data2[0:data2.find('v')])
sdata3 = str(data3[0:data3.find('v')])
sdata4 = str(data4[0:data4.find('v')])
sdata5 = str(data5[0:data5.find('v')])
sdata6 = str(data6[0:data6.find('v')])
sdata7 = str(data7[0:data7.find('v')])
sdata8 = str(data8[0:data8.find('v')])
sdata9 = str(data9[0:data9.find('v')])
sdata10 = str(data10[0:data10.find('v')])
sdata11 = str(data11[0:data11.find('v')])
sdata12 = str(data12[0:data12.find('v')])


rdata1 = pd.read_csv(data1, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata2 = pd.read_csv(data2, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata3 = pd.read_csv(data3, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata4 = pd.read_csv(data4, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata5 = pd.read_csv(data5, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata6 = pd.read_csv(data6, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata7 = pd.read_csv(data7, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata8 = pd.read_csv(data8, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata9 = pd.read_csv(data9, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata10 = pd.read_csv(data10, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata11 = pd.read_csv(data11, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata12 = pd.read_csv(data12, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])

#
sigs = rdata1['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata1.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata1.to_csv(sdata1 + "_FDR" + sfdr + ".csv")

sigs = rdata2['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata2.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata2.to_csv(sdata2 + "_FDR" + sfdr + ".csv")

sigs = rdata3['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata3.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata3.to_csv(sdata3 + "_FDR" + sfdr + ".csv")

sigs = rdata4['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata4.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata4.to_csv(sdata4 + "_FDR" + sfdr + ".csv")

sigs = rdata5['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata5.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata5.to_csv(sdata5 + "_FDR" + sfdr + ".csv")

sigs = rdata6['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata6.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata6.to_csv(sdata6 + "_FDR" + sfdr + ".csv")

sigs = rdata7['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata7.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata7.to_csv(sdata7 + "_FDR" + sfdr + ".csv")

sigs = rdata8['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata8.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata8.to_csv(sdata8 + "_FDR" + sfdr + ".csv")

sigs = rdata9['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata9.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata9.to_csv(sdata9 + "_FDR" + sfdr + ".csv")

sigs = rdata10['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata10.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata10.to_csv(sdata10 + "_FDR" + sfdr + ".csv")

sigs = rdata11['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata11.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata11.to_csv(sdata11 + "_FDR" + sfdr + ".csv")

sigs = rdata12['casTLE p-value']
multi = smm.multipletests(sigs, alpha=fdr, method=fdrmeth, is_sorted=False, returnsorted=False)
rdata12.insert(loc=9, column="FDR_" + sfdr, value=multi[1])
fincsv = rdata12.to_csv(sdata12 + "_FDR" + sfdr + ".csv")
