#concat screen files for PCA
#make sure you are in the directory containing csv casTLE files: cd Desktop/Dropbox/TR_Cas9Paper_CodeCheck/Results
#Script must also be placed in directory

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np

output = "Results/SC3456_sharedsigs/"

data1 = 'Results/SC3_KpvLib.csv'
data2 = 'Results/SC3_RAvLib.csv'
data3 = 'Results/SC4_KpvLib.csv'
data4 = 'Results/SC4_RAvLib.csv'
data5 = 'Results/SC5_KpvLib.csv'
data6 = 'Results/SC5_RAvLib.csv'
data7 = 'Results/SC6_KpvLib.csv'
data8 = 'Results/SC6_RAvLib.csv'


stdata1 = str(data1[data1.find('SC'):data1.find('v')])
stdata2 = str(data2[data1.find('SC'):data2.find('v')])
stdata3 = str(data3[data1.find('SC'):data3.find('v')])
stdata4 = str(data4[data1.find('SC'):data4.find('v')])
stdata5 = str(data5[data1.find('SC'):data5.find('v')])
stdata6 = str(data6[data1.find('SC'):data6.find('v')])
stdata7 = str(data7[data1.find('SC'):data7.find('v')])
stdata8 = str(data8[data1.find('SC'):data8.find('v')])



rdata1 = pd.read_csv(data1, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata2 = pd.read_csv(data2, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata3 = pd.read_csv(data3, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata4 = pd.read_csv(data4, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata5 = pd.read_csv(data5, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata6 = pd.read_csv(data6, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata7 = pd.read_csv(data7, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata8 = pd.read_csv(data8, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])


concat = pd.concat([rdata1['Symbol'], (rdata1['casTLE Score'].rename(stdata1)), (rdata2['casTLE Score'].rename(stdata2)), (rdata3['casTLE Score'].rename(stdata3)), (rdata4['casTLE Score'].rename(stdata4)), (rdata5['casTLE Score'].rename(stdata5)),
(rdata6['casTLE Score'].rename(stdata6)), (rdata7['casTLE Score'].rename(stdata7)), (rdata8['casTLE Score'].rename(stdata8))], axis=1, ignore_index=False, join='inner')
final = concat.to_csv(output + 'SC3456_castlescores.csv')
