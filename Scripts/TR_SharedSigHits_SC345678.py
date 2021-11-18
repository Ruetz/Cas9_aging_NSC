# Make a list of significant hits (specifying FDR) in 2 or more screens from SC345678
# cd Desktop/Dropbox/TR_Cas9_CodeCheck
# run: python SC345678_sharedsigs_020520.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame, Series


#specify fdr
fdr = float(0.1)
sfdr = str("_fdr" + str(fdr))

#specify output folder directory
output = "Results/SC345678_sharedsigs/"

#data import
data1 = 'Results/SC3_Kp_FDR0.1.csv'
data2 = 'Results/SC3_RA_FDR0.1.csv'
data3 = 'Results/SC4_Kp_FDR0.1.csv'
data4 = 'Results/SC4_RA_FDR0.1.csv'
data5 = 'Results/SC5_Kp_FDR0.1.csv'
data6 = 'Results/SC5_RA_FDR0.1.csv'
data7 = 'Results/SC6_Kp_FDR0.1.csv'
data8 = 'Results/SC6_RA_FDR0.1.csv'
data9 = 'Results/Gene_list.csv'
data10 = 'Results/SC7_Kp_FDR0.1.csv'
data11 = 'Results/SC7_RA_FDR0.1.csv'
data12 = 'Results/SC8_Kp_FDR0.1.csv'
data13 = 'Results/SC8_RA_FDR0.1.csv'



#Capture file names for columns
stdata1 = str(data1[8:data1.find('F')])
stdata2 = str(data2[8:data2.find('F')])
stdata3 = str(data3[8:data3.find('F')])
stdata4 = str(data4[8:data4.find('F')])
stdata5 = str(data5[8:data5.find('F')])
stdata6 = str(data6[8:data6.find('F')])
stdata7 = str(data7[8:data7.find('F')])
stdata8 = str(data8[8:data8.find('F')])
stdata10 = str(data10[8:data10.find('F')])
stdata11 = str(data11[8:data11.find('F')])
stdata12 = str(data12[8:data12.find('F')])
stdata13 = str(data13[8:data13.find('F')])


##input files
rdata1 = pd.read_csv(data1, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata2 = pd.read_csv(data2, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata3 = pd.read_csv(data3, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata4 = pd.read_csv(data4, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata5 = pd.read_csv(data5, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata6 = pd.read_csv(data6, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata7 = pd.read_csv(data7, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata8 = pd.read_csv(data8, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata9 = pd.read_csv(data9, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata10 = pd.read_csv(data10, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata11 = pd.read_csv(data11, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata12 = pd.read_csv(data12, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata13 = pd.read_csv(data13, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])


#narrow to only FDR_X significant gene list
d1pval = rdata1[rdata1['FDR_0.1'] <= fdr]
d2pval = rdata2[rdata2['FDR_0.1'] <= fdr]
d3pval = rdata3[rdata3['FDR_0.1'] <= fdr]
d4pval = rdata4[rdata4['FDR_0.1'] <= fdr]
d5pval = rdata5[rdata5['FDR_0.1'] <= fdr]
d6pval = rdata6[rdata6['FDR_0.1'] <= fdr]
d7pval = rdata7[rdata7['FDR_0.1'] <= fdr]
d8pval = rdata8[rdata8['FDR_0.1'] <= fdr]
d10pval = rdata10[rdata10['FDR_0.1'] <= fdr]
d11pval = rdata11[rdata11['FDR_0.1'] <= fdr]
d12pval = rdata12[rdata12['FDR_0.1'] <= fdr]
d13pval = rdata13[rdata13['FDR_0.1'] <= fdr]


#Write to csv file
sigsd1 = d1pval.to_csv(output + stdata1 + sfdr + ".csv")
sigsd2 = d2pval.to_csv(output + stdata2 + sfdr + ".csv")
sigsd3 = d3pval.to_csv(output + stdata3 + sfdr + ".csv")
sigsd4 = d4pval.to_csv(output + stdata4 + sfdr + ".csv")
sigsd5 = d5pval.to_csv(output + stdata5 + sfdr + ".csv")
sigsd6 = d6pval.to_csv(output + stdata6 + sfdr + ".csv")
sigsd7 = d7pval.to_csv(output + stdata7 + sfdr + ".csv")
sigsd8 = d8pval.to_csv(output + stdata8 + sfdr + ".csv")
sigsd10 = d10pval.to_csv(output + stdata10 + sfdr + ".csv")
sigsd11 = d11pval.to_csv(output + stdata11 + sfdr + ".csv")
sigsd12 = d12pval.to_csv(output + stdata12 + sfdr + ".csv")
sigsd13 = d13pval.to_csv(output + stdata13 + sfdr + ".csv")


##first read csv files
rsigsd1 = pd.read_csv(output + stdata1 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd2 = pd.read_csv(output + stdata2 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd3 = pd.read_csv(output + stdata3 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd4 = pd.read_csv(output + stdata4 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd5 = pd.read_csv(output + stdata5 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd6 = pd.read_csv(output + stdata6 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd7 = pd.read_csv(output + stdata7 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd8 = pd.read_csv(output + stdata8 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd10 = pd.read_csv(output + stdata10 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd11 = pd.read_csv(output + stdata11 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd12 = pd.read_csv(output + stdata12 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd13 = pd.read_csv(output + stdata13 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')


#narrow to only enriched values
rsigsd1_enr = rsigsd1[rsigsd1['casTLE Effect'] > 0]
rsigsd2_enr = rsigsd2[rsigsd2['casTLE Effect'] > 0]
rsigsd3_enr = rsigsd3[rsigsd3['casTLE Effect'] > 0]
rsigsd4_enr = rsigsd4[rsigsd4['casTLE Effect'] > 0]
rsigsd5_enr = rsigsd5[rsigsd5['casTLE Effect'] > 0]
rsigsd6_enr = rsigsd6[rsigsd6['casTLE Effect'] > 0]
rsigsd7_enr = rsigsd7[rsigsd7['casTLE Effect'] > 0]
rsigsd8_enr = rsigsd8[rsigsd8['casTLE Effect'] > 0]
rsigsd10_enr = rsigsd10[rsigsd10['casTLE Effect'] > 0]
rsigsd11_enr = rsigsd11[rsigsd11['casTLE Effect'] > 0]
rsigsd12_enr = rsigsd12[rsigsd12['casTLE Effect'] > 0]
rsigsd13_enr = rsigsd13[rsigsd13['casTLE Effect'] > 0]

###then concat
#Old
concat3Kp_6Kp = pd.merge(rsigsd1_enr, rsigsd7_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3Kp_6RA = pd.merge(rsigsd1_enr, rsigsd8_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_6Kp = pd.merge(rsigsd2_enr, rsigsd7_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_6RA = pd.merge(rsigsd2_enr, rsigsd8_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3Kp_8Kp = pd.merge(rsigsd1_enr, rsigsd12_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3Kp_8RA = pd.merge(rsigsd1_enr, rsigsd13_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_8Kp = pd.merge(rsigsd2_enr, rsigsd12_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_8RA = pd.merge(rsigsd2_enr, rsigsd13_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat6Kp_8Kp = pd.merge(rsigsd7_enr, rsigsd12_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat6Kp_8RA = pd.merge(rsigsd7_enr, rsigsd13_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat6RA_8Kp = pd.merge(rsigsd8_enr, rsigsd12_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat6RA_8RA = pd.merge(rsigsd8_enr, rsigsd13_enr, on='Symbol', how='inner').drop_duplicates(keep='first')

#Young
concat4Kp_5Kp = pd.merge(rsigsd3_enr, rsigsd5_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4Kp_5RA = pd.merge(rsigsd3_enr, rsigsd6_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_5Kp = pd.merge(rsigsd4_enr, rsigsd5_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_5RA = pd.merge(rsigsd4_enr, rsigsd6_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4Kp_7Kp = pd.merge(rsigsd3_enr, rsigsd10_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4Kp_7RA = pd.merge(rsigsd3_enr, rsigsd11_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_7Kp = pd.merge(rsigsd4_enr, rsigsd10_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_7RA = pd.merge(rsigsd4_enr, rsigsd11_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat5Kp_7Kp = pd.merge(rsigsd5_enr, rsigsd10_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat5Kp_7RA = pd.merge(rsigsd5_enr, rsigsd11_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat5RA_7Kp = pd.merge(rsigsd6_enr, rsigsd10_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat5RA_7RA = pd.merge(rsigsd6_enr, rsigsd11_enr, on='Symbol', how='inner').drop_duplicates(keep='first')


##write individual concats to file
wconcat3Kp_6Kp = concat3Kp_6Kp.to_csv(output + stdata1 + stdata7 + sfdr + ".csv")
wconcat3Kp_6RA = concat3Kp_6RA.to_csv(output + stdata1 + stdata8 + sfdr + ".csv")
wconcat3RA_6Kp = concat3RA_6Kp.to_csv(output + stdata2 + stdata7 + sfdr + ".csv")
wconcat3RA_6RA = concat3RA_6RA.to_csv(output + stdata2 + stdata8 + sfdr + ".csv")
wconcat3Kp_8Kp = concat3Kp_8Kp.to_csv(output + stdata1 + stdata12 + sfdr + ".csv")
wconcat3Kp_8RA = concat3Kp_8RA.to_csv(output + stdata1 + stdata13 + sfdr + ".csv")
wconcat3RA_8Kp = concat3RA_8Kp.to_csv(output + stdata2 + stdata12 + sfdr + ".csv")
wconcat3RA_8RA = concat3RA_8RA.to_csv(output + stdata2 + stdata13 + sfdr + ".csv")
wconcat6Kp_8Kp = concat6Kp_8Kp.to_csv(output + stdata7 + stdata12 + sfdr + ".csv")
wconcat6Kp_8RA = concat6Kp_8RA.to_csv(output + stdata7 + stdata13 + sfdr + ".csv")
wconcat6RA_8Kp = concat6RA_8Kp.to_csv(output + stdata8 + stdata12 + sfdr + ".csv")
wconcat6RA_8RA = concat6RA_8RA.to_csv(output + stdata8 + stdata13 + sfdr + ".csv")

wconcat4Kp_5Kp = concat4Kp_5Kp.to_csv(output + stdata3 + stdata5 + sfdr + ".csv")
wconcat4Kp_5RA = concat4Kp_5RA.to_csv(output + stdata3 + stdata6 + sfdr + ".csv")
wconcat4RA_5Kp = concat4RA_5Kp.to_csv(output + stdata4 + stdata5 + sfdr + ".csv")
wconcat4RA_5RA = concat4RA_5RA.to_csv(output + stdata4 + stdata6 + sfdr + ".csv")
wconcat4Kp_7Kp = concat4Kp_7Kp.to_csv(output + stdata3 + stdata10 + sfdr + ".csv")
wconcat4Kp_7RA = concat4Kp_7RA.to_csv(output + stdata3 + stdata11 + sfdr + ".csv")
wconcat4RA_7Kp = concat4RA_7Kp.to_csv(output + stdata4 + stdata10 + sfdr + ".csv")
wconcat4RA_7RA = concat4RA_7RA.to_csv(output + stdata4 + stdata11 + sfdr + ".csv")
wconcat5Kp_7Kp = concat5Kp_7Kp.to_csv(output + stdata5 + stdata10 + sfdr + ".csv")
wconcat5Kp_7RA = concat5Kp_7RA.to_csv(output + stdata5 + stdata11 + sfdr + ".csv")
wconcat5RA_7Kp = concat5RA_7Kp.to_csv(output + stdata6 + stdata10 + sfdr + ".csv")
wconcat5RA_7RA = concat5RA_7RA.to_csv(output + stdata6 + stdata11 + sfdr + ".csv")


##read files
rwconcat3Kp_6Kp = pd.read_csv(output + stdata1 + stdata7 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3Kp_6RA = pd.read_csv(output + stdata1 + stdata8 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_6Kp = pd.read_csv(output + stdata2 + stdata7 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_6RA = pd.read_csv(output + stdata2 + stdata8 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3Kp_8Kp = pd.read_csv(output + stdata1 + stdata12 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3Kp_8RA = pd.read_csv(output + stdata1 + stdata13 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_8Kp = pd.read_csv(output + stdata2 + stdata12 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_8RA = pd.read_csv(output + stdata2 + stdata13 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat6Kp_8Kp = pd.read_csv(output + stdata7 + stdata12 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat6Kp_8RA = pd.read_csv(output + stdata7 + stdata13 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat6RA_8Kp = pd.read_csv(output + stdata8 + stdata12 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat6RA_8RA = pd.read_csv(output + stdata8 + stdata13 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

rwconcat4Kp_5Kp = pd.read_csv(output + stdata3 + stdata5 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4Kp_5RA = pd.read_csv(output + stdata3 + stdata6 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_5Kp = pd.read_csv(output + stdata4 + stdata5 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_5RA = pd.read_csv(output + stdata4 + stdata6 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4Kp_7Kp = pd.read_csv(output + stdata3 + stdata10 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4Kp_7RA = pd.read_csv(output + stdata3 + stdata11 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_7Kp = pd.read_csv(output + stdata4 + stdata10 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_7RA = pd.read_csv(output + stdata4 + stdata11 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat5Kp_7Kp = pd.read_csv(output + stdata5 + stdata10 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat5Kp_7RA = pd.read_csv(output + stdata5 + stdata11 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat5RA_7Kp = pd.read_csv(output + stdata6 + stdata10 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat5RA_7RA = pd.read_csv(output + stdata6 + stdata11 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#make a master sigs files
#df1.merge(df2,on='name').merge(df3,on='name')
concat345678 = rwconcat3Kp_6Kp.merge(rwconcat3Kp_6RA, on='Symbol', how='outer').merge(rwconcat3RA_6Kp, on='Symbol', how='outer').merge(rwconcat3RA_6RA,
on='Symbol', how='outer').merge(rwconcat4Kp_5Kp, on='Symbol', how='outer').merge(rwconcat4Kp_5RA, on='Symbol',
how='outer').merge(rwconcat4RA_5Kp, on='Symbol', how='outer').merge(rwconcat4RA_5RA, on='Symbol', how='outer').merge(rwconcat3Kp_8Kp, on='Symbol', how='outer').merge(rwconcat3Kp_8RA, on='Symbol',
how='outer').merge(rwconcat3RA_8Kp, on='Symbol', how='outer').merge(rwconcat3RA_8RA,
on='Symbol', how='outer').merge(rwconcat6Kp_8Kp, on='Symbol', how='outer').merge(rwconcat6Kp_8RA, on='Symbol', how='outer').merge(rwconcat6RA_8Kp, on='Symbol',
how='outer').merge(rwconcat6RA_8RA, on='Symbol', how='outer').merge(rwconcat4Kp_7Kp, on='Symbol', how='outer').merge(rwconcat4Kp_7RA, on='Symbol',
how='outer').merge(rwconcat4RA_7Kp, on='Symbol', how='outer').merge(rwconcat4RA_7RA, on='Symbol', how='outer').merge(rwconcat5Kp_7Kp, on='Symbol',
how='outer').merge(rwconcat5Kp_7RA, on='Symbol', how='outer').merge(rwconcat5RA_7Kp, on='Symbol', how='outer').merge(rwconcat5RA_7RA, on='Symbol', how='outer')
#concat3456 = pd.concat([, , , , , , , ], axis=1, join='outer').drop_duplicates(keep='first')
Sharedlist = concat345678.to_csv(output + "SC345678_2ormorescreens_all_enriched" + sfdr + ".csv", encoding='utf-8')
rsigs345678 = pd.read_csv(output + "SC345678_2ormorescreens_all_enriched" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#concat3456_old = pd.concat([rwconcat3Kp_6Kp, rwconcat3Kp_6RA, rwconcat3RA_6Kp, rwconcat3RA_6RA], axis=1, join='outer').drop_duplicates(keep='first')
concat345678_old = rwconcat3Kp_6Kp.merge(rwconcat3Kp_6RA, on='Symbol', how='outer').merge(rwconcat3RA_6Kp, on='Symbol', how='outer').merge(rwconcat3RA_6RA, on='Symbol', how='outer').merge(rwconcat3Kp_8Kp, on='Symbol', how='outer').merge(rwconcat3Kp_8RA, on='Symbol',
how='outer').merge(rwconcat3RA_8Kp, on='Symbol', how='outer').merge(rwconcat3RA_8RA, on='Symbol', how='outer').merge(rwconcat6Kp_8Kp, on='Symbol', how='outer').merge(rwconcat6Kp_8RA, on='Symbol', how='outer').merge(rwconcat6RA_8Kp, on='Symbol',
how='outer').merge(rwconcat6RA_8RA, on='Symbol', how='outer')
Sharedlist_old = concat345678_old.to_csv(output + "SC345678_2ormorescreens_all_enriched_OLD_" + sfdr + ".csv", encoding='utf-8')
rsigs345678_old = pd.read_csv(output + "SC345678_2ormorescreens_all_enriched_OLD_" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#concat3456_young = pd.concat([rwconcat4Kp_5Kp, rwconcat4Kp_5RA, rwconcat4RA_5Kp, rwconcat4RA_5RA], axis=1, join='outer').drop_duplicates(keep='first')
concat345678_young = rwconcat4Kp_5Kp.merge(rwconcat4Kp_5RA, on='Symbol', how='outer').merge(rwconcat4RA_5Kp, on='Symbol', how='outer').merge(rwconcat4RA_5RA, on='Symbol', how='outer').merge(rwconcat4Kp_7Kp, on='Symbol', how='outer').merge(rwconcat4Kp_7RA, on='Symbol',
how='outer').merge(rwconcat4RA_7Kp, on='Symbol', how='outer').merge(rwconcat4RA_7RA, on='Symbol', how='outer').merge(rwconcat5Kp_7Kp, on='Symbol', how='outer').merge(rwconcat5Kp_7RA, on='Symbol', how='outer').merge(rwconcat5RA_7Kp, on='Symbol',
how='outer').merge(rwconcat5RA_7RA, on='Symbol', how='outer')
Sharedlist_young = concat345678_young.to_csv(output + "SC345678_2ormorescreens_all_enriched_YOUNG_" + sfdr + ".csv", encoding='utf-8')
rsigs345678_young = pd.read_csv(output + "SC345678_2ormorescreens_all_enriched_YOUNG_" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

##make final lists by forming a union of young and old, and then joing union with young or old as outer merge
#make union
union_list = pd.merge(rsigs345678_young, rsigs345678_old, on='Symbol', how='inner')
writefinal_union = union_list.to_csv(output + "SC345678_Union_enriched" + sfdr + ".csv", encoding='utf-8')
readfinal_union = pd.read_csv(output + "SC345678_Union_enriched" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

finallist_old = pd.merge(readfinal_union, rsigs345678_old, on='Symbol', how='outer', indicator=True).query('_merge=="right_only"')
writefinal_old = finallist_old.to_csv(output + "SC345678_Oldspecific_all_enriched" + sfdr + ".csv", encoding='utf-8')
readfinal_old = pd.read_csv(output + "SC345678_Oldspecific_all_enriched" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

finallist_young = pd.merge(readfinal_union, rsigs345678_young, on='Symbol', how='outer', indicator=True).query('_merge=="right_only"')
writefinal_young = finallist_young.to_csv(output + "SC345678_Youngspecific_all_enriched" + sfdr + ".csv", encoding='utf-8')
readfinal_young = pd.read_csv(output + "SC345678_Youngspecific_all_enriched" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#############################################
#############################################
#depleted list
#narrow to only enriched values
rsigsd1_enr = rsigsd1[rsigsd1['casTLE Effect'] < 0]
rsigsd2_enr = rsigsd2[rsigsd2['casTLE Effect'] < 0]
rsigsd3_enr = rsigsd3[rsigsd3['casTLE Effect'] < 0]
rsigsd4_enr = rsigsd4[rsigsd4['casTLE Effect'] < 0]
rsigsd5_enr = rsigsd5[rsigsd5['casTLE Effect'] < 0]
rsigsd6_enr = rsigsd6[rsigsd6['casTLE Effect'] < 0]
rsigsd7_enr = rsigsd7[rsigsd7['casTLE Effect'] < 0]
rsigsd8_enr = rsigsd8[rsigsd8['casTLE Effect'] < 0]
rsigsd10_enr = rsigsd10[rsigsd10['casTLE Effect'] < 0]
rsigsd11_enr = rsigsd11[rsigsd11['casTLE Effect'] < 0]
rsigsd12_enr = rsigsd12[rsigsd12['casTLE Effect'] < 0]
rsigsd13_enr = rsigsd13[rsigsd13['casTLE Effect'] < 0]

###then concat
#Old
concat3Kp_6Kp = pd.merge(rsigsd1_enr, rsigsd7_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3Kp_6RA = pd.merge(rsigsd1_enr, rsigsd8_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_6Kp = pd.merge(rsigsd2_enr, rsigsd7_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_6RA = pd.merge(rsigsd2_enr, rsigsd8_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3Kp_8Kp = pd.merge(rsigsd1_enr, rsigsd12_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3Kp_8RA = pd.merge(rsigsd1_enr, rsigsd13_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_8Kp = pd.merge(rsigsd2_enr, rsigsd12_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_8RA = pd.merge(rsigsd2_enr, rsigsd13_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat6Kp_8Kp = pd.merge(rsigsd7_enr, rsigsd12_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat6Kp_8RA = pd.merge(rsigsd7_enr, rsigsd13_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat6RA_8Kp = pd.merge(rsigsd8_enr, rsigsd12_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat6RA_8RA = pd.merge(rsigsd8_enr, rsigsd13_enr, on='Symbol', how='inner').drop_duplicates(keep='first')

#Young
concat4Kp_5Kp = pd.merge(rsigsd3_enr, rsigsd5_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4Kp_5RA = pd.merge(rsigsd3_enr, rsigsd6_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_5Kp = pd.merge(rsigsd4_enr, rsigsd5_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_5RA = pd.merge(rsigsd4_enr, rsigsd6_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4Kp_7Kp = pd.merge(rsigsd3_enr, rsigsd10_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4Kp_7RA = pd.merge(rsigsd3_enr, rsigsd11_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_7Kp = pd.merge(rsigsd4_enr, rsigsd10_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_7RA = pd.merge(rsigsd4_enr, rsigsd11_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat5Kp_7Kp = pd.merge(rsigsd5_enr, rsigsd10_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat5Kp_7RA = pd.merge(rsigsd5_enr, rsigsd11_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat5RA_7Kp = pd.merge(rsigsd6_enr, rsigsd10_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat5RA_7RA = pd.merge(rsigsd6_enr, rsigsd11_enr, on='Symbol', how='inner').drop_duplicates(keep='first')


##write individual concats to file
wconcat3Kp_6Kp = concat3Kp_6Kp.to_csv(output + stdata1 + stdata7 + sfdr + ".csv")
wconcat3Kp_6RA = concat3Kp_6RA.to_csv(output + stdata1 + stdata8 + sfdr + ".csv")
wconcat3RA_6Kp = concat3RA_6Kp.to_csv(output + stdata2 + stdata7 + sfdr + ".csv")
wconcat3RA_6RA = concat3RA_6RA.to_csv(output + stdata2 + stdata8 + sfdr + ".csv")
wconcat3Kp_8Kp = concat3Kp_8Kp.to_csv(output + stdata1 + stdata12 + sfdr + ".csv")
wconcat3Kp_8RA = concat3Kp_8RA.to_csv(output + stdata1 + stdata13 + sfdr + ".csv")
wconcat3RA_8Kp = concat3RA_8Kp.to_csv(output + stdata2 + stdata12 + sfdr + ".csv")
wconcat3RA_8RA = concat3RA_8RA.to_csv(output + stdata2 + stdata13 + sfdr + ".csv")
wconcat6Kp_8Kp = concat6Kp_8Kp.to_csv(output + stdata7 + stdata12 + sfdr + ".csv")
wconcat6Kp_8RA = concat6Kp_8RA.to_csv(output + stdata7 + stdata13 + sfdr + ".csv")
wconcat6RA_8Kp = concat6RA_8Kp.to_csv(output + stdata8 + stdata12 + sfdr + ".csv")
wconcat6RA_8RA = concat6RA_8RA.to_csv(output + stdata8 + stdata13 + sfdr + ".csv")

wconcat4Kp_5Kp = concat4Kp_5Kp.to_csv(output + stdata3 + stdata5 + sfdr + ".csv")
wconcat4Kp_5RA = concat4Kp_5RA.to_csv(output + stdata3 + stdata6 + sfdr + ".csv")
wconcat4RA_5Kp = concat4RA_5Kp.to_csv(output + stdata4 + stdata5 + sfdr + ".csv")
wconcat4RA_5RA = concat4RA_5RA.to_csv(output + stdata4 + stdata6 + sfdr + ".csv")
wconcat4Kp_7Kp = concat4Kp_7Kp.to_csv(output + stdata3 + stdata10 + sfdr + ".csv")
wconcat4Kp_7RA = concat4Kp_7RA.to_csv(output + stdata3 + stdata11 + sfdr + ".csv")
wconcat4RA_7Kp = concat4RA_7Kp.to_csv(output + stdata4 + stdata10 + sfdr + ".csv")
wconcat4RA_7RA = concat4RA_7RA.to_csv(output + stdata4 + stdata11 + sfdr + ".csv")
wconcat5Kp_7Kp = concat5Kp_7Kp.to_csv(output + stdata5 + stdata10 + sfdr + ".csv")
wconcat5Kp_7RA = concat5Kp_7RA.to_csv(output + stdata5 + stdata11 + sfdr + ".csv")
wconcat5RA_7Kp = concat5RA_7Kp.to_csv(output + stdata6 + stdata10 + sfdr + ".csv")
wconcat5RA_7RA = concat5RA_7RA.to_csv(output + stdata6 + stdata11 + sfdr + ".csv")


##read files
rwconcat3Kp_6Kp = pd.read_csv(output + stdata1 + stdata7 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3Kp_6RA = pd.read_csv(output + stdata1 + stdata8 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_6Kp = pd.read_csv(output + stdata2 + stdata7 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_6RA = pd.read_csv(output + stdata2 + stdata8 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3Kp_8Kp = pd.read_csv(output + stdata1 + stdata12 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3Kp_8RA = pd.read_csv(output + stdata1 + stdata13 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_8Kp = pd.read_csv(output + stdata2 + stdata12 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_8RA = pd.read_csv(output + stdata2 + stdata13 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat6Kp_8Kp = pd.read_csv(output + stdata7 + stdata12 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat6Kp_8RA = pd.read_csv(output + stdata7 + stdata13 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat6RA_8Kp = pd.read_csv(output + stdata8 + stdata12 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat6RA_8RA = pd.read_csv(output + stdata8 + stdata13 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

rwconcat4Kp_5Kp = pd.read_csv(output + stdata3 + stdata5 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4Kp_5RA = pd.read_csv(output + stdata3 + stdata6 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_5Kp = pd.read_csv(output + stdata4 + stdata5 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_5RA = pd.read_csv(output + stdata4 + stdata6 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4Kp_7Kp = pd.read_csv(output + stdata3 + stdata10 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4Kp_7RA = pd.read_csv(output + stdata3 + stdata11 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_7Kp = pd.read_csv(output + stdata4 + stdata10 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_7RA = pd.read_csv(output + stdata4 + stdata11 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat5Kp_7Kp = pd.read_csv(output + stdata5 + stdata10 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat5Kp_7RA = pd.read_csv(output + stdata5 + stdata11 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat5RA_7Kp = pd.read_csv(output + stdata6 + stdata10 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat5RA_7RA = pd.read_csv(output + stdata6 + stdata11 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#make a master sigs files
#df1.merge(df2,on='name').merge(df3,on='name')
concat345678 = rwconcat3Kp_6Kp.merge(rwconcat3Kp_6RA, on='Symbol', how='outer').merge(rwconcat3RA_6Kp, on='Symbol', how='outer').merge(rwconcat3RA_6RA,
on='Symbol', how='outer').merge(rwconcat4Kp_5Kp, on='Symbol', how='outer').merge(rwconcat4Kp_5RA, on='Symbol',
how='outer').merge(rwconcat4RA_5Kp, on='Symbol', how='outer').merge(rwconcat4RA_5RA, on='Symbol', how='outer').merge(rwconcat3Kp_8Kp, on='Symbol', how='outer').merge(rwconcat3Kp_8RA, on='Symbol',
how='outer').merge(rwconcat3RA_8Kp, on='Symbol', how='outer').merge(rwconcat3RA_8RA,
on='Symbol', how='outer').merge(rwconcat6Kp_8Kp, on='Symbol', how='outer').merge(rwconcat6Kp_8RA, on='Symbol', how='outer').merge(rwconcat6RA_8Kp, on='Symbol',
how='outer').merge(rwconcat6RA_8RA, on='Symbol', how='outer').merge(rwconcat4Kp_7Kp, on='Symbol', how='outer').merge(rwconcat4Kp_7RA, on='Symbol',
how='outer').merge(rwconcat4RA_7Kp, on='Symbol', how='outer').merge(rwconcat4RA_7RA, on='Symbol', how='outer').merge(rwconcat5Kp_7Kp, on='Symbol',
how='outer').merge(rwconcat5Kp_7RA, on='Symbol', how='outer').merge(rwconcat5RA_7Kp, on='Symbol', how='outer').merge(rwconcat5RA_7RA, on='Symbol', how='outer')
#concat3456 = pd.concat([, , , , , , , ], axis=1, join='outer').drop_duplicates(keep='first')
Sharedlist = concat345678.to_csv(output + "SC345678_2ormorescreens_all_depleted" + sfdr + ".csv", encoding='utf-8')
rsigs345678 = pd.read_csv(output + "SC345678_2ormorescreens_all_depleted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#concat3456_old = pd.concat([rwconcat3Kp_6Kp, rwconcat3Kp_6RA, rwconcat3RA_6Kp, rwconcat3RA_6RA], axis=1, join='outer').drop_duplicates(keep='first')
concat345678_old = rwconcat3Kp_6Kp.merge(rwconcat3Kp_6RA, on='Symbol', how='outer').merge(rwconcat3RA_6Kp, on='Symbol', how='outer').merge(rwconcat3RA_6RA, on='Symbol', how='outer').merge(rwconcat3Kp_8Kp, on='Symbol', how='outer').merge(rwconcat3Kp_8RA, on='Symbol',
how='outer').merge(rwconcat3RA_8Kp, on='Symbol', how='outer').merge(rwconcat3RA_8RA, on='Symbol', how='outer').merge(rwconcat6Kp_8Kp, on='Symbol', how='outer').merge(rwconcat6Kp_8RA, on='Symbol', how='outer').merge(rwconcat6RA_8Kp, on='Symbol',
how='outer').merge(rwconcat6RA_8RA, on='Symbol', how='outer')
Sharedlist_old = concat345678_old.to_csv(output + "SC345678_2ormorescreens_all_depleted_OLD_" + sfdr + ".csv", encoding='utf-8')
rsigs345678_old = pd.read_csv(output + "SC345678_2ormorescreens_all_depleted_OLD_" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#concat3456_young = pd.concat([rwconcat4Kp_5Kp, rwconcat4Kp_5RA, rwconcat4RA_5Kp, rwconcat4RA_5RA], axis=1, join='outer').drop_duplicates(keep='first')
concat345678_young = rwconcat4Kp_5Kp.merge(rwconcat4Kp_5RA, on='Symbol', how='outer').merge(rwconcat4RA_5Kp, on='Symbol', how='outer').merge(rwconcat4RA_5RA, on='Symbol', how='outer').merge(rwconcat4Kp_7Kp, on='Symbol', how='outer').merge(rwconcat4Kp_7RA, on='Symbol',
how='outer').merge(rwconcat4RA_7Kp, on='Symbol', how='outer').merge(rwconcat4RA_7RA, on='Symbol', how='outer').merge(rwconcat5Kp_7Kp, on='Symbol', how='outer').merge(rwconcat5Kp_7RA, on='Symbol', how='outer').merge(rwconcat5RA_7Kp, on='Symbol',
how='outer').merge(rwconcat5RA_7RA, on='Symbol', how='outer')
Sharedlist_young = concat345678_young.to_csv(output + "SC345678_2ormorescreens_all_depleted_YOUNG_" + sfdr + ".csv", encoding='utf-8')
rsigs345678_young = pd.read_csv(output + "SC345678_2ormorescreens_all_depleted_YOUNG_" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

##make final lists by forming a union of young and old, and then joing union with young or old as outer merge
#make union
union_list = pd.merge(rsigs345678_young, rsigs345678_old, on='Symbol', how='inner')
writefinal_union = union_list.to_csv(output + "SC345678_Union_depleted" + sfdr + ".csv", encoding='utf-8')
readfinal_union = pd.read_csv(output + "SC345678_Union_depleted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

finallist_old = pd.merge(readfinal_union, rsigs345678_old, on='Symbol', how='outer', indicator=True).query('_merge=="right_only"')
writefinal_old = finallist_old.to_csv(output + "SC345678_Oldspecific_all_depleted" + sfdr + ".csv", encoding='utf-8')
readfinal_old = pd.read_csv(output + "SC345678_Oldspecific_all_depleted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

finallist_young = pd.merge(readfinal_union, rsigs345678_young, on='Symbol', how='outer', indicator=True).query('_merge=="right_only"')
writefinal_young = finallist_young.to_csv(output + "SC345678_Youngspecific_all_depleted" + sfdr + ".csv", encoding='utf-8')
readfinal_young = pd.read_csv(output + "SC345678_Youngspecific_all_depleted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
