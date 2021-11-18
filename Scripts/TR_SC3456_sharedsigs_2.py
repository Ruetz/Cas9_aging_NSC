##Make a list of shared sig hits from SC3456
#Start in cd Desktop/Dropbox/TR_Cas9Paper_CodeCheck
#PATH=$PATH/usr/local/bin:/usr/local/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/Users/tysonruetz/Desktop/Dropbox/TR_Cas9Paper_CodeCheck_BACKUP/bowtie

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame, Series

#specify fdr
fdr = float(0.1)
sfdr = str("_fdr" + str(fdr))

#specify output folder directory
output = "Results/SC3456_sharedsigs/"


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


#Capture file names for columns
stdata1 = str(data1[data1.find('SC'):data1.find('F')])
stdata2 = str(data2[data1.find('SC'):data2.find('F')])
stdata3 = str(data3[data1.find('SC'):data3.find('F')])
stdata4 = str(data4[data1.find('SC'):data4.find('F')])
stdata5 = str(data5[data1.find('SC'):data5.find('F')])
stdata6 = str(data6[data1.find('SC'):data6.find('F')])
stdata7 = str(data7[data1.find('SC'):data7.find('F')])
stdata8 = str(data8[data1.find('SC'):data8.find('F')])

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


#narrow to only FDR_X significant gene list
d1pval = rdata1[rdata1['FDR_0.1'] <= fdr]
d2pval = rdata2[rdata2['FDR_0.1'] <= fdr]
d3pval = rdata3[rdata3['FDR_0.1'] <= fdr]
d4pval = rdata4[rdata4['FDR_0.1'] <= fdr]
d5pval = rdata5[rdata5['FDR_0.1'] <= fdr]
d6pval = rdata6[rdata6['FDR_0.1'] <= fdr]
d7pval = rdata7[rdata7['FDR_0.1'] <= fdr]
d8pval = rdata8[rdata8['FDR_0.1'] <= fdr]


#Write to csv file
sigsd1 = d1pval.to_csv(output + stdata1 + sfdr + ".csv")
sigsd2 = d2pval.to_csv(output + stdata2 + sfdr + ".csv")
sigsd3 = d3pval.to_csv(output + stdata3 + sfdr + ".csv")
sigsd4 = d4pval.to_csv(output + stdata4 + sfdr + ".csv")
sigsd5 = d5pval.to_csv(output + stdata5 + sfdr + ".csv")
sigsd6 = d6pval.to_csv(output + stdata6 + sfdr + ".csv")
sigsd7 = d7pval.to_csv(output + stdata7 + sfdr + ".csv")
sigsd8 = d8pval.to_csv(output + stdata8 + sfdr + ".csv")


##first read csv files
rsigsd1 = pd.read_csv(output + stdata1 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd2 = pd.read_csv(output + stdata2 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd3 = pd.read_csv(output + stdata3 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd4 = pd.read_csv(output + stdata4 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd5 = pd.read_csv(output + stdata5 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd6 = pd.read_csv(output + stdata6 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd7 = pd.read_csv(output + stdata7 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd8 = pd.read_csv(output + stdata8 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

##########################################################
#Making enriched list
#narrow to only enriched values
rsigsd1_enr = rsigsd1[rsigsd1['casTLE Effect'] > 0]
rsigsd2_enr = rsigsd2[rsigsd2['casTLE Effect'] > 0]
rsigsd3_enr = rsigsd3[rsigsd3['casTLE Effect'] > 0]
rsigsd4_enr = rsigsd4[rsigsd4['casTLE Effect'] > 0]
rsigsd5_enr = rsigsd5[rsigsd5['casTLE Effect'] > 0]
rsigsd6_enr = rsigsd6[rsigsd6['casTLE Effect'] > 0]
rsigsd7_enr = rsigsd7[rsigsd7['casTLE Effect'] > 0]
rsigsd8_enr = rsigsd8[rsigsd8['casTLE Effect'] > 0]

###then concat
concat3Kp_6Kp = pd.merge(rsigsd1_enr, rsigsd7_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3Kp_6RA = pd.merge(rsigsd1_enr, rsigsd8_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_6Kp = pd.merge(rsigsd2_enr, rsigsd7_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_6RA = pd.merge(rsigsd2_enr, rsigsd8_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4Kp_5Kp = pd.merge(rsigsd3_enr, rsigsd5_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4Kp_5RA = pd.merge(rsigsd3_enr, rsigsd6_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_5Kp = pd.merge(rsigsd4_enr, rsigsd5_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_5RA = pd.merge(rsigsd4_enr, rsigsd6_enr, on='Symbol', how='inner').drop_duplicates(keep='first')


##write individual concats to file
wconcat3Kp_6Kp = concat3Kp_6Kp.to_csv(output + stdata1 + stdata7 + sfdr + ".csv")
wconcat3Kp_6RA = concat3Kp_6RA.to_csv(output + stdata1 + stdata8 + sfdr + ".csv")
wconcat3RA_6Kp = concat3RA_6Kp.to_csv(output + stdata2 + stdata7 + sfdr + ".csv")
wconcat3RA_6RA = concat3RA_6RA.to_csv(output + stdata2 + stdata8 + sfdr + ".csv")
wconcat4Kp_5Kp = concat4Kp_5Kp.to_csv(output + stdata3 + stdata5 + sfdr + ".csv")
wconcat4Kp_5RA = concat4Kp_5RA.to_csv(output + stdata3 + stdata6 + sfdr + ".csv")
wconcat4RA_5Kp = concat4RA_5Kp.to_csv(output + stdata4 + stdata5 + sfdr + ".csv")
wconcat4RA_5RA = concat4RA_5RA.to_csv(output + stdata4 + stdata6 + sfdr + ".csv")

##read files
rwconcat3Kp_6Kp = pd.read_csv(output + stdata1 + stdata7 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3Kp_6RA = pd.read_csv(output + stdata1 + stdata8 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_6Kp = pd.read_csv(output + stdata2 + stdata7 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_6RA = pd.read_csv(output + stdata2 + stdata8 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4Kp_5Kp = pd.read_csv(output + stdata3 + stdata5 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4Kp_5RA = pd.read_csv(output + stdata3 + stdata6 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_5Kp = pd.read_csv(output + stdata4 + stdata5 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_5RA = pd.read_csv(output + stdata4 + stdata6 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')


#make a master sigs files
#df1.merge(df2,on='name').merge(df3,on='name')
concat3456 = rwconcat3Kp_6Kp.merge(rwconcat3Kp_6RA, on='Symbol', how='outer').merge(rwconcat3RA_6Kp, on='Symbol', how='outer').merge(rwconcat3RA_6RA,
on='Symbol', how='outer').merge(rwconcat4Kp_5Kp, on='Symbol', how='outer').merge(rwconcat4Kp_5RA, on='Symbol', how='outer').merge(rwconcat4RA_5Kp, on='Symbol', how='outer').merge(rwconcat4RA_5RA, on='Symbol', how='outer')
#concat3456 = pd.concat([, , , , , , , ], axis=1, join='outer').drop_duplicates(keep='first')
Sharedlist = concat3456.to_csv(output + "SC3456_2ormorescreens_all_enriched" + sfdr + ".csv", encoding='utf-8')
rsigs3456 = pd.read_csv(output + "SC3456_2ormorescreens_all_enriched" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#concat3456_old = pd.concat([rwconcat3Kp_6Kp, rwconcat3Kp_6RA, rwconcat3RA_6Kp, rwconcat3RA_6RA], axis=1, join='outer').drop_duplicates(keep='first')
concat3456_old = rwconcat3Kp_6Kp.merge(rwconcat3Kp_6RA, on='Symbol', how='outer').merge(rwconcat3RA_6Kp, on='Symbol', how='outer').merge(rwconcat3RA_6RA, on='Symbol', how='outer')
Sharedlist_old = concat3456_old.to_csv(output + "SC3465_2ormorescreens_all_enriched_OLD_" + sfdr + ".csv", encoding='utf-8')
rsigs3456_old = pd.read_csv(output + "SC3465_2ormorescreens_all_enriched_OLD_" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#concat3456_young = pd.concat([rwconcat4Kp_5Kp, rwconcat4Kp_5RA, rwconcat4RA_5Kp, rwconcat4RA_5RA], axis=1, join='outer').drop_duplicates(keep='first')
concat3456_young = rwconcat4Kp_5Kp.merge(rwconcat4Kp_5RA, on='Symbol', how='outer').merge(rwconcat4RA_5Kp, on='Symbol', how='outer').merge(rwconcat4RA_5RA, on='Symbol', how='outer')
Sharedlist_young = concat3456_young.to_csv(output + "SC34_65_2ormorescreens_all_enriched_YOUNG_" + sfdr + ".csv", encoding='utf-8')
rsigs3456_young = pd.read_csv(output + "SC34_65_2ormorescreens_all_enriched_YOUNG_" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

##make final lists by forming a union of young and old, and then joing union with young or old as outer merge
#make union
union_list = pd.merge(rsigs3456_young, rsigs3456_old, on='Symbol', how='inner')
writefinal_union = union_list.to_csv(output + "SC3456_Union_enriched" + sfdr + ".csv", encoding='utf-8')
readfinal_union = pd.read_csv(output + "SC3456_Union_enriched" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

finallist_old = pd.merge(readfinal_union, rsigs3456_old, on='Symbol', how='outer', indicator=True).query('_merge=="right_only"')
writefinal_old = finallist_old.to_csv(output + "SC3456_Oldspecific_all_enriched" + sfdr + ".csv", encoding='utf-8')
readfinal_old = pd.read_csv(output + "SC3456_Oldspecific_all_enriched" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

finallist_young = pd.merge(readfinal_union, rsigs3456_young, on='Symbol', how='outer', indicator=True).query('_merge=="right_only"')
writefinal_young = finallist_young.to_csv(output + "SC3456_Youngspecific_all_enriched" + sfdr + ".csv", encoding='utf-8')
readfinal_young = pd.read_csv(output + "SC3456_Youngspecific_all_enriched" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

##########################################################
#Making depleted list
rsigsd1_enr = rsigsd1[rsigsd1['casTLE Effect'] < 0]
rsigsd2_enr = rsigsd2[rsigsd2['casTLE Effect'] < 0]
rsigsd3_enr = rsigsd3[rsigsd3['casTLE Effect'] < 0]
rsigsd4_enr = rsigsd4[rsigsd4['casTLE Effect'] < 0]
rsigsd5_enr = rsigsd5[rsigsd5['casTLE Effect'] < 0]
rsigsd6_enr = rsigsd6[rsigsd6['casTLE Effect'] < 0]
rsigsd7_enr = rsigsd7[rsigsd7['casTLE Effect'] < 0]
rsigsd8_enr = rsigsd8[rsigsd8['casTLE Effect'] < 0]

###then concat
concat3Kp_6Kp = pd.merge(rsigsd1_enr, rsigsd7_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3Kp_6RA = pd.merge(rsigsd1_enr, rsigsd8_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_6Kp = pd.merge(rsigsd2_enr, rsigsd7_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat3RA_6RA = pd.merge(rsigsd2_enr, rsigsd8_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4Kp_5Kp = pd.merge(rsigsd3_enr, rsigsd5_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4Kp_5RA = pd.merge(rsigsd3_enr, rsigsd6_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_5Kp = pd.merge(rsigsd4_enr, rsigsd5_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
concat4RA_5RA = pd.merge(rsigsd4_enr, rsigsd6_enr, on='Symbol', how='inner').drop_duplicates(keep='first')


##write individual concats to file
wconcat3Kp_6Kp = concat3Kp_6Kp.to_csv(output + stdata1 + stdata7 + sfdr + ".csv")
wconcat3Kp_6RA = concat3Kp_6RA.to_csv(output + stdata1 + stdata8 + sfdr + ".csv")
wconcat3RA_6Kp = concat3RA_6Kp.to_csv(output + stdata2 + stdata7 + sfdr + ".csv")
wconcat3RA_6RA = concat3RA_6RA.to_csv(output + stdata2 + stdata8 + sfdr + ".csv")
wconcat4Kp_5Kp = concat4Kp_5Kp.to_csv(output + stdata3 + stdata5 + sfdr + ".csv")
wconcat4Kp_5RA = concat4Kp_5RA.to_csv(output + stdata3 + stdata6 + sfdr + ".csv")
wconcat4RA_5Kp = concat4RA_5Kp.to_csv(output + stdata4 + stdata5 + sfdr + ".csv")
wconcat4RA_5RA = concat4RA_5RA.to_csv(output + stdata4 + stdata6 + sfdr + ".csv")

##read files
rwconcat3Kp_6Kp = pd.read_csv(output + stdata1 + stdata7 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3Kp_6RA = pd.read_csv(output + stdata1 + stdata8 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_6Kp = pd.read_csv(output + stdata2 + stdata7 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat3RA_6RA = pd.read_csv(output + stdata2 + stdata8 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4Kp_5Kp = pd.read_csv(output + stdata3 + stdata5 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4Kp_5RA = pd.read_csv(output + stdata3 + stdata6 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_5Kp = pd.read_csv(output + stdata4 + stdata5 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rwconcat4RA_5RA = pd.read_csv(output + stdata4 + stdata6 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')


#make a master sigs files
#df1.merge(df2,on='name').merge(df3,on='name')
concat3456 = rwconcat3Kp_6Kp.merge(rwconcat3Kp_6RA, on='Symbol', how='outer').merge(rwconcat3RA_6Kp, on='Symbol', how='outer').merge(rwconcat3RA_6RA,
on='Symbol', how='outer').merge(rwconcat4Kp_5Kp, on='Symbol', how='outer').merge(rwconcat4Kp_5RA, on='Symbol', how='outer').merge(rwconcat4RA_5Kp, on='Symbol', how='outer').merge(rwconcat4RA_5RA, on='Symbol', how='outer')
#concat3456 = pd.concat([, , , , , , , ], axis=1, join='outer').drop_duplicates(keep='first')
Sharedlist = concat3456.to_csv(output + "SC3456_2ormorescreens_all_depleted" + sfdr + ".csv", encoding='utf-8')
rsigs3456 = pd.read_csv(output + "SC3456_2ormorescreens_all_depleted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#concat3456_old = pd.concat([rwconcat3Kp_6Kp, rwconcat3Kp_6RA, rwconcat3RA_6Kp, rwconcat3RA_6RA], axis=1, join='outer').drop_duplicates(keep='first')
concat3456_old = rwconcat3Kp_6Kp.merge(rwconcat3Kp_6RA, on='Symbol', how='outer').merge(rwconcat3RA_6Kp, on='Symbol', how='outer').merge(rwconcat3RA_6RA, on='Symbol', how='outer')
Sharedlist_old = concat3456_old.to_csv(output + "SC3465_2ormorescreens_all_depleted_OLD_" + sfdr + ".csv", encoding='utf-8')
rsigs3456_old = pd.read_csv(output + "SC3465_2ormorescreens_all_depleted_OLD_" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#concat3456_young = pd.concat([rwconcat4Kp_5Kp, rwconcat4Kp_5RA, rwconcat4RA_5Kp, rwconcat4RA_5RA], axis=1, join='outer').drop_duplicates(keep='first')
concat3456_young = rwconcat4Kp_5Kp.merge(rwconcat4Kp_5RA, on='Symbol', how='outer').merge(rwconcat4RA_5Kp, on='Symbol', how='outer').merge(rwconcat4RA_5RA, on='Symbol', how='outer')
Sharedlist_young = concat3456_young.to_csv(output + "SC34_65_2ormorescreens_all_depleted_YOUNG_" + sfdr + ".csv", encoding='utf-8')
rsigs3456_young = pd.read_csv(output + "SC34_65_2ormorescreens_all_depleted_YOUNG_" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

##make final lists by forming a union of young and old, and then joing union with young or old as outer merge
union_list = pd.merge(rsigs3456_young, rsigs3456_old, on='Symbol', how='inner')
writefinal_union = union_list.to_csv(output + "SC3456_Union_depleted" + sfdr + ".csv", encoding='utf-8')
readfinal_union = pd.read_csv(output + "SC3456_Union_depleted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

finallist_old = pd.merge(readfinal_union, rsigs3456_old, on='Symbol', how='outer', indicator=True).query('_merge=="right_only"')
writefinal_old = finallist_old.to_csv(output + "SC3456_Oldspecific_all_depleted" + sfdr + ".csv", encoding='utf-8')
readfinal_old = pd.read_csv(output + "SC3456_Oldspecific_all_depleted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

finallist_young = pd.merge(readfinal_union, rsigs3456_young, on='Symbol', how='outer', indicator=True).query('_merge=="right_only"')
writefinal_young = finallist_young.to_csv(output + "SC3456_Youngspecific_all_depleted" + sfdr + ".csv", encoding='utf-8')
readfinal_young = pd.read_csv(output + "SC3456_Youngspecific_all_depleted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#######################
# make a unique depleted list by removing any gene depleted in SC3456 Q libraries
# Make Q libraries
import subprocess
subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC4_Y_Q_FKDL190724827-1a-31_HWHG3CCXY_L2_1.fq.gz SC4_Y_Q mm-Cas9-10,', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC5_Y_Q_FKDL190724827-1a-3_HWHG3CCXY_L2_1.fq.gz SC5_Y_Q mm-Cas9-10,', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC36_O_Q_FKDL190724827-1a-8_HWHG3CCXY_L2_1.fq.gz SC36_O_Q mm-Cas9-10,', shell=True)

subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC4_Y_Q_counts.csv SC4_QvLib', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC5_Y_Q_counts.csv SC5_QvLib', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC36_O_Q_counts.csv SC36_QvLib', shell=True)
subprocess.call('python Scripts/addPermutations.py Results/SC4_QvLib.csv 100000', shell=True)
subprocess.call('python Scripts/addPermutations.py Results/SC5_QvLib.csv 100000', shell=True)
subprocess.call('python Scripts/addPermutations.py Results/SC36_QvLib.csv 100000', shell=True)

##Q library FDR adjust
import scipy
import statsmodels.stats.multitest as smm
import csv

fdr = float('0.1')
sfdr = str(fdr)
fdrmeth = 'fdr_bh'
sfdrmeth = str(fdrmeth)

data10 = 'Results/SC4_QvLib.csv'
data11 = 'Results/SC5_QvLib.csv'
data12 = 'Results/SC36_QvLib.csv'

sdata10 = str(data10[0:data10.find('.')])
sdata11 = str(data11[0:data11.find('.')])
sdata12 = str(data12[0:data12.find('.')])

rdata10 = pd.read_csv(data10, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata11 = pd.read_csv(data11, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata12 = pd.read_csv(data12, encoding='utf-8', header=None, skiprows=1, index_col='Symbol', names=['#GeneID', 'Symbol', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'Minimum Effect Estimate', 'Maximum Effect Estimate', 'Individual Elements',
'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
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

##Q library depleted list
data10 = 'Results/SC4_QvLib_FDR0.1.csv'
data11 = 'Results/SC5_QvLib_FDR0.1.csv'
data12 = 'Results/SC36_QvLib_FDR0.1.csv'
stdata10 = str(data10[data10.find('SC'):data10.find('F')])
stdata11 = str(data11[data11.find('SC'):data11.find('F')])
stdata12 = str(data12[data12.find('SC'):data12.find('F')])
rdata10 = pd.read_csv(data10, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata11 = pd.read_csv(data11, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata12 = pd.read_csv(data12, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate',
'Maximum Effect Estimate', 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
d10pval = rdata10[rdata10['FDR_0.1'] <= fdr]
d11pval = rdata11[rdata11['FDR_0.1'] <= fdr]
d12pval = rdata12[rdata12['FDR_0.1'] <= fdr]
sigsd10 = d10pval.to_csv(output + stdata10 + sfdr + ".csv")
sigsd11 = d11pval.to_csv(output + stdata11 + sfdr + ".csv")
sigsd12 = d12pval.to_csv(output + stdata12 + sfdr + ".csv")
rsigsd10 = pd.read_csv(output + stdata10 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd11 = pd.read_csv(output + stdata11 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd12 = pd.read_csv(output + stdata12 + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
rsigsd10_enr = rsigsd10[rsigsd10['casTLE Effect'] < 0]
rsigsd11_enr = rsigsd11[rsigsd11['casTLE Effect'] < 0]
rsigsd12_enr = rsigsd12[rsigsd12['casTLE Effect'] < 0]

#readfinal_old = pd.read_csv(output + "SC3456_Oldspecific_all_depleted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#make a master Q sigs files
concat3456Q = rsigsd10_enr.merge(rsigsd11_enr, on='Symbol', how='outer').merge(rsigsd12_enr,on='Symbol', how='outer')
wconcat3456Q = concat3456Q.to_csv(output + "SC3456_Q_all_depleted" + sfdr + ".csv", encoding='utf-8')
rwsigs3456Q = pd.read_csv(output + "SC3456_Q_all_depleted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
##make final lists for Old Screen depleted with Q depleted subtracted
finallist_oldQ = pd.merge(rwsigs3456Q, readfinal_old, on='Symbol', how='outer', indicator='merge2').query('merge2=="right_only"')
writefinal_oldQ = finallist_oldQ.to_csv(output + "SC3456_Oldspecific_all_depleted_Q_subtracted" + sfdr + ".csv", encoding='utf-8')
readfinal_oldQ = pd.read_csv(output + "SC3456_Oldspecific_all_depleted_Q_subtracted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')
##make final lists for Young Screen depleted with Q depleted subtracted
finallist_YoungQ = pd.merge(rwsigs3456Q, readfinal_young, on='Symbol', how='outer', indicator='merge2').query('merge2=="right_only"')
writefinal_youngQ = finallist_YoungQ.to_csv(output + "SC3456_Youngspecific_all_depleted_Q_subtracted" + sfdr + ".csv", encoding='utf-8')
readfinal_YoungQ = pd.read_csv(output + "SC3456_Youngspecific_all_depleted_Q_subtracted" + sfdr + ".csv", index_col='Symbol', encoding='utf-8')

#############################
#######################
#make a unique depleted list using only CasTLE p-values by removing any gene depleted in SC3456 Q libraries
#Make Q libraries
#adjust p-value as needed
# pval = float('0.05')
# spval = str(pval)
#
# #Narrow to pval significant hits
# d1pval = rdata1[rdata1['casTLE p-value'] <= pval]
# d2pval = rdata2[rdata2['casTLE p-value'] <= pval]
# d3pval = rdata3[rdata3['casTLE p-value'] <= pval]
# d4pval = rdata4[rdata4['casTLE p-value'] <= pval]
# d5pval = rdata5[rdata5['casTLE p-value'] <= pval]
# d6pval = rdata6[rdata6['casTLE p-value'] <= pval]
# d7pval = rdata7[rdata7['casTLE p-value'] <= pval]
# d8pval = rdata8[rdata8['casTLE p-value'] <= pval]
#
# #Write to csv file
# sigsd1 = d1pval.to_csv(output + stdata1 + spval + ".csv")
# sigsd2 = d2pval.to_csv(output + stdata2 + spval + ".csv")
# sigsd3 = d3pval.to_csv(output + stdata3 + spval + ".csv")
# sigsd4 = d4pval.to_csv(output + stdata4 + spval + ".csv")
# sigsd5 = d5pval.to_csv(output + stdata5 + spval + ".csv")
# sigsd6 = d6pval.to_csv(output + stdata6 + spval + ".csv")
# sigsd7 = d7pval.to_csv(output + stdata7 + spval + ".csv")
# sigsd8 = d8pval.to_csv(output + stdata8 + spval + ".csv")
#
# ## read csv files
# rsigsd1 = pd.read_csv(output + stdata1 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rsigsd2 = pd.read_csv(output + stdata2 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rsigsd3 = pd.read_csv(output + stdata3 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rsigsd4 = pd.read_csv(output + stdata4 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rsigsd5 = pd.read_csv(output + stdata5 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rsigsd6 = pd.read_csv(output + stdata6 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rsigsd7 = pd.read_csv(output + stdata7 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rsigsd8 = pd.read_csv(output + stdata8 + spval + ".csv", index_col='Symbol', encoding='utf-8')
#
# #Making enriched list
# #narrow to only enriched values
# rsigsd1_enr = rsigsd1[rsigsd1['casTLE Effect'] < 0]
# rsigsd2_enr = rsigsd2[rsigsd2['casTLE Effect'] < 0]
# rsigsd3_enr = rsigsd3[rsigsd3['casTLE Effect'] < 0]
# rsigsd4_enr = rsigsd4[rsigsd4['casTLE Effect'] < 0]
# rsigsd5_enr = rsigsd5[rsigsd5['casTLE Effect'] < 0]
# rsigsd6_enr = rsigsd6[rsigsd6['casTLE Effect'] < 0]
# rsigsd7_enr = rsigsd7[rsigsd7['casTLE Effect'] < 0]
# rsigsd8_enr = rsigsd8[rsigsd8['casTLE Effect'] < 0]
#
# ###then concat
# concat3Kp_6Kp = pd.merge(rsigsd1_enr, rsigsd7_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
# concat3Kp_6RA = pd.merge(rsigsd1_enr, rsigsd8_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
# concat3RA_6Kp = pd.merge(rsigsd2_enr, rsigsd7_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
# concat3RA_6RA = pd.merge(rsigsd2_enr, rsigsd8_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
# concat4Kp_5Kp = pd.merge(rsigsd3_enr, rsigsd5_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
# concat4Kp_5RA = pd.merge(rsigsd3_enr, rsigsd6_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
# concat4RA_5Kp = pd.merge(rsigsd4_enr, rsigsd5_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
# concat4RA_5RA = pd.merge(rsigsd4_enr, rsigsd6_enr, on='Symbol', how='inner').drop_duplicates(keep='first')
#
#
# ##write individual concats to file
# wconcat3Kp_6Kp = concat3Kp_6Kp.to_csv(output + stdata1 + stdata7 + spval + ".csv")
# wconcat3Kp_6RA = concat3Kp_6RA.to_csv(output + stdata1 + stdata8 + spval + ".csv")
# wconcat3RA_6Kp = concat3RA_6Kp.to_csv(output + stdata2 + stdata7 + spval + ".csv")
# wconcat3RA_6RA = concat3RA_6RA.to_csv(output + stdata2 + stdata8 + spval + ".csv")
# wconcat4Kp_5Kp = concat4Kp_5Kp.to_csv(output + stdata3 + stdata5 + spval + ".csv")
# wconcat4Kp_5RA = concat4Kp_5RA.to_csv(output + stdata3 + stdata6 + spval + ".csv")
# wconcat4RA_5Kp = concat4RA_5Kp.to_csv(output + stdata4 + stdata5 + spval + ".csv")
# wconcat4RA_5RA = concat4RA_5RA.to_csv(output + stdata4 + stdata6 + spval + ".csv")
#
# ##read files
# rwconcat3Kp_6Kp = pd.read_csv(output + stdata1 + stdata7 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rwconcat3Kp_6RA = pd.read_csv(output + stdata1 + stdata8 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rwconcat3RA_6Kp = pd.read_csv(output + stdata2 + stdata7 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rwconcat3RA_6RA = pd.read_csv(output + stdata2 + stdata8 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rwconcat4Kp_5Kp = pd.read_csv(output + stdata3 + stdata5 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rwconcat4Kp_5RA = pd.read_csv(output + stdata3 + stdata6 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rwconcat4RA_5Kp = pd.read_csv(output + stdata4 + stdata5 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rwconcat4RA_5RA = pd.read_csv(output + stdata4 + stdata6 + spval + ".csv", index_col='Symbol', encoding='utf-8')
#
# #make a master sigs files
# #df1.merge(df2,on='name').merge(df3,on='name')
# concat3456 = rwconcat3Kp_6Kp.merge(rwconcat3Kp_6RA, on='Symbol', how='outer').merge(rwconcat3RA_6Kp, on='Symbol', how='outer').merge(rwconcat3RA_6RA,
# on='Symbol', how='outer').merge(rwconcat4Kp_5Kp, on='Symbol', how='outer').merge(rwconcat4Kp_5RA, on='Symbol', how='outer').merge(rwconcat4RA_5Kp, on='Symbol', how='outer').merge(rwconcat4RA_5RA, on='Symbol', how='outer')
# #concat3456 = pd.concat([, , , , , , , ], axis=1, join='outer').drop_duplicates(keep='first')
# Sharedlist = concat3456.to_csv(output + "SC3456_2ormorescreens_all_depleted_PVAL" + spval + ".csv", encoding='utf-8')
# rsigs3456 = pd.read_csv(output + "SC3456_2ormorescreens_all_depleted_PVAL" + spval + ".csv", index_col='Symbol', encoding='utf-8')
#
# #concat3456_old = pd.concat([rwconcat3Kp_6Kp, rwconcat3Kp_6RA, rwconcat3RA_6Kp, rwconcat3RA_6RA], axis=1, join='outer').drop_duplicates(keep='first')
# concat3456_old = rwconcat3Kp_6Kp.merge(rwconcat3Kp_6RA, on='Symbol', how='outer').merge(rwconcat3RA_6Kp, on='Symbol', how='outer').merge(rwconcat3RA_6RA, on='Symbol', how='outer')
# Sharedlist_old = concat3456_old.to_csv(output + "SC3465_2ormorescreens_all_depleted_OLD_PVAL" + spval + ".csv", encoding='utf-8')
# rsigs3456_old = pd.read_csv(output + "SC3465_2ormorescreens_all_depleted_OLD_PVAL" + spval + ".csv", index_col='Symbol', encoding='utf-8')
#
# #concat3456_young = pd.concat([rwconcat4Kp_5Kp, rwconcat4Kp_5RA, rwconcat4RA_5Kp, rwconcat4RA_5RA], axis=1, join='outer').drop_duplicates(keep='first')
# concat3456_young = rwconcat4Kp_5Kp.merge(rwconcat4Kp_5RA, on='Symbol', how='outer').merge(rwconcat4RA_5Kp, on='Symbol', how='outer').merge(rwconcat4RA_5RA, on='Symbol', how='outer')
# Sharedlist_young = concat3456_young.to_csv(output + "SC34_65_2ormorescreens_all_depleted_YOUNG_PVAL" + spval + ".csv", encoding='utf-8')
# rsigs3456_young = pd.read_csv(output + "SC34_65_2ormorescreens_all_depleted_YOUNG_PVAL" + spval + ".csv", index_col='Symbol', encoding='utf-8')
#
# ##make final lists by forming a union of young and old, and then joing union with young or old as outer merge
# union_list = pd.merge(rsigs3456_young, rsigs3456_old, on='Symbol', how='inner')
# writefinal_union = union_list.to_csv(output + "SC3456_Union_depleted_PVAL" + spval + ".csv", encoding='utf-8')
# readfinal_union = pd.read_csv(output + "SC3456_Union_depleted_PVAL" + spval + ".csv", index_col='Symbol', encoding='utf-8')
#
# finallist_old = pd.merge(readfinal_union, rsigs3456_old, on='Symbol', how='outer', indicator=True).query('_merge=="right_only"')
# writefinal_old = finallist_old.to_csv(output + "SC3456_Oldspecific_all_depleted_PVAL" + spval + ".csv", encoding='utf-8')
# readfinal_old = pd.read_csv(output + "SC3456_Oldspecific_all_depleted_PVAL" + spval + ".csv", index_col='Symbol', encoding='utf-8')
#
# finallist_young = pd.merge(readfinal_union, rsigs3456_young, on='Symbol', how='outer', indicator=True).query('_merge=="right_only"')
# writefinal_young = finallist_young.to_csv(output + "SC3456_Youngspecific_all_depleted_PVAL" + spval + ".csv", encoding='utf-8')
# readfinal_young = pd.read_csv(output + "SC3456_Youngspecific_all_depleted_PVAL" + spval + ".csv", index_col='Symbol', encoding='utf-8')
#
#
# ##############################################################################################################
# ##Q library depleted list
#
# sigsd10 = d10pval.to_csv(output + stdata10 + spval + ".csv")
# sigsd11 = d11pval.to_csv(output + stdata11 + spval + ".csv")
# sigsd12 = d12pval.to_csv(output + stdata12 + spval + ".csv")
# rsigsd10 = pd.read_csv(output + stdata10 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rsigsd11 = pd.read_csv(output + stdata11 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rsigsd12 = pd.read_csv(output + stdata12 + spval + ".csv", index_col='Symbol', encoding='utf-8')
# rsigsd10_enr = rsigsd10[rsigsd10['casTLE Effect'] < 0]
# rsigsd11_enr = rsigsd11[rsigsd11['casTLE Effect'] < 0]
# rsigsd12_enr = rsigsd12[rsigsd12['casTLE Effect'] < 0]
#
# #make a master Q sigs files
# concat3456Q = rsigsd10_enr.merge(rsigsd11_enr, on='Symbol', how='outer').merge(rsigsd12_enr,on='Symbol', how='outer')
# Sharedlist = concat3456Q.to_csv(output + "SC3456_Q_all_depleted_pval" + spval + ".csv", encoding='utf-8')
# rsigs3456Q = pd.read_csv(output + "SC3456_Q_all_depleted_pval" + spval + ".csv", index_col='Symbol', encoding='utf-8')
# ##make final lists for Screen depleted with Q depleted subtracted
# finallist_oldQ = pd.merge(rsigs3456Q, readfinal_old, on='Symbol', how='outer', indicator='merge2').query('merge2=="right_only"')
# writefinal_oldQ = finallist_oldQ.to_csv(output + "SC3456_Oldspecific_all_depleted_Q_subtracted_pval" + spval + ".csv", encoding='utf-8')
# readfinal_oldQ = pd.read_csv(output + "SC3456_Oldspecific_all_depleted_Q_subtracted_pval" + spval + ".csv", index_col='Symbol', encoding='utf-8')
