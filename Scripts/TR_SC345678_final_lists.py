#Make final lists of genes with gene scores from all screens

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame, Series

output = str('Results/Final_gene_lists/')

list1 = 'Results/SC345678_sharedsigs/SC345678_Oldspecific_all_enriched_fdr0.1.csv'
list2 = 'Results/PCA/SC345678_allgenescores_d4.csv'
list3 = 'Results/PCA/SC345678_allgenescores_d14.csv'
list4 = 'Results/SC345678_sharedsigs/SC345678_Youngspecific_all_enriched_fdr0.1.csv'
list5 = 'Results/SC345678_sharedsigs/SC345678_Union_enriched_fdr0.1.csv'
list6 = 'Results/SC345678_sharedsigs/SC345678_Oldspecific_all_depleted_fdr0.1.csv'
list7 = 'Results/SC345678_sharedsigs/SC345678_Youngspecific_all_depleted_fdr0.1.csv'
list8 = 'Results/SC345678_sharedsigs/SC345678_Union_depleted_fdr0.1.csv'
list9 = 'Results/SC3456_sharedsigs/SC3456_Oldspecific_all_depleted_Q_subtracted0.1.csv'
list10 = 'Results/SC3456_sharedsigs/SC3456_Youngspecific_all_depleted_Q_subtracted0.1.csv'
list11 = 'Results/SC3456_sharedsigs/SC3456_Union_all_depleted_Q_subtracted0.1.csv'

######################
#Read in files
rlist1 = pd.read_csv(list1, encoding='utf-8', header=0)
rlist1b = rlist1['Symbol']
wrlist1 = rlist1b.to_csv(output + "Old_enriched_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist1 = pd.read_csv(output + "Old_enriched_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])
rlist2 = pd.read_csv(list2, encoding='utf-8', header=0, index_col='Symbol').drop_duplicates(keep='first')
rlist3 = pd.read_csv(list3, encoding='utf-8', header=0, index_col='Symbol').drop_duplicates(keep='first')
rlist4 = pd.read_csv(list4, encoding='utf-8', header=0)
rlist4b = rlist4['Symbol']
wrlist4 = rlist4b.to_csv(output + "Young_enriched_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist4 = pd.read_csv(output + "Young_enriched_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

rlist5 = pd.read_csv(list5, encoding='utf-8', header=0)
rlist5b = rlist5['Symbol']
wrlist5 = rlist5b.to_csv(output + "Union_enriched_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist5 = pd.read_csv(output + "Union_enriched_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])


rlist6 = pd.read_csv(list6, encoding='utf-8', header=0)
rlist6b = rlist6['Symbol']
wrlist6 = rlist6b.to_csv(output + "Old_depleted_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist6 = pd.read_csv(output + "Old_depleted_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

rlist7 = pd.read_csv(list7, encoding='utf-8', header=0)
rlist7b = rlist7['Symbol']
wrlist7 = rlist7b.to_csv(output + "Young_depleted_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist7 = pd.read_csv(output + "Young_depleted_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

rlist8 = pd.read_csv(list8, encoding='utf-8', header=0)
rlist8b = rlist8['Symbol']
wrlist8 = rlist8b.to_csv(output + "Union_depleted_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist8 = pd.read_csv(output + "Union_depleted_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

rlist9 = pd.read_csv(list9, encoding='utf-8', header=0)
rlist9b = rlist9['Symbol']
wrlist9 = rlist9b.to_csv(output + "SC3456_Old_depleted_Q_subtracted" + ".csv", encoding='utf-8', index=False)
rwrlist9 = pd.read_csv(output + "SC3456_Old_depleted_Q_subtracted" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

rlist10 = pd.read_csv(list10, encoding='utf-8', header=0)
rlist10b = rlist10['Symbol']
wrlist10 = rlist10b.to_csv(output + "SC3456_Young_depleted_Q_subtracted" + ".csv", encoding='utf-8', index=False)
rwrlist10 = pd.read_csv(output + "SC3456_Young_depleted_Q_subtracted" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])


######################
#Make master old enriched lists
union_list1 = pd.merge(rwrlist1, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union1 = union_list1.to_csv(output + "SC345678_Oldspecific_all_enriched_SCORESd4" + ".csv", encoding='utf-8')

union_list2 = pd.merge(rwrlist1, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union2 = union_list2.to_csv(output + "SC345678_Oldspecific_all_enriched_SCORESd14" + ".csv", encoding='utf-8')

union_list3 = rwrlist1.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union3 = union_list3.to_csv(output + "SC345678_Oldspecific_all_enriched_SCORES_alltimepoints" + ".csv", encoding='utf-8')
#
# ###################'
#Make master young enriched lists
union_list4 = pd.merge(rwrlist4, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union4 = union_list4.to_csv(output + "SC345678_Youngspecific_all_enriched_SCORESd4" + ".csv", encoding='utf-8')

union_list5 = pd.merge(rwrlist4, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union5 = union_list5.to_csv(output + "SC345678_Youngspecific_all_enriched_SCORESd14" + ".csv", encoding='utf-8')

union_list6 = rwrlist4.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union6 = union_list6.to_csv(output + "SC345678_Youngspecific_all_enriched_SCORES_alltimepoints" + ".csv", encoding='utf-8')
#
# ###################
#Make Union enriched lists
union_list7 = pd.merge(rwrlist5, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union7 = union_list7.to_csv(output + "SC345678_Union_all_enriched_SCORESd4" + ".csv", encoding='utf-8')

union_list8 = pd.merge(rwrlist5, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union8 = union_list8.to_csv(output + "SC345678_Union_all_enriched_SCORESd14" + ".csv", encoding='utf-8')

union_list9 = rwrlist5.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union9 = union_list9.to_csv(output + "SC345678_Union_all_enriched_SCORES_alltimepoints" + ".csv", encoding='utf-8')

###################
#Make Old depleted lists
union_list1 = pd.merge(rwrlist6, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union1 = union_list1.to_csv(output + "SC345678_Oldspecific_all_depleted_SCORESd4" + ".csv", encoding='utf-8')

union_list2 = pd.merge(rwrlist6, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union2 = union_list2.to_csv(output + "SC345678_Oldspecific_all_depleted_SCORESd14" + ".csv", encoding='utf-8')

union_list3 = rwrlist6.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union3 = union_list3.to_csv(output + "SC345678_Oldspecific_all_depleted_SCORES_alltimepoints" + ".csv", encoding='utf-8')
#
# ###################
#Make young depleted lists
union_list4 = pd.merge(rwrlist7, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union4 = union_list4.to_csv(output + "SC345678_Youngspecific_all_depleted_SCORESd4" + ".csv", encoding='utf-8')

union_list5 = pd.merge(rwrlist7, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union5 = union_list5.to_csv(output + "SC345678_Youngspecific_all_depleted_SCORESd14" + ".csv", encoding='utf-8')

union_list6 = rwrlist7.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union6 = union_list6.to_csv(output + "SC345678_Youngspecific_all_depleted_SCORES_alltimepoints" + ".csv", encoding='utf-8')

###################
#Make union depleted list
union_list7 = pd.merge(rwrlist8, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union7 = union_list7.to_csv(output + "SC345678_Union_all_depleted_SCORESd4" + ".csv", encoding='utf-8')

union_list8 = pd.merge(rwrlist8, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union8 = union_list8.to_csv(output + "SC345678_Union_all_depleted_SCORESd14" + ".csv", encoding='utf-8')

union_list9 = rwrlist8.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union9 = union_list9.to_csv(output + "SC345678_Union_all_depleted_SCORES_alltimepoints" + ".csv", encoding='utf-8')



########################################################################################################
########################################################################################################
#Make same lists for first 2 screens only
output = str('Results/Final_gene_lists/SC3456_final_gene_lists/')

list1 = 'Results/SC3456_sharedsigs/SC3456_Oldspecific_all_enriched_fdr0.1.csv'
list2 = 'Results/PCA/SC345678_allgenescores_d4.csv'
list3 = 'Results/PCA/SC345678_allgenescores_d14.csv'
list4 = 'Results/SC3456_sharedsigs/SC3456_Youngspecific_all_enriched_fdr0.1.csv'
list5 = 'Results/SC3456_sharedsigs/SC3456_Union_enriched_fdr0.1.csv'
list6 = 'Results/SC3456_sharedsigs/SC3456_Oldspecific_all_depleted_fdr0.1.csv'
list7 = 'Results/SC3456_sharedsigs/SC3456_Youngspecific_all_depleted_fdr0.1.csv'
list8 = 'Results/SC3456_sharedsigs/SC3456_Union_depleted_fdr0.1.csv'
######################
#Read in files
rlist1 = pd.read_csv(list1, encoding='utf-8', header=0)
rlist1b = rlist1['Symbol']
wrlist1 = rlist1b.to_csv(output + "Old_enriched_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist1 = pd.read_csv(output + "Old_enriched_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])
rlist2 = pd.read_csv(list2, encoding='utf-8', header=0, index_col='Symbol').drop_duplicates(keep='first')
rlist3 = pd.read_csv(list3, encoding='utf-8', header=0, index_col='Symbol').drop_duplicates(keep='first')
rlist4 = pd.read_csv(list4, encoding='utf-8', header=0)
rlist4b = rlist4['Symbol']
wrlist4 = rlist4b.to_csv(output + "Young_enriched_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist4 = pd.read_csv(output + "Young_enriched_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

rlist5 = pd.read_csv(list5, encoding='utf-8', header=0)
rlist5b = rlist5['Symbol']
wrlist5 = rlist5b.to_csv(output + "Union_enriched_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist5 = pd.read_csv(output + "Union_enriched_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])


rlist6 = pd.read_csv(list6, encoding='utf-8', header=0)
rlist6b = rlist6['Symbol']
wrlist6 = rlist6b.to_csv(output + "Old_depleted_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist6 = pd.read_csv(output + "Old_depleted_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

rlist7 = pd.read_csv(list7, encoding='utf-8', header=0)
rlist7b = rlist4['Symbol']
wrlist7 = rlist7b.to_csv(output + "Young_depleted_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist7 = pd.read_csv(output + "Young_depleted_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

rlist8 = pd.read_csv(list8, encoding='utf-8', header=0)
rlist8b = rlist5['Symbol']
wrlist8 = rlist8b.to_csv(output + "Union_depleted_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist8 = pd.read_csv(output + "Union_depleted_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])


######################
#Make master old enriched lists
union_list1 = pd.merge(rwrlist1, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union1 = union_list1.to_csv(output + "SC3456_Oldspecific_all_enriched_SCORESd4" + ".csv", encoding='utf-8')

union_list2 = pd.merge(rwrlist1, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union2 = union_list2.to_csv(output + "SC3456_Oldspecific_all_enriched_SCORESd14" + ".csv", encoding='utf-8')

union_list3 = rwrlist1.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union3 = union_list3.to_csv(output + "SC3456_Oldspecific_all_enriched_SCORES_alltimepoints" + ".csv", encoding='utf-8')
#
# ###################'
#Make master young enriched lists
union_list4 = pd.merge(rwrlist4, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union4 = union_list4.to_csv(output + "SC3456_Youngspecific_all_enriched_SCORESd4" + ".csv", encoding='utf-8')

union_list5 = pd.merge(rwrlist4, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union5 = union_list5.to_csv(output + "SC3456_Youngspecific_all_enriched_SCORESd14" + ".csv", encoding='utf-8')

union_list6 = rwrlist4.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union6 = union_list6.to_csv(output + "SC3456_Youngspecific_all_enriched_SCORES_alltimepoints" + ".csv", encoding='utf-8')
#
# ###################
#Make Union enriched lists
union_list7 = pd.merge(rwrlist5, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union7 = union_list7.to_csv(output + "SC3456_Union_all_enriched_SCORESd4" + ".csv", encoding='utf-8')

union_list8 = pd.merge(rwrlist5, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union8 = union_list8.to_csv(output + "SC3456_Union_all_enriched_SCORESd14" + ".csv", encoding='utf-8')

union_list9 = rwrlist5.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union9 = union_list9.to_csv(output + "SC3456_Union_all_enriched_SCORES_alltimepoints" + ".csv", encoding='utf-8')

###################
#Make Old depleted lists
union_list1 = pd.merge(rwrlist6, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union1 = union_list1.to_csv(output + "SC3456_Oldspecific_all_depleted_SCORESd4" + ".csv", encoding='utf-8')

union_list2 = pd.merge(rwrlist6, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union2 = union_list2.to_csv(output + "SC3456_Oldspecific_all_depleted_SCORESd14" + ".csv", encoding='utf-8')

union_list3 = rwrlist6.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union3 = union_list3.to_csv(output + "SC3456_Oldspecific_all_depleted_SCORES_alltimepoints" + ".csv", encoding='utf-8')
#
# ###################
#Make young depleted lists
union_list4 = pd.merge(rwrlist7, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union4 = union_list4.to_csv(output + "SC3456_Youngspecific_all_depleted_SCORESd4" + ".csv", encoding='utf-8')

union_list5 = pd.merge(rwrlist7, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union5 = union_list5.to_csv(output + "SC3456_Youngspecific_all_depleted_SCORESd14" + ".csv", encoding='utf-8')

union_list6 = rwrlist7.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union6 = union_list6.to_csv(output + "SC3456_Youngspecific_all_depleted_SCORES_alltimepoints" + ".csv", encoding='utf-8')

###################
#Make union depleted list
union_list7 = pd.merge(rwrlist8, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union7 = union_list7.to_csv(output + "SC3456_Union_all_depleted_SCORESd4" + ".csv", encoding='utf-8')

union_list8 = pd.merge(rwrlist8, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union8 = union_list8.to_csv(output + "SC3456_Union_all_depleted_SCORESd14" + ".csv", encoding='utf-8')

union_list9 = rwrlist8.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union9 = union_list9.to_csv(output + "SC3456_Union_all_depleted_SCORES_alltimepoints" + ".csv", encoding='utf-8')

###################
#Make SC3456 Old depleted with Q subtracted list
union_list10 = pd.merge(rwrlist9, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union10 = union_list10.to_csv(output + "SC3456_Old_depleted_Q_subtracted_SCORESd4" + ".csv", encoding='utf-8')

union_list11 = pd.merge(rwrlist9, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union11 = union_list11.to_csv(output + "SC3456_Old_depleted_Q_subtracted_SCORESd14" + ".csv", encoding='utf-8')

union_list12 = rwrlist9.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union12 = union_list12.to_csv(output + "SC3456_Old_depleted_Q_subtracted_SCORES_alltimepoints" + ".csv", encoding='utf-8')

###################
#Make SC3456 Old depleted with Q subtracted list
union_list10 = pd.merge(rwrlist10, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union10 = union_list10.to_csv(output + "SC3456_Young_depleted_Q_subtracted_SCORESd4" + ".csv", encoding='utf-8')

union_list11 = pd.merge(rwrlist10, rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union11 = union_list11.to_csv(output + "SC3456_Young_depleted_Q_subtracted_SCORESd14" + ".csv", encoding='utf-8')

union_list12 = rwrlist10.merge(rlist2, on='Symbol', how='inner').drop_duplicates(keep='first').merge(rlist3, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union12 = union_list12.to_csv(output + "SC3456_Young_depleted_Q_subtracted_SCORES_alltimepoints" + ".csv", encoding='utf-8')
