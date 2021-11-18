#cd Desktop/Dropbox/TR_Cas9Paper_CodeCheck


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas import DataFrame, Series

output = "Results/SC3456_sharedsigs/"

list1 = 'Results/SC3456_sharedsigs/SC3456_Oldspecific_all_enriched_fdr0.1.csv'
list2 = 'Results/SC3456_sharedsigs/SC3456_castlescores.csv'
list4 = 'Results/SC3456_sharedsigs/SC3456_Youngspecific_all_enriched_fdr0.1.csv'
#list5 = 'Results/SC3456_sharedsigs/depleted/SC3456_Union_depleted_fdr0.1.csv'

list6 = 'Results/SC3456_sharedsigs/SC3456_Oldspecific_all_depleted_fdr0.1.csv'
list7 = 'Results/SC3456_sharedsigs/SC3456_Youngspecific_all_depleted_fdr0.1.csv'
list8 = 'Results/SC3456_sharedsigs/depleted/SC3456_Union_depleted_fdr0.1.csv'

rlist1 = pd.read_csv(list1, encoding='utf-8', header=0)
rlist1b = rlist1['Symbol']
wrlist1 = rlist1b.to_csv(output + "Old_enriched_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist1 = pd.read_csv(output + "Old_enriched_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])
rlist2 = pd.read_csv(list2, encoding='utf-8', header=0, index_col='Symbol').drop_duplicates(keep='first')
#rlist3 = pd.read_csv(list3, encoding='utf-8', header=0, index_col='Symbol').drop_duplicates(keep='first')
rlist4 = pd.read_csv(list4, encoding='utf-8', header=0)
rlist4b = rlist4['Symbol']
wrlist4 = rlist4b.to_csv(output + "Young_enriched_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist4 = pd.read_csv(output + "Young_enriched_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

rlist6 = pd.read_csv(list6, encoding='utf-8', header=0)
rlist6b = rlist6['Symbol']
wrlist6 = rlist6b.to_csv(output + "Old_depleted_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist6 = pd.read_csv(output + "Old_depleted_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

rlist7 = pd.read_csv(list7, encoding='utf-8', header=0)
rlist7b = rlist7['Symbol']
wrlist7 = rlist7b.to_csv(output + "Young_depleted_symbols" + ".csv", encoding='utf-8', index=False)
rwrlist7 = pd.read_csv(output + "Young_depleted_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])
# rlist5 = pd.read_csv(list5, encoding='utf-8', header=0)
# rlist5b = rlist5['Symbol']
# wrlist5 = rlist5b.to_csv(output + "Union_enriched_symbols" + ".csv", encoding='utf-8', index=False)
# rwrlist5 = pd.read_csv(output + "Union_enriched_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])


# rlist1 = pd.read_csv(list1, encoding='utf-8', header=0)
# rlist1b = rlist1['Symbol']
# wrlist1 = rlist1b.to_csv(output + "Old_depleted_symbols" + ".csv", encoding='utf-8', index=False)
# rwrlist1 = pd.read_csv(output + "Old_depleted_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])
# rlist2 = pd.read_csv(list2, encoding='utf-8', header=0, index_col='Symbol').drop_duplicates(keep='first')
# rlist3 = pd.read_csv(list3, encoding='utf-8', header=0, index_col='Symbol').drop_duplicates(keep='first')
# rlist4 = pd.read_csv(list4, encoding='utf-8', header=0)
# rlist4b = rlist4['Symbol']
# wrlist4 = rlist4b.to_csv(output + "Young_depleted_symbols" + ".csv", encoding='utf-8', index=False)
# rwrlist4 = pd.read_csv(output + "Young_depleted_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

# rlist5 = pd.read_csv(list5, encoding='utf-8', header=0)
# rlist5b = rlist5['Symbol']
# wrlist5 = rlist5b.to_csv(output + "Union_depleted_symbols" + ".csv", encoding='utf-8', index=False)
# rwrlist5 = pd.read_csv(output + "Union_depleted_symbols" + ".csv", encoding='utf-8', index_col='Symbol', header=None, names=['Symbol'])

union_list1 = pd.merge(rwrlist1, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union1 = union_list1.to_csv(output + "SC3456_Oldspecific_all_enriched_SCORES" + ".csv", encoding='utf-8')

union_list4 = pd.merge(rwrlist4, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union4 = union_list4.to_csv(output + "SC3456_Youngspecific_all_enriched_SCORES" + ".csv", encoding='utf-8')

#######################
#Depleted
union_list6 = pd.merge(rwrlist6, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union6 = union_list6.to_csv(output + "SC3456_Oldspecific_all_depleted_SCORES" + ".csv", encoding='utf-8')

# ###################
union_list7 = pd.merge(rwrlist7, rlist2, on='Symbol', how='inner').drop_duplicates(keep='first')
writefinal_union7 = union_list7.to_csv(output + "SC3456_Youngspecific_all_depleted_SCORES" + ".csv", encoding='utf-8')
