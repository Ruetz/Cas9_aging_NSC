##concat screen files and run PCA
##Working Directory: cd Desktop/Dropbox/TR_Cas9Paper_CodeCheck

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np

output = "Results/PCA/"

data1 = 'Results/SC3_RA_FDR0.1.csv'
data2 = 'Results/SC3_Kp_FDR0.1.csv'
data3 = 'Results/SC4_RA_FDR0.1.csv'
data4 = 'Results/SC4_Kp_FDR0.1.csv'
data5 = 'Results/SC5_RA_FDR0.1.csv'
data6 = 'Results/SC5_Kp_FDR0.1.csv'
data7 = 'Results/SC6_RA_FDR0.1.csv'
data8 = 'Results/SC6_Kp_FDR0.1.csv'
data9 = 'Results/SC7_RA_FDR0.1.csv'
data10 = 'Results/SC7_Kp_FDR0.1.csv'
data11 = 'Results/SC8_RA_FDR0.1.csv'
data12 = 'Results/SC8_Kp_FDR0.1.csv'


# data1 = 'SC3_RAvLib.csv'
# data2 = 'SC3_Kpvlib.csv'
# data3 = 'SC4_RAvLib.csv'
# data4 = 'SC4_KpvLib.csv'
# data5 = 'SC5_RAvlib.csv'
# data6 = 'SC5_KpvLib.csv'
# data7 = 'SC6_RAvLib.csv'
# data8 = 'SC6_KpvLib.csv'
# data9 = 'SC7_RAvLib.csv'
# data10 = 'SC7_KpvLib.csv'
# data11 = 'SC8_RAvLib.csv'
# data12 = 'SC8_KpvLib.csv'

stdata1 = str('Old1')
stdata2 = str('Old1')
stdata3 = str('Young1')
stdata4 = str('Young1')
stdata5 = str('Young2')
stdata6 = str('Young2')
stdata7 = str('Old2')
stdata8 = str('Old2')
stdata9 = str('Young3')
stdata10 = str('Young3')
stdata11 = str('Old3')
stdata12 = str('Old3')



rdata1 = pd.read_csv(data1, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate',
'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata2 = pd.read_csv(data2, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate'
, 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata3 = pd.read_csv(data3, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate'
, 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata4 = pd.read_csv(data4, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate'
, 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata5 = pd.read_csv(data5, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate'
, 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata6 = pd.read_csv(data6, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate'
, 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata7 = pd.read_csv(data7, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate'
, 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata8 = pd.read_csv(data8, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate'
, 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata9 = pd.read_csv(data9, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate'
, 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata10 = pd.read_csv(data10, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate'
, 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata11 = pd.read_csv(data11, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate'
, 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])
rdata12 = pd.read_csv(data12, encoding='utf-8', header=None, skiprows=1, index_col='#GeneID', names=['Symbol', '#GeneID', 'GeneInfo', 'Localization', 'Process', 'Function', 'Element #', 'casTLE Effect', 'casTLE Score', 'casTLE p-value', 'FDR_0.1', 'Minimum Effect Estimate', 'Maximum Effect Estimate'
, 'Individual Elements', 'Unnamed: 13', 'Unnamed: 14', 'Unnamed: 15', 'Unnamed: 16', 'Unnamed: 17', 'Unnamed: 18', 'Unnamed: 19', 'Unnamed: 20', 'Unnamed: 21'])



##Plotting d14 PCA's
concat = pd.concat([rdata1['Symbol'], (rdata1['casTLE Score'].rename(stdata1)), (rdata3['casTLE Score'].rename(stdata3)), (rdata5['casTLE Score'].rename(stdata5)), (rdata7['casTLE Score'].rename(stdata7)), (rdata9['casTLE Score'].rename(stdata9)),
(rdata11['casTLE Score'].rename(stdata11))], axis=1, ignore_index=False, join='inner')
final = concat.to_csv(output + 'castlescores.csv')

concat = pd.concat([rdata1['Symbol'], (rdata1['casTLE Effect'].rename(stdata1 + ' casTLE Score')), (rdata3['casTLE Effect'].rename(stdata3)), (rdata5['casTLE Effect'].rename(stdata5)), (rdata7['casTLE Effect'].rename(stdata7)), (rdata9['casTLE Effect'].rename(stdata9)),
(rdata11['casTLE Effect'].rename(stdata11))], axis=1, ignore_index=False, join='inner')
final = concat.to_csv(output + 'castleEffect.csv')

#make a master list of genes with score, effect, FDR p-value
concat = pd.concat([rdata1['Symbol'], (rdata3['casTLE Score'].rename(stdata3 + ' casTLE Score')), (rdata3['casTLE Effect'].rename(stdata3 + ' casTLE Effect')), (rdata3['FDR_0.1'].rename(stdata3 + ' FDR p-value')), (rdata5['casTLE Score'].rename(stdata5 + ' casTLE Score')), (rdata5['casTLE Effect'].rename(stdata5 + ' casTLE Effect')), (rdata5['FDR_0.1'].rename(stdata5 + ' FDR p-value')), (rdata9['casTLE Score'].rename(stdata9 + ' casTLE Score')), (rdata9['casTLE Effect'].rename(stdata9 + ' casTLE Effect')), (rdata9['FDR_0.1'].rename(stdata9 + ' FDR p-value')), (rdata1['casTLE Score'].rename(stdata1 + ' casTLE Score')), (rdata1['casTLE Effect'].rename(stdata1 + ' casTLE Effect')), (rdata1['FDR_0.1'].rename(stdata1 + ' FDR p-value')), (rdata7['casTLE Score'].rename(stdata7 + ' casTLE Score')), (rdata7['casTLE Effect'].rename(stdata7 + ' casTLE Effect')), (rdata7['FDR_0.1'].rename(stdata7 + ' FDR p-value')), (rdata11['casTLE Score'].rename(stdata11 + ' casTLE Score')), (rdata11['casTLE Effect'].rename(stdata11 + ' casTLE Effect')), (rdata11['FDR_0.1'].rename(stdata11 + ' FDR p-value'))], axis=1, ignore_index=False, join='inner')
final = concat.to_csv(output + 'SC345678_allgenescores_d14.csv')


Adata = pd.read_csv(output + 'castlescores.csv', header=None, skiprows=1, names=['#GeneID', 'Symbol', stdata1, stdata3, stdata5, stdata7, stdata9, stdata11])
tAdata = Adata.transpose()
wtAdata = tAdata.to_csv(output + 'castlescores_transpose.csv')
rdata = pd.read_csv(output + 'castlescores_transpose.csv', skiprows=1, header=1)
print(rdata.head())

targets = [stdata1, stdata3, stdata5, stdata7, stdata9, stdata11]
features = tAdata.loc['Symbol', :]
x = rdata.loc[:, features].values
y = rdata.loc[:, ['Symbol']].values

y = StandardScaler().fit_transform(x)

pca = PCA(n_components=4)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2', 'principal component 3', 'principal component 4'])

print(principalDf)

finalDf = pd.concat([principalDf, rdata[['Symbol']]], axis = 1)

print(finalDf)
explained_variance = str(pca.explained_variance_ratio_)
print(explained_variance)

fig = plt.figure(figsize = (14,14))
ax = fig.add_subplot(1,1,1)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.set_xlabel('Principal Component 1', fontsize = 30)
ax.set_ylabel('Principal Component 2', fontsize = 30)
ax.set_xlabel('Principal Component 1' + ' (' + explained_variance[3:5] + '%)', fontsize = 30)
ax.set_ylabel('Principal Component 2' + ' (' + explained_variance[14:16] + '%)', fontsize = 30)
ax.set_title('PCA Castle Scores Day 14', fontsize = 30)

# Label to color dict (automatic)
label_color_dict = {target:idx for idx,target in enumerate(np.unique(targets))}
colors = ['xkcd:scarlet', 'xkcd:cyan', 'xkcd:cyan', 'xkcd:scarlet', 'xkcd:cyan', 'xkcd:scarlet']

for target, color in zip(targets, colors):
    indicesToKeep = finalDf['Symbol'] == target
    ax.scatter(finalDf['principal component 1'], finalDf['principal component 2'], c=colors, s=100, linewidths=4)
    plt.text(finalDf.loc[indicesToKeep, 'principal component 1'], finalDf.loc[indicesToKeep, 'principal component 2'], target, fontsize=30)

##ax.grid()
plt.savefig(output + 'SC345678_PCA_d14_age_PC12_score', dpi=300)

# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot(1,1,1)
# ax.tick_params(axis='both', which='major', labelsize=20)
# ax.set_xlabel('Principal Component 3', fontsize = 30)
# ax.set_ylabel('Principal Component 4', fontsize = 30)
# ax.set_xlabel('Principal Component 3' + ' (' + explained_variance[25:27] + '%)', fontsize = 30)
# ax.set_ylabel('Principal Component 4' + ' (' + explained_variance[36:38] + '%)', fontsize = 30)
# ax.set_title('PCA Castle Scores Day 14', fontsize = 30)

# # Label to color dict (automatic)
# label_color_dict = {target:idx for idx,target in enumerate(np.unique(targets))}
# colors = ['xkcd:scarlet', 'xkcd:cyan', 'xkcd:cyan', 'xkcd:scarlet', 'xkcd:cyan', 'xkcd:scarlet']
#
# for target, color in zip(targets, colors):
#     indicesToKeep = finalDf['Symbol'] == target
#     ax.scatter(finalDf['principal component 3'], finalDf['principal component 4'], c=colors, s=100, linewidths=4)
#     #plt.text(finalDf.loc[indicesToKeep, 'principal component 3'], finalDf.loc[indicesToKeep, 'principal component 4'], target, fontsize=30)
#
# #ax.grid()
# plt.savefig(output + 'SC345678_PCA_d14_age_PC34_score', dpi=300)

#############################
##plotting d14 effect PCA
# Adata = pd.read_csv(output + 'castleEffect.csv', header=None, skiprows=1, names=['#GeneID', 'Symbol', stdata1, stdata3, stdata5, stdata7, stdata9, stdata11])
# tAdata = Adata.transpose()
# wtAdata = tAdata.to_csv(output + 'castlescores_transpose.csv')
# rdata = pd.read_csv(output + 'castlescores_transpose.csv', skiprows=1, header=1)
# print(rdata.head())
#
# targets = [stdata1, stdata3, stdata5, stdata7, stdata9, stdata11]
# features = tAdata.loc['Symbol', :]
# x = rdata.loc[:, features].values
# y = rdata.loc[:, ['Symbol']].values
#
# y = StandardScaler().fit_transform(x)
#
# pca = PCA(n_components=4)
# principalComponents = pca.fit_transform(x)
# principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2', 'principal component 3', 'principal component 4'])
#
# print(principalDf)
#
# finalDf = pd.concat([principalDf, rdata[['Symbol']]], axis = 1)
#
# print(finalDf)
# explained_variance = str(pca.explained_variance_ratio_)
# print(explained_variance)
#
# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot(1,1,1)
# ax.tick_params(axis='both', which='major', labelsize=20)
# # ax.set_xlabel('Principal Component 1', fontsize = 30)
# # ax.set_ylabel('Principal Component 2', fontsize = 30)
# # ax.set_xlabel('Principal Component 1' + ' (' + explained_variance[3:5] + '%)', fontsize = 30)
# # ax.set_ylabel('Principal Component 2' + ' (' + explained_variance[14:16] + '%)', fontsize = 30)
# # ax.set_title('PCA Castle Effect Day 14', fontsize = 30)
#
# # Label to color dict (automatic)
# label_color_dict = {target:idx for idx,target in enumerate(np.unique(targets))}
#
# colors = ['xkcd:scarlet', 'xkcd:cyan', 'xkcd:cyan', 'xkcd:scarlet', 'xkcd:cyan', 'xkcd:scarlet']
#
# for target, color in zip(targets, colors):
#     indicesToKeep = finalDf['Symbol'] == target
#     ax.scatter(finalDf['principal component 1'], finalDf['principal component 2'], c=colors, s=100, linewidths=4)
#     #plt.text(finalDf.loc[indicesToKeep, 'principal component 1'], finalDf.loc[indicesToKeep, 'principal component 2'], target, fontsize=30)
#
#
#
# #ax.grid()
# plt.savefig(output + 'SC345678_PCA_d14_age_PC12_effect', dpi=300)
#
# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot(1,1,1)
# ax.tick_params(axis='both', which='major', labelsize=20)
# # ax.set_xlabel('Principal Component 3', fontsize = 30)
# # ax.set_ylabel('Principal Component 4', fontsize = 30)
# # ax.set_xlabel('Principal Component 3' + ' (' + explained_variance[25:27] + '%)', fontsize = 30)
# # ax.set_ylabel('Principal Component 4' + ' (' + explained_variance[36:38] + '%)', fontsize = 30)
# # ax.set_title('PCA Castle Effect Day 14', fontsize = 30)
#
#
# for target, color in zip(targets, colors):
#     indicesToKeep = finalDf['Symbol'] == target
#     ax.scatter(finalDf['principal component 3'], finalDf['principal component 4'], c=colors, s=100, linewidths=4)
#     #plt.text(finalDf.loc[indicesToKeep, 'principal component 3'], finalDf.loc[indicesToKeep, 'principal component 4'], target, fontsize=30)
#
# #ax.grid()
# plt.savefig(output + 'SC345678_PCA_d14_age_PC34_effect', dpi=300)

#####################
#####################
##NOW d4 PCA/list production
concat = pd.concat([rdata2['Symbol'], (rdata2['casTLE Score'].rename(stdata2)), (rdata4['casTLE Score'].rename(stdata4)), (rdata6['casTLE Score'].rename(stdata6)), (rdata8['casTLE Score'].rename(stdata8)), (rdata10['casTLE Score'].rename(stdata10)),
(rdata12['casTLE Score'].rename(stdata12))], axis=1, ignore_index=False, join='inner')
final = concat.to_csv(output + 'castlescores.csv')

concat = pd.concat([rdata2['Symbol'], (rdata2['casTLE Effect'].rename(stdata2 + ' casTLE Score')), (rdata4['casTLE Effect'].rename(stdata4)), (rdata6['casTLE Effect'].rename(stdata6)), (rdata8['casTLE Effect'].rename(stdata8)), (rdata10['casTLE Effect'].rename(stdata10)),
(rdata12['casTLE Effect'].rename(stdata12))], axis=1, ignore_index=False, join='inner')
final = concat.to_csv(output + 'castleEffect.csv')

concat = pd.concat([rdata2['Symbol'], (rdata4['casTLE Score'].rename(stdata4 + ' casTLE Score')), (rdata4['casTLE Effect'].rename(stdata4 + ' casTLE Effect')), (rdata4['FDR_0.1'].rename(stdata4 + ' FDR p-value')), (rdata6['casTLE Score'].rename(stdata6 + ' casTLE Score')), (rdata6['casTLE Effect'].rename(stdata6 + ' casTLE Effect')), (rdata6['FDR_0.1'].rename(stdata6 + ' FDR p-value')), (rdata10['casTLE Score'].rename(stdata10 + ' casTLE Score')), (rdata10['casTLE Effect'].rename(stdata10 + ' casTLE Effect')), (rdata10['FDR_0.1'].rename(stdata10 + ' FDR p-value')), (rdata2['casTLE Score'].rename(stdata2 + ' casTLE Score')), (rdata2['casTLE Effect'].rename(stdata2 + ' casTLE Effect')), (rdata2['FDR_0.1'].rename(stdata2 + ' FDR p-value')), (rdata8['casTLE Score'].rename(stdata8 + ' casTLE Score')), (rdata8['casTLE Effect'].rename(stdata8 + ' casTLE Effect')), (rdata8['FDR_0.1'].rename(stdata8 + ' FDR p-value')), (rdata12['casTLE Score'].rename(stdata12 + ' casTLE Score')), (rdata12['casTLE Effect'].rename(stdata12 + ' casTLE Effect')), (rdata12['FDR_0.1'].rename(stdata12 + ' FDR p-value'))], axis=1, ignore_index=False, join='inner')
final = concat.to_csv(output + 'SC345678_allgenescores_d4.csv')

##scores
Adata = pd.read_csv(output + 'castlescores.csv', header=None, skiprows=1, names=['#GeneID', 'Symbol', stdata2, stdata4, stdata6, stdata8, stdata10, stdata12])
tAdata = Adata.transpose()
wtAdata = tAdata.to_csv(output + 'castlescores_transpose.csv')
rdata = pd.read_csv(output + 'castlescores_transpose.csv', skiprows=1, header=1)
print(rdata.head())

targets = [stdata2, stdata4, stdata6, stdata8, stdata10, stdata12]
features = tAdata.loc['Symbol', :]
x = rdata.loc[:, features].values
y = rdata.loc[:, ['Symbol']].values

y = StandardScaler().fit_transform(x)

pca = PCA(n_components=4)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2', 'principal component 3', 'principal component 4'])

print(principalDf)

finalDf = pd.concat([principalDf, rdata[['Symbol']]], axis = 1)

print(finalDf)
explained_variance = str(pca.explained_variance_ratio_)
print(explained_variance)

fig = plt.figure(figsize = (14,14))
ax = fig.add_subplot(1,1,1)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.set_xlabel('Principal Component 1', fontsize = 30)
ax.set_ylabel('Principal Component 2', fontsize = 30)
ax.set_xlabel('Principal Component 1' + ' (' + explained_variance[3:5] + '%)', fontsize = 30)
ax.set_ylabel('Principal Component 2' + ' (' + explained_variance[14:16] + '%)', fontsize = 30)
ax.set_title('PCA Castle Scores Day 4', fontsize = 30)

# Label to color dict (automatic)
label_color_dict = {target:idx for idx,target in enumerate(np.unique(targets))}
colors = ['xkcd:scarlet', 'xkcd:cyan', 'xkcd:cyan', 'xkcd:scarlet', 'xkcd:cyan', 'xkcd:scarlet']
for target, color in zip(targets, colors):
    indicesToKeep = finalDf['Symbol'] == target
    ax.scatter(finalDf['principal component 1'], finalDf['principal component 2'], c=colors, s=100, linewidths=4)
    plt.text(finalDf.loc[indicesToKeep, 'principal component 1'], finalDf.loc[indicesToKeep, 'principal component 2'], target, fontsize=30)

#ax.grid()
plt.savefig(output + 'SC345678_PCA_d4_age_PC12_score', dpi=300)

fig = plt.figure(figsize = (14,14))
ax = fig.add_subplot(1,1,1)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.set_xlabel('Principal Component 3', fontsize = 30)
ax.set_ylabel('Principal Component 4', fontsize = 30)
ax.set_xlabel('Principal Component 3' + ' (' + explained_variance[25:27] + '%)', fontsize = 30)
ax.set_ylabel('Principal Component 4' + ' (' + explained_variance[36:38] + '%)', fontsize = 30)
ax.set_title('PCA Castle Effect Day 4', fontsize = 30)

for target, color in zip(targets, colors):
    indicesToKeep = finalDf['Symbol'] == target
    ax.scatter(finalDf['principal component 3'], finalDf['principal component 4'], c=colors, s=100, linewidths=4)
    plt.text(finalDf.loc[indicesToKeep, 'principal component 3'], finalDf.loc[indicesToKeep, 'principal component 4'], target, fontsize=30)

#ax.grid()
plt.savefig(output + 'SC345678_PCA_d4_age_PC34_score', dpi=300)

#############################
# Adata = pd.read_csv(output + 'castleEffect.csv', header=None, skiprows=1, names=['#GeneID', 'Symbol', stdata2, stdata4, stdata6, stdata8, stdata10, stdata12])
# tAdata = Adata.transpose()
# wtAdata = tAdata.to_csv(output + 'castlescores_transpose.csv')
# rdata = pd.read_csv(output + 'castlescores_transpose.csv', skiprows=1, header=1)
# print(rdata.head())
#
# targets = [stdata2, stdata4, stdata6, stdata8, stdata10, stdata12]
# features = tAdata.loc['Symbol', :]
# x = rdata.loc[:, features].values
# y = rdata.loc[:, ['Symbol']].values
#
# y = StandardScaler().fit_transform(x)
#
# pca = PCA(n_components=4)
# principalComponents = pca.fit_transform(x)
# principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2', 'principal component 3', 'principal component 4'])
#
# print(principalDf)
#
# finalDf = pd.concat([principalDf, rdata[['Symbol']]], axis = 1)
#
# print(finalDf)
# explained_variance = str(pca.explained_variance_ratio_)
# print(explained_variance)
#
# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot(1,1,1)
# ax.tick_params(axis='both', which='major', labelsize=20)
# # ax.set_xlabel('Principal Component 1', fontsize = 30)
# # ax.set_ylabel('Principal Component 2', fontsize = 30)
# # ax.set_xlabel('Principal Component 1' + ' (' + explained_variance[3:5] + '%)', fontsize = 30)
# # ax.set_ylabel('Principal Component 2' + ' (' + explained_variance[14:16] + '%)', fontsize = 30)
# # ax.set_title('PCA Castle Effect Day 4', fontsize = 30)
#
# # Label to color dict (automatic)
# label_color_dict = {target:idx for idx,target in enumerate(np.unique(targets))}
# colors = ['xkcd:scarlet', 'xkcd:cyan', 'xkcd:cyan', 'xkcd:scarlet', 'xkcd:cyan', 'xkcd:scarlet']
# for target, color in zip(targets, colors):
#     indicesToKeep = finalDf['Symbol'] == target
#     ax.scatter(finalDf['principal component 1'], finalDf['principal component 2'], c=colors, s=100, linewidths=4)
#     #plt.text(finalDf.loc[indicesToKeep, 'principal component 1'], finalDf.loc[indicesToKeep, 'principal component 2'], target, fontsize=30)
#
# #ax.grid()
# plt.savefig(output + 'SC345678_PCA_d4_age_PC12_effect', dpi=300)

exit()
