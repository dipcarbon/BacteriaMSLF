import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text



my_width = 0.5
my_fontsize = 'large'
blue = '#003D79'
red = '#930000'

gene_name = pd.read_csv('plot_data/weight.csv', ).iloc[:10]
gene_name.columns = ['name','mean','std']
gene_name['name'] = gene_name['name'].apply(lambda x: x[:-1]+x[-1].upper())
gene_name['std'] = gene_name['std'].apply(lambda x: float(x))
gene_name['mean'] = gene_name['mean'].apply(lambda x: float(x))

fig, ax = plt.subplots()
ax.set_axisbelow(True)
ax.yaxis.grid(True,linestyle='--',linewidth=1)
#ax.grid(linestyle = "--",zorder = 1) 
rects1 = ax.bar(range(10), gene_name['mean'] , my_width, color=blue, yerr=gene_name['std'])
plt.xticks(range(10), gene_name['name'], rotation='horizontal', fontsize = my_fontsize)
#ax.set_xticklabels(gene_name['name'])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_ylabel('Relative Gene Weights', fontsize = my_fontsize)
plt.savefig('1.tif', bbox_inches='tight')


match = pd.read_csv('plot_data/matched.csv', )
ax.set_axisbelow(True)
match.columns = ['name','log','mat']
match['name'] = match['name'].apply(lambda x: x[0].upper()+x[1:])
fig, ax = plt.subplots()
rects1 = ax.scatter(match['mat'],match['log'], s=10, color=blue, )
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.plot([-1,50],[10,10], linestyle='--',linewidth=1, color = 'grey')
ax.set_xlabel('Matched Peaks(Avg.)', fontsize = my_fontsize)
ax.set_ylabel('Log2 of Reviewed Protein Numbers', fontsize = my_fontsize)
plt.xlim(-1, 50)
plt.savefig('2.tif', bbox_inches='tight')


cross = pd.read_csv('plot_data/cross.csv', )[['Accuracy', 'Cross-validation', 'Double Cross-validation',]]
fig, ax = plt.subplots()
bar2= [i+my_width/2 for i in range(11)]
rects1 = ax.bar(range(11), cross['Cross-validation'], my_width/2, color=blue,)
rects2 = ax.bar(bar2, cross['Double Cross-validation'], my_width/2, color=red,)
plt.xticks(range(11), cross['Accuracy'], rotation='horizontal', fontsize = my_fontsize)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Frequency', fontsize = my_fontsize)
ax.set_ylabel('Accuracy', fontsize = my_fontsize)
plt.xlim(-1, 10)
ax.yaxis.grid(True,linestyle='--',linewidth=1)
ax.set_axisbelow(True)
ax.legend((rects1[0], rects2[0]), ('Cross-validation', 'Double Cross-validation'))
plt.savefig('3.tif', bbox_inches='tight')

all_file = [i for i in list(cat['Genera'])]
cat = pd.read_csv('plot_data/cat.csv', )
ax= cat.plot(kind='bar', stacked=True, color = [blue,red,'grey'], fontsize = my_fontsize,)
ax.set_axisbelow(True)
ax.yaxis.grid(True,linestyle='--',linewidth=1)
ax.set_xlabel('Genera', fontsize = my_fontsize)
ax.set_ylabel('Number of spectra', fontsize = my_fontsize)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xticks(range(11), all_file, rotation=45, fontsize = my_fontsize)
plt.savefig('4.tif', bbox_inches='tight')


all_file = [i.lower() for i in list(cat['Genera'])]
def getName(my_name):
    uc = pd.read_csv('uniprot/'+my_name+'.csv', )
    rpmh = uc[uc.gn=='rpmh']
    return rpmh
    return list(rpmh['mw'])

all_file = [i.lower() for i in list(cat['Genera'])]
del all_file[-1]
c = pd.concat([getName(i) for i in all_file])
c.columns = ['Unnamed: 0', 'species', 'gn', 'mw', 'genera']
c['genera'] = c['genera'].apply(lambda x: x[0].upper()+x[1:])
ax = c[['mw','genera']].boxplot(by='genera',fontsize = my_fontsize,rot = 45)
#plt.xticks(range(9), all_file, rotation=45, fontsize = my_fontsize)
ax.set_axisbelow(True)
ax.yaxis.grid(False)
ax.xaxis.grid(True,linestyle='--',linewidth=1)
ax.set_xlabel('Bacteria', fontsize = my_fontsize)
ax.set_ylabel('Molecular Weight', fontsize = my_fontsize)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_title('rpmH')

ocr = pd.read_csv('plot_data/ocr.csv', )
ocr.index= ocr.TPR
ocr = ocr[['rpmC', 'rpmJ', 'rpmH', 'Top 3', 'Top 10']]
ax = ocr.plot(color = [blue, red, '#64A600', '#737300', '#007979'])

ax.set_xlabel('False Positive Rate (FPR)', fontsize = my_fontsize)
ax.set_ylabel('True Positive Rate (TPR)', fontsize = my_fontsize)
ax.plot([0,0],[1,1], linestyle='--',linewidth=2, color = 'grey')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('5.tif', bbox_inches='tight')

