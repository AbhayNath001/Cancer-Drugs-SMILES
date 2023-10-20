import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from os import listdir
from os.path import isdir, join
import matplotlib.pyplot as plt
import random
import warnings
warnings.filterwarnings('ignore')

folders_upper = []

with open('prohibitive_restrictions.txt', 'r') as f:
    folders_upper.extend(f.readlines())

with open('aggregating_restrictions.txt', 'r') as f:
    folders_upper.extend(f.readlines())

with open('no_restrictions.txt', 'r') as f:
    folders_upper.extend(f.readlines())
    
folders_upper = [x.strip() for x in folders_upper]

folders_all = []

for fold in folders_upper:
    all_obj = listdir(fold)
    folders_cur = []
    for q in all_obj:
        if isdir(join(fold, q)):
            folders_cur.append(join(fold, q))
    if len(folders_cur) == 0:
        folders_cur.append(fold)
    folders_all.extend(folders_cur)
    
number_of_colors = len(folders_all)

color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

total_data = pd.read_csv('clustering_results/total_data.csv', header=None)
filenames = total_data[0]
total_data[0] = [0] + [np.nan] * (len(total_data) - 1)

folders_colors = dict()
for x in range(len(folders_all)):
    folders_colors[folders_all[x]] = color[x]
    
clusters_colors = dict()

for x in range(len(filenames)):
    wo_parsed = filenames[x].replace('_parsed_data', '')
    for k, v in folders_colors.items():
        if k in wo_parsed:
            clusters_colors[x] = v
            break
            
cur_cluster = len(total_data)

Z_list = []

cols_reversed = list(reversed(total_data.columns))
prev = total_data[cols_reversed[0]].values
cur_cluster = len(total_data)

d = dict()
cols = dict()
dist = 1

for x in range(len(total_data)):
    d[x] = x
    cols[x] = 1

for x in cols_reversed[1:]:
    if dist % 1000 == 0:
        print(str(dist) + ' objects have been processed')
    cur = total_data[x].values
    fst = -1
    snd = -1
    for k in list(zip(prev, cur)):
        if not np.isnan(k[1]):
            fst = int(k[1])
        if not np.isnan(k[0]) and np.isnan(k[1]):
            snd = int(k[0])
            break
    fst_real = d[fst]
    snd_real = d[snd]
    cols_new = cols[fst_real] + cols[snd_real]
    del d[snd]
    d[fst] = cur_cluster
    if cols[fst_real] > cols[snd_real]:
        clusters_colors[cur_cluster] = clusters_colors[fst_real]
    else:
        clusters_colors[cur_cluster] = clusters_colors[snd_real]
    cols[cur_cluster] = cols_new
    del cols[fst_real]
    del cols[snd_real]
    cur_cluster += 1
    
    Z_list.append([fst_real, snd_real, dist + 0.0, cols_new])
    dist += 1
    prev = cur

print('building dendrogram...')
plt.figure(figsize=(150,150))
xx = dendrogram(Z_list, link_color_func=lambda x: clusters_colors[x])
plt.savefig('dendrogram.png', format='png', bbox_inches='tight')
print("Done! Dendrogram saved to file 'dendrogram.png'")