import pandas as pd
import numpy as np
import os
import sys

def zapoln(x):
    cur = None
    res = []
    for k in x:
        if not np.isnan(k):
            res.append(k)
            cur = k
        else:
            res.append(cur)
    return res


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("Path to genes or None required!")
        exit()

    genes = None

    if sys.argv[1] != 'None':
        with open(sys.argv[1], 'r') as f:
            genes = [x.strip() for x in f.readlines()]

    xx = pd.read_csv('total_data.csv', header=None)

    for col in xx.columns[1:]:
        xx[col] = zapoln(xx[col])
        
    d = dict()
    d2 = set()

    for k in xx.values:
        hier = set()
        for q in k[1:]:
            if np.isnan(q):
                break
            hier.add(int(q))
        d[k[0].split('\\')[-1]] = ' '.join([str(len(hier))] + [str(x) for x in sorted(list(hier))])
        d2.add(k[0])
        
    all_genes = set()

    for k in d2:
        aa = pd.read_csv(k)
        aa.columns = ['gene', k.split('\\')[-1]]
        all_genes.update(aa.gene)

    table = dict()
    for x in all_genes:
        table[x] = [None] * len(d2)
        
    merged = pd.DataFrame(sorted(list(all_genes)), columns = ['gene'])

    cols = []

    for i, k in enumerate(d2):
        print(i)
        aa = pd.read_csv(k)
        for t in aa.values:
            if t[0] not in table.keys():
                continue
            table[t[0]][i] = t[1]

        cols.append(k.split('\\')[-1])
        aa.columns = ['gene', k.split('\\')[-1]]

    vals = []
    for k, v in table.items():
        lst = [k]
        lst.extend(v)
        vals.append(lst)

    merged = pd.DataFrame(vals, columns=['gene']+list(cols))
    merged.fillna(0, inplace=True)

    genes_needed = []

    for x in merged.values:
        if sum(x[1:]) != 0:
            genes_needed.append(x[0])

    merged_needed = pd.DataFrame(genes_needed, columns=['gene']).merge(merged)
    if genes != None:
        merged_needed = merged_needed[merged_needed['gene'].isin(genes)]

    with open('hierarchy_new.txt', 'w') as f:
        f.write('\n'.join(pd.DataFrame(merged_needed.columns[1:]).merge(pd.DataFrame(d.items()))[1].values))

    with open('columns_new.txt', 'w') as f:
        f.write('\n'.join(merged_needed.columns[1:]))
        
    with open('genes_new.txt', 'w') as f:
        f.write('\n'.join(merged_needed['gene']))
        
    dropped = merged_needed.drop('gene', axis=1)
    dropped = dropped.astype(int)
    dropped.to_csv('values_new.txt', header=None, index=False, sep=' ')

    os.system("g++ count_variances.cpp")
    os.system("./a.out")
