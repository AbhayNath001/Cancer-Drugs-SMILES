import pandas as pd
import numpy as np
import random
import sys


if __name__ == '__main__':

    vals = pd.read_csv('values_new.txt', sep=' ', header=None)
    columns = pd.read_csv('columns_new.txt', header=None)
    vals.columns = columns[0]

    genes = pd.read_csv('genes_new.txt', header=None)
    genes.columns = ['gene']

    vals_total = pd.concat([genes, vals], axis=1)
    rules_list = []

    with open('final_rules.txt', 'r') as f:
        res = f.readlines()
        for i, r in enumerate(res):
            if 'Rule' in r:
                rules_list.append(res[i+2])
                
    while True:
        rule_num = input('Enter rule number to consider, from {} to {}: '.format(1, len(rules_list)))
        try:
            rule_num_int = int(rule_num)
            if rule_num_int < 1 or rule_num_int > len(rules_list):
                print('Rule number out of range, please try once again.')
                continue
            else:
                break
        except:
            print('Rule number is not an integer, please try once again.')

    start_index = None
    for i, x in enumerate(res):
        if rules_list[rule_num_int - 1] in x:
            start_index = i
            break
    
    clusters_list = []
    while start_index < len(res):
        if '----------' in res[start_index]:
            break
        if 'Cluster' in res[start_index]:
            clusters_list.append(res[start_index].split(' ')[1].replace(':', ''))
        start_index += 1
    
    with open('hierarchy_new.txt', 'r') as f:
        res_hier = f.readlines()
        hier_lists = []
        for i, r in enumerate(res_hier):
            hier_lists.append((i, r.strip().split(' ')[1:]))

    files_order = []
    for k in clusters_list:
        for l in hier_lists:
            if k not in l[1]:
                continue
            if l[0] in files_order:
                continue
            files_order.append(l[0])
    
    files_order_dict = dict()
    for i, k in enumerate(files_order):
        files_order_dict[k] = i
    
    needed_rules = rules_list[rule_num_int - 1].strip().replace('IF ', '').replace('THEN:', '').split(' AND ')
    needed_rules_true = []
    genes = []
    signs = []
    values = []
    for k in needed_rules:
        needed_rule = k.split(' ')
        needed_rules_true.append(needed_rule)
        genes.append(needed_rule[0])
        signs.append(needed_rule[1])
        values.append(float(needed_rule[2]))

    good_list = []
    cols = vals_total.columns
    vals = vals_total[vals_total.gene.isin(genes)].values

    for i in range(len(cols)):
        if i == 0:
            continue
        ok = True
        for j in range(0, len(vals)):
            if signs[j] == '<' and int(vals[j][i]) >= values[j]:
                ok = False
            if signs[j] == '>' and int(vals[j][i]) <= values[j]:
                ok = False
            if signs[j] == '=' and int(vals[j][i]) != values[j]:
                ok = False
        if ok:
            good_list.append((files_order_dict[i-1], cols[i]))
            
    with open('rule_{}.txt'.format(rule_num_int), 'w') as f:
        f.write('Rule #{}\n'.format(rule_num_int))
        f.write('Samples that satisfy the rule: \n')
        for k in list(sorted(good_list)):
            f.write(k[1] + '\n')

