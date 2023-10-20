import pandas as pd
import numpy as np
import random
import sys

# Prompt user for input values
MAX_TERMS = int(input("Enter the maximum number of terms in signature: "))
ANT_COLONY_SIZE = int(input("Enter Ant Colony Size (the number of ants): "))
CONVERGENCE = int(input("Enter the number of iterations yielding the same best signature for identifying signature convergence: "))
MIN_SAMPLES = int(input("Enter the minimum number of samples covered by a signature: "))
MAX_UNCOVERED = int(input("Enter the maximum number of uncovered samples to terminate the process: "))
MAX_WEIGHT = 0.99
MIN_WEIGHT = 0.01

covered_objects = []
cannot_clusters = []
threshold = None


def choose_rules(max_terms):
    
    global variances_table
    
    curres = np.random.uniform(0, 1, max_terms)
    arr_all = variances_table.prob
    arr_cumulative = [arr_all[0]]

    for x in arr_all[1:]:
        arr_cumulative.append(arr_cumulative[-1] + x)
        
    rules_chosen = []

    for x in curres:
        l = 0
        r = len(arr_all) - 1
        while r - l > 1:
            m = (l + r) // 2
            if arr_cumulative[m] < x:
                l = m
            else:
                r = m
        rules_chosen.append(variances_table[variances_table.index == l].values[0][:3])
    
    return rules_chosen


def get_objects(rules):
    global vals_total
    global covered_objects
    
    s = set()

    for i, r in enumerate(rules):
        pp = vals_total[vals_total.gene == r[0]].drop('gene', axis=1).T
        pp.columns = ['value']
        if r[2] == 0:
            total = pp[pp.value == r[1]].reset_index()
        elif r[2] == 1:
            total = pp[pp.value < r[1]].reset_index()
        else:
            total = pp[pp.value > r[1]].reset_index()
        if i == 0:
            s.update(total['index'])
        else:
            s = s.intersection(total['index'])
            
    s = s - set(covered_objects)
    return list(s)


def get_output(objects):
    global hier_total
    
    h_conc_needed = hier_total[hier_total[0].isin(objects)]
    
    total_list = []
    for x in list(h_conc_needed.hierarchy):
        total_list.extend(x)
        
    output = dict((i, total_list.count(i)) for i in total_list)

    res = []
    for k, v in output.items():
        res.append((v * 1.0 / len(objects), k))

    return list(reversed(sorted(res)))


def get_variance(hier):
    
    global columns
    global MAX_WEIGHT
    global MIN_WEIGHT
    
    col = len(hier[0])
    col_all = len(columns[0])
    
    vals_total = []

    for x in hier.hierarchy:
        vals_total.extend(x)

    raspr_avg = dict((int(i), vals_total.count(i) * 1.0 / col) for i in vals_total)

    var_total = 0

    for x in hier.hierarchy:
        raspr_cur = {int(k) for k in x}
        diff = {k: v for k, v in raspr_avg.items()}
        for k in raspr_cur:
            if k in raspr_avg.keys():
                diff[k] -= 1
            else:
                diff[k] = -1

        var_sq = 0
        for k, v in diff.items():
            cur_weight = MIN_WEIGHT + (MAX_WEIGHT - MIN_WEIGHT) * (col_all - k - 1) * 1.0 / (col_all - 1)
            var_sq += cur_weight * v * v
        var_total += var_sq

    return var_total / col


def get_var_gain(objects):
    
    if len(objects) == 0:
        return None

    global hier_total
    
    var_all = get_variance(hier_total)
    
    covered = hier_total[hier_total[0].isin(objects)]
    uncovered = hier_total[~hier_total[0].isin(objects)]
    
    var_covered = get_variance(covered)
    
    if len(uncovered) == 0:
        var_uncovered = 0
    else:
        var_uncovered = get_variance(uncovered)
    
    return var_all - len(covered) * var_covered / len(hier_total) - len(uncovered) * var_uncovered / len(hier_total)


def choose_rules_smart(max_terms):
    
    global MAX_TERMS
    global MIN_SAMPLES
    global cannot_clusters
    global threshold
    
    for iter in range(100):
        rules = choose_rules(MAX_TERMS)
        objects = get_objects(rules)

        cannot_ok = False
        if threshold is None:
            cannot_ok = True

        while True:
            if len(rules) == 0:
                break

            if len(objects) >= MIN_SAMPLES and cannot_ok:
                break

            rules = rules[:-1]
            objects = get_objects(rules)
            output = get_output(objects)

            cur_ok = True
            if threshold is None:
                continue
            for k, v in output:
                if v in cannot_clusters and k > threshold:
                    cur_ok = False
                    break
            cannot_ok = cur_ok

        if len(rules) > 0:
            return rules
        
    return None


def prune_rules(rules):
    
    objects = get_objects(rules)
    max_var_gain = get_var_gain(objects)
    best_rules = rules

    while True:

        if len(best_rules) == 1:
            break
        is_change = False
        for i in range(len(best_rules)):
            next_rules = []
            for j in range(len(best_rules)):
                if i != j:
                    next_rules.append(best_rules[j])

            next_objects = get_objects(next_rules)
            next_var_gain = get_var_gain(next_objects)

            if next_var_gain >= max_var_gain:
                max_var_gain = next_var_gain
                best_rules = next_rules
                is_change = True

        if not is_change:
            break
    
    return best_rules


def count_rules_quality(rules):
    global hier_total
    global columns
    
    objects = get_objects(rules)
    output = get_output(objects)
    
    vals_total = []

    for x in hier_total.hierarchy:
        vals_total.extend(x)
    
    raspr_total = dict((int(i), vals_total.count(i)) for i in vals_total)
    output_obj = dict()

    for x in output:
        output_obj[int(x[1])] = x[0] * len(objects)
        
    best_Q = 0
    best_class = None

    for i, x in output_obj.items():
        TP = x
        FP = len(objects) - x
        FN = raspr_total[i] - x
        TN = len(columns) - TP - FP - FN
        cur_Q = (TP * 1.0 / (TP + FN)) * (TN * 1.0 / (FP + TN))
        if cur_Q > best_Q:
            best_Q = cur_Q
            best_class = i
            
    return best_Q, best_class


def make_iteration():
    global variances_table
    global MAX_TERMS
    
    rules = choose_rules_smart(MAX_TERMS)
    
    if rules is None:
        return None
    
    best_rules = prune_rules(rules)
    best_objects = get_objects(best_rules)
    best_var_gain = get_var_gain(best_objects)
    
    rules_quality = count_rules_quality(best_rules)
    
    df = pd.DataFrame(best_rules)
    df['delta'] = -rules_quality[0]
    
    merged = variances_table.merge(df, how='left')
    merged.fillna(0, inplace=True)
    
    merged['delta'] *= merged['pheromone']
    merged['delta_2'] = merged['delta'].sum() / (len(merged) - len(best_rules))
    merged['delta_2'] = merged[['delta', 'delta_2']].min(axis=1)
    merged['delta_2'] -= 2 * merged['delta']
    
    merged['pheromone'] += merged['delta_2']
    merged['prob'] = merged.pheromone * merged.quality / sum(merged.pheromone * merged.quality)
    
    variances_table = merged.drop(['delta', 'delta_2'], axis=1)
    
    return best_rules, best_var_gain


def get_sign(num):
    if num == 0:
        return '='
    elif num == 1:
        return '<'
    else:
        return '>'

if __name__ == '__main__':
    all_ok = False
    while True:
        path = input('Enter path to file with cannot clusters or "None" without quotes:\n')
        if path == 'None':
            break
        else:
            try:
                with open(path, 'r') as f:
                    cannot_clusters = [x.strip() for x in ' '.join(f.readlines()).split(',')]
            except:
                print('File doesn\'t exist, please retry.\n')
                continue
            try:
                res = [int(x) for x in cannot_clusters]
            except:
                print('Cannot clusters don\'t satisfy the given format, please update file and retry.\n')
                continue
        while True:
            thres = input('Enter the inclusion threshold for cannot clusters:\n')
            try:
                threshold = float(thres)
            except:
                print('Inclusion threshold is not a valid float number. Please retry.\n')
                continue

            if threshold < 0 or threshold > 1:
                print('Inclusion threshold should be between 0 and 1. Please retry.\n')
                continue

            all_ok = True
            break

        if all_ok:
            break

    variances_table_unclean = pd.read_csv('variances_new.csv', header=None, sep=' ')
    variances_table_list = []

    for i, k in enumerate(variances_table_unclean.values):
        try:
            _ = float(k[1])
            _ = float(k[2])
            _ = float(k[3])
            variances_table_list.append(k)
        except:
            pass

    variances_table = pd.DataFrame(variances_table_list, columns = variances_table_unclean.columns)

    for k in range(1, 4):
        variances_table[k] = variances_table[k].astype(float)

    variances_table['quality'] = (-variances_table[3] + variances_table[3].max() + variances_table[3].min()) / (variances_table[3].max() + variances_table[3].min())
    variances_table['pheromone'] = 1.0 / len(variances_table[variances_table[2] == 0])
    variances_table['prob'] = variances_table.pheromone * variances_table.quality / sum(variances_table.pheromone * variances_table.quality)

    vals = pd.read_csv('values_new.txt', sep=' ', header=None)
    columns = pd.read_csv('columns_new.txt', header=None)
    vals.columns = columns[0]

    genes = pd.read_csv('genes_new.txt', header=None)
    genes.columns = ['gene']

    vals_total = pd.concat([genes, vals], axis=1)

    hier = pd.read_csv('hierarchy_new.txt', header=None)
    hier.columns = ['hierarchy']

    hier_total = pd.concat([columns, hier], axis=1)
    hier_total['hierarchy'] = hier_total.hierarchy.apply(lambda x: x.split(' ')[1:])

    rules_all = []

    iter_tot = 1

    while True:
        if len(hier_total) <= MAX_UNCOVERED:
            break
        
        best_gain = 0
        best_rule = None
        cnt = 0
        
        convergence = False

        for x in range(ANT_COLONY_SIZE):
            print('iteration ' + str(x))
            it = make_iteration()
            if it is None:
                continue
            if it[1] > best_gain:
                best_gain = it[1]
                best_rule = it[0]
                cnt = 0
            else:
                cnt += 1
                if cnt >= CONVERGENCE:
                    convergence = True
                    break
                
        if best_rule is None:
            break
        
        objects = get_objects(best_rule)
        covered_objects.extend(objects)
        rules_all.append((best_rule, convergence))
        
        hier_total = hier_total[~hier_total[0].isin(covered_objects)]
        print('Rule #' + str(iter_tot) + ' defined.')
        iter_tot += 1

    covered_objects = []
    hier_total = pd.concat([columns, hier], axis=1)
    hier_total['hierarchy'] = hier_total.hierarchy.apply(lambda x: x.split(' ')[1:])
    hier_total_copy = hier_total.copy(deep=True)

    cnt = 0
    for k in hier_total_copy.values:
        if len(set(k[1]).intersection(set(cannot_clusters))):
            cnt += 1

    with open('final_rules.txt', 'w') as f:
        for i, (rule, conv) in enumerate(rules_all):
            objects = get_objects(rule)
            output = get_output(objects)
            f.write('Rule #' + str(i+1) + ':\n')
            if conv:
                f.write('(converged)\n')
            else:
                f.write('(max iterations)\n')
            rule_str = 'IF ' + rule[0][0] + ' ' + get_sign(rule[0][2]) + ' ' + str(rule[0][1])
            for x in rule[1:]:
                rule_str += ' AND ' + x[0] + ' ' + get_sign(x[2]) + ' ' + str(x[1])
            rule_str += ' THEN:\n'
            
            covered_objects.extend(objects)
            hier_total = hier_total[~hier_total[0].isin(covered_objects)]
            f.write(rule_str)
            cont = False
            f.write('Can clusters:\n')
            for x in output:
                if int(x[1]) not in [int(k) for k in cannot_clusters]:
                    f.write('Cluster ' + x[1] + ': ' + str(int(x[0] * 1000) * 0.1) + '% covered\n')
                else:
                    cont = True
            if len(cannot_clusters) > 0:
                f.write('...............\nCannot clusters:\n')
                for x in output:
                    if int(x[1]) in [int(k) for k in cannot_clusters]:
                        f.write('Cluster ' + x[1] + ': ' + str(int(x[0] * 1000) * 0.1) + '% covered\n')
                for k in rule:
                    f.write('----------\n')
                    for step in range(10, 110, 10):
                        mul = 1.0
                        if k[2] == 1:
                            mul *= step
                        elif k[2] == 2:
                            mul /= step
                        objects = get_objects([[k[0], k[1] * mul, (3 - k[2]) % 3]])
                        objects_in = hier_total_copy[hier_total_copy[0].isin(objects)]
                        cur_cnt = 0
                        for kk in objects_in.values:
                            if len(set(kk[1]).intersection(set(cannot_clusters))):
                                cur_cnt += 1
                        prob = 0
                        for t in output:
                            if int(t[1]) in cannot_clusters:
                                prob += t[0]
                        f.write('Scale = ' + str(step) + ', ' + k[0] + ' ' + get_sign(3 - k[2]) + ' ' + str(k[1] * mul) + ': ' + str(int(cur_cnt * 1000.0 / cnt) * 0.1) + '% covered in cannot clusters\n')
                    f.write('----------\n')
            f.write('\n')
    print('Done. ' + str(len(hier_total)) + ' of ' + str(len(hier_total_copy)) + ' objects out of signature.\n')
