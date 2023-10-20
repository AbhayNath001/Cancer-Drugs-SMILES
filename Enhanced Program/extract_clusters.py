import pandas as pd
import sys

def equals(column_1, column_2):
    if len(column_1) != len(column_2):
        return False
    for i in range(len(column_1)):
        try:
            val_1 = int(float(column_1[i]))
            val_2 = int(float(column_2[i]))
            if val_1 != val_2:
                return False
        except:
            if str(column_1[i]) != str(column_2[i]):
                return False
    return True

if __name__ == "__main__":

    data = pd.read_csv('total_data.csv', header=None)
    data.columns = [str(x) for x in data.columns]

    while True:
        cluster_num_str = input("Enter number of a cluster or 'Exit' (without quotes) to end the program:\n")
        if cluster_num_str == 'Exit':
            break
        try:
            cluster_num = int(cluster_num_str)
            if cluster_num < 0 or cluster_num >= len(data.columns) - 1:
                print("Cluster number out of range - should be between 0 and " + str(len(data.columns) - 2))
                continue
        except:
            print("Cluster number should be integer!")
            continue

        column_needed = data[str(cluster_num+1)]
        start_index, end_index = None, None
        started = False
        for i, k in enumerate(column_needed):
            if started:
                end_index = i
            if str(k) == 'nan':
                continue
            if int(k) == cluster_num:
                if start_index is None:
                    start_index = i
                started = True
            else:
                started = False
        end_index = max(end_index, start_index+1)
        data_fragment = data[start_index:end_index]
        data_needed = data_fragment[['0'] + list(data_fragment.columns[cluster_num+1:-1])]
        del data_fragment
        columns_all = [data_needed['0'], data_needed[data_needed.columns[1]]]
        for k in data_needed.columns[2:]:
            if equals(data_needed[k].values, columns_all[-1].values):
                 continue
            else:
                columns_all.append(data_needed[k])
        data_filtered = pd.DataFrame(columns_all).T
        while True:
            res = sorted([int(k) for k in columns_all[-1].values if str(k).lower() != 'nan'])
            with open("cluster_" + cluster_num_str + ".txt", 'w') as f:
                f.write(', '.join([str(k) for k in res]))
            print("Cluster " + str(cluster_num) + " has " + str(len(columns_all)-1) + " subcluster(s). They are listed in file cluster_" + cluster_num_str + ".txt")
            level_str = input("Enter number to specify the hierarchy level:\n")
            try:
                level = int(level_str)
                if (level < 1) or (level > len(data_filtered)):
                    print("Hierarchy level should be a number in range between 1 and " + str(len(columns_all)-1) + ". Please try once again.\n")
                    continue
                else:
                    values_to_output = data_filtered[[data_filtered.columns[0]] + [data_filtered.columns[level]]].values
            except:
                print("Hierarchy level should be an integer! Please try once again.\n")
                continue
            with open("cluster_" + cluster_num_str + "_hier_" + level_str + ".txt", 'w') as f:
                str_to_write = ""
                for k in values_to_output:
                    if str(k[1]) != 'nan':
                        str_to_write += 'Cluster #'+ str(int(float(k[1]))) + ':\n'
                    str_to_write += k[0] + '\n'
                f.write(str_to_write)
            yn = input("Result was written to file cluster_" + cluster_num_str + "_hier_" + level_str + ".txt\nDo you want to continue with cluster " + str(cluster_num) + "?(Y/N)")
            if yn.lower() == 'y':
                continue
            else:
                break
