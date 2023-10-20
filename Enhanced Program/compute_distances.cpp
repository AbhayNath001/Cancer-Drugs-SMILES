#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <unordered_map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

const int MAX_GENES_COUNT = 65000;
const int MAX_OBJECT_FILE_SIZE = 2000000;

vector<string> filenames;
std::vector<std::vector<unsigned short> > expressions_matrix;
std::vector<std::vector<float> > distance_matrix;

unordered_map<string, int> mp;
vector<unsigned short> cur_clustering;
vector<unsigned short> files_count_by_folder(4);

char buf[MAX_OBJECT_FILE_SIZE];
string line;
int total_index = 0;
unsigned short max_value = 0;
unsigned short max_cluster = 0;
bool first_time = true;
bool first_time_2 = true;

vector<float> sum_dists;
vector<float> dists_to_new;
vector<unsigned short> counts;
vector<vector<unsigned short> > cur_distrib;

int get_folder_num(int file_index) {
    int second_folder_start = files_count_by_folder[0];
    int third_folder_start = files_count_by_folder[0] + files_count_by_folder[1];
    int forth_folder_start = files_count_by_folder[0] + files_count_by_folder[1] + files_count_by_folder[2];
    if (file_index < second_folder_start) return 1;
    else if (file_index < third_folder_start) return 2;
    else if (file_index < forth_folder_start) return 3;
    else return 4;
}

void update_stats() {
    sum_dists.clear();
    sum_dists.resize(filenames.size());
    counts.clear();
    counts.resize(filenames.size());
    dists_to_new.clear();
    dists_to_new.resize(filenames.size());
    for (int i = 0; i < cur_distrib.size(); ++i) cur_distrib[i].clear();
    cur_distrib.resize(max_cluster + 1);

    for (int i = 0; i < filenames.size(); ++i) {
        cur_distrib[cur_clustering[i]].push_back(i);
    }

    for (int t = 0; t < cur_distrib.size(); ++t) {
        counts[t] = cur_distrib[t].size();
        for (int i1 = 0; i1 < cur_distrib[t].size(); ++i1) {
            for (int j1 = 0; j1 < i1; ++j1) {
                int i = cur_distrib[t][i1];
                int j = cur_distrib[t][j1];
                sum_dists[i] += distance_matrix[i][j];
                sum_dists[j] += distance_matrix[i][j];
            }
        }
    }
}

void extract_element(int element_to_extract) {
    int cluster_to_take = cur_clustering[element_to_extract];
    counts[cluster_to_take]--;
    counts[max_cluster]++;
    cur_clustering[element_to_extract] = max_cluster;

    for (int i1 = 0; i1 < cur_distrib[cluster_to_take].size(); ++i1) {
        int i = cur_distrib[cluster_to_take][i1];
        if (cur_clustering[i] != cluster_to_take) continue;
        sum_dists[i] -= distance_matrix[max(i, element_to_extract)][min(i, element_to_extract)];
        dists_to_new[i] += distance_matrix[max(i, element_to_extract)][min(i, element_to_extract)];
    }
}

vector<string> read_restriction_folders(string restriction_file) {
    ifstream myfile(restriction_file);
    vector<string> restriction_folders;
    if (myfile.is_open()) {
        while (getline(myfile, line)) {
            restriction_folders.push_back(line);
        }
        myfile.close();
    }
    return restriction_folders;
}

void read_filenames_from_folder(string filenames_file, int index) {
    ifstream myfile(filenames_file);
    if (myfile.is_open()) {
        while (getline(myfile, line)) {
            filenames.push_back(line);
            files_count_by_folder[index]++;
        }
        myfile.close();
    }
}

void write_file_to_matrix(int file_index) {
    FILE* fp;
    fp = fopen(filenames[file_index].c_str(), "r");
    int cnt = fread(buf, sizeof(char), MAX_OBJECT_FILE_SIZE, fp);
    bool first_part = true;
    string cur_str = "";
    int index = -1;
    for (int k = 16; k < cnt; ++k) {
        if (buf[k] == ',') {
            if (mp.find(cur_str) == mp.end()) {
                mp[cur_str] = total_index++;
            }
            index = mp[cur_str];
            first_part = false;
            cur_str = "";
            continue;
        }
        if (buf[k] == '\n') {
            expressions_matrix[file_index][index] = atoi(cur_str.c_str());
            max_value = max(max_value, expressions_matrix[file_index][index]);
            first_part = true;
            cur_str = "";
            continue;
        }
        cur_str += buf[k];
    }
    fclose(fp);
}

void count_euclidean_distances() {
    for (int i = 0; i < filenames.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            for (int k = 0; k < total_index; ++k) {
                float cur = expressions_matrix[i][k];
                cur -= expressions_matrix[j][k];
                distance_matrix[i][j] += cur * cur;
            }
            distance_matrix[i][j] = sqrt(distance_matrix[i][j]);
        }
        cout << "All distances for object number " << i << " have been computed." << endl;
    }
}

void count_pearson_distances() {
    vector<float> avg(filenames.size());
    for (int i = 0; i < filenames.size(); ++i) {
        for (int k = 0; k < total_index; ++k) {
            avg[i] += expressions_matrix[i][k];
        }
        avg[i] /= total_index;
    }
    for (int i = 0; i < filenames.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            float sum1 = 0, sum2 = 0, sum3 = 0;
            for (int k = 0; k < total_index; ++k) {
                sum1 += (expressions_matrix[i][k] - avg[i]) * (expressions_matrix[j][k] - avg[j]);
                sum2 += (expressions_matrix[i][k] - avg[i]) * (expressions_matrix[i][k] - avg[i]);
                sum3 += (expressions_matrix[j][k] - avg[j]) * (expressions_matrix[j][k] - avg[j]);
            }
            distance_matrix[i][j] = 1 - sum1 / (sqrt(sum2) * sqrt(sum3));
        }
        cout << "All distances for object number " << i << " have been computed." << endl;
    }
}

void count_spearman_distances() {
    vector<float> avg(filenames.size());
    for (int i = 0; i < filenames.size(); ++i) {
        for (int k = 0; k < total_index; ++k) {
            avg[i] += expressions_matrix[i][k];
        }
        avg[i] /= total_index;
    }
    for (int i = 0; i < filenames.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            float sum1 = 0, sum2 = 0, sum3 = 0;
            for (int k = 0; k < total_index; ++k) {
                sum1 += (expressions_matrix[i][k] - avg[i]) * (expressions_matrix[j][k] - avg[j]);
                sum2 += (expressions_matrix[i][k] - avg[i]) * (expressions_matrix[i][k] - avg[i]);
                sum3 += (expressions_matrix[j][k] - avg[j]) * (expressions_matrix[j][k] - avg[j]);
            }
            distance_matrix[i][j] = 1 - sum1 / (sqrt(sum2) * sqrt(sum3));
        }
    }
    count_pearson_distances();
}

void count_manhattan_distances() {
    for (int i = 0; i < filenames.size(); ++i) {
        for (int j = 0; j < i; ++j) {
            for (int k = 0; k < total_index; ++k) {
                distance_matrix[i][j] += abs(expressions_matrix[i][k] - expressions_matrix[j][k]);
            }
        }
        cout << "All distances for object number " << i << " have been computed." << endl;
    }
}

void create_new_cluster() {

    update_stats();
    float max_dist = -1;
    int cluster_to_take = -1;

    for (int t = 0; t < cur_distrib.size(); ++t) {
        for (int i1 = 0; i1 < cur_distrib[t].size(); ++i1) {
            for (int j1 = 0; j1 < i1; ++j1) {
                int i = cur_distrib[t][i1];
                int j = cur_distrib[t][j1];
                if (distance_matrix[i][j] > max_dist) {
                    max_dist = distance_matrix[i][j];
                    cluster_to_take = t;
                }
            }
        }
    }

    bool is_from_first_folder = false;
    bool is_from_second_folder = false;

    for (int i1 = 0; i1 < cur_distrib[cluster_to_take].size(); ++i1) {
        int x = get_folder_num(cur_distrib[cluster_to_take][i1]);
        if (x == 1) is_from_first_folder = true;
        if (x == 2) is_from_second_folder = true;
    }

    bool constraint = false;
    if (is_from_first_folder && is_from_second_folder) constraint = true;

    float max_sum_dist = -1;
    int element_to_extract = -1;
    for (int i1 = 0; i1 < cur_distrib[cluster_to_take].size(); ++i1) {
        int i = cur_distrib[cluster_to_take][i1];
        if (constraint && (get_folder_num(i) > 2)) continue;
        if (sum_dists[i] > max_sum_dist) {
            max_sum_dist = sum_dists[i];
            element_to_extract = i;
        }
    }

    max_cluster++;

    if (get_folder_num(element_to_extract) == 3 && first_time) {
        for (int i1 = 0; i1 < cur_distrib[cluster_to_take].size(); ++i1) {
            int i = cur_distrib[cluster_to_take][i1];
            if (i == element_to_extract) continue;
            if (get_folder_num(i) == 3) {
                extract_element(i);
            }
        }
        update_stats();
        first_time = false;
    }

    if (get_folder_num(element_to_extract) <= 2 && first_time_2) {
        for (int i1 = 0; i1 < cur_distrib[cluster_to_take].size(); ++i1) {
            int i = cur_distrib[cluster_to_take][i1];
            if (i == element_to_extract) continue;
            if (get_folder_num(i) == get_folder_num(element_to_extract)) {
                extract_element(i);
            }
        }
        update_stats();
        first_time_2 = false;
    }

    int element_to_extract_old = element_to_extract;
    dists_to_new.resize(filenames.size());
    while (element_to_extract >= 0) {
        extract_element(element_to_extract);
        float max_diff = 0;
        element_to_extract = -1;

        for (int i1 = 0; i1 < cur_distrib[cluster_to_take].size(); ++i1) {
            int i = cur_distrib[cluster_to_take][i1];
            if (cur_clustering[i] != cluster_to_take) continue;
            if (get_folder_num(i) + get_folder_num(element_to_extract_old) == 3) continue;
            float cur_diff = sum_dists[i] / counts[cluster_to_take] - dists_to_new[i] / counts[max_cluster];

            if (cur_diff > max_diff) {
                max_diff = cur_diff;
                element_to_extract = i;
            }
        }
    }
}

string data_folder = "C:\\Users\\user\\Documents\\science\\pipeline\\sandbox\\fortesting\\";
string distance_type = "2";

int main(int argc, char* argv[]) {
    vector<string> prohibitive = read_restriction_folders(data_folder + "prohibitive_restrictions.txt");
    if (prohibitive.size() == 2) {
        read_filenames_from_folder(data_folder + prohibitive[0] + "_parsed_files.txt", 0);
        read_filenames_from_folder(data_folder + prohibitive[1] + "_parsed_files.txt", 1);
    }

    vector<string> aggregating = read_restriction_folders(data_folder + "aggregating_restrictions.txt");

    if (aggregating.size() == 1) {
        read_filenames_from_folder(data_folder + aggregating[0] + "_parsed_files.txt", 2);
    }

    vector<string> random = read_restriction_folders(data_folder + "no_restrictions.txt");

    if (random.size() == 1) {
        read_filenames_from_folder(data_folder + random[0] + "_parsed_files.txt", 3);
    }

    expressions_matrix.resize(filenames.size());
    for (int i = 0; i < expressions_matrix.size(); ++i) {
        expressions_matrix[i].resize(MAX_GENES_COUNT);
    }
    for (int i = 0; i < filenames.size(); ++i) {
        write_file_to_matrix(i);
        cout << "Object number " << i << " has been parsed." << endl;
    }
    distance_matrix.resize(filenames.size());
    for (int i = 0; i < filenames.size(); ++i) {
        distance_matrix[i].resize(i);
    }

    if (distance_type[0] == '1') {
        count_euclidean_distances();
    }
    else if (distance_type[0] == '2') {
        count_pearson_distances();
    }
    else if (distance_type[0] == '3') {
        count_spearman_distances();
    }
    else if (distance_type[0] == '4') {
        count_manhattan_distances();
    }
    cur_clustering.resize(filenames.size());
    while (max_cluster < filenames.size() - 1) {
        create_new_cluster();
        ofstream myfile;
        myfile.open(data_folder + "clustering_results//" + "iter" + to_string(max_cluster) + ".txt");
        for (int i = 0; i < cur_clustering.size(); ++i) {
            myfile << cur_clustering[i] << " ";
        }
        myfile.close();
        cout << "Clustering iteration number " << max_cluster << " completed." << endl;
    }
    ofstream myfile;
    myfile.open(data_folder + "clustering_results//" + "filenames.txt");
    for (int i = 0; i < filenames.size(); ++i) {
        myfile << filenames[i] << endl;
    }
    myfile.close();
    return 0;
}