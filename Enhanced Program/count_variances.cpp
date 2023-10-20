#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

unordered_map<string, int> genes_dict;
unordered_map<string, int> files_dict;

vector<string> genes_list;
vector<string> files_list;

vector<string> values_matrix_str;
vector<vector<int> > values_matrix;
vector<vector<int> > unique_values;

vector<vector<int> > hierarchy;
vector<pair<pair<int, int>, int> > rules;

vector<int> rule_elements;

void read_genes_file(string genes_file) {
    ifstream myfile(genes_file.c_str());
    string line;
    if (myfile.is_open()) {
        int cnt = 0;
        while (getline(myfile, line)) {
            genes_dict[line] = cnt++;
            genes_list.push_back(line);
        }
        myfile.close();
    }
}

void read_columns_file(string columns_file) {
    ifstream myfile(columns_file.c_str());
    string line;
    if (myfile.is_open()) {
        int cnt = 0;
        while (getline(myfile, line)) {
            files_dict[line] = cnt++;
            files_list.push_back(line);
        }
        myfile.close();
    }
}

void read_values_file(string values_file) {
    ifstream myfile(values_file.c_str());
    string line;
    if (myfile.is_open()) {
        values_matrix.resize(genes_list.size());
        for (int i=0;i<genes_list.size();++i) {
            getline(myfile, line);
            values_matrix_str.push_back(line);
        }
        for (int i=0;i<genes_list.size();++i) {
            stringstream temp(values_matrix_str[i]);
            values_matrix[i].resize(files_list.size());
            for (int j=0;j<files_list.size();++j) {
                temp >> values_matrix[i][j];
            }
        }
        myfile.close();
    }
}

void read_hierarchy_file(string hierarchy_file) {
    ifstream myfile(hierarchy_file.c_str());
    string line;
    hierarchy.resize(files_list.size());
    
    if (myfile.is_open()) {
        for (int i=0;i<files_list.size();++i) {
            int sz;
            myfile >> sz;
            hierarchy[i].resize(sz);
            for (int j=0;j<sz;++j) {
                myfile >> hierarchy[i][j];
            }
        }
        myfile.close();
    }
}  

void get_unique_values() {
    for (int i=0;i<values_matrix.size();++i) {
        set<int> s;
        for (int j=0;j<values_matrix[i].size();++j) {
            s.insert(values_matrix[i][j]);
        }
        vector<int> cur_unique_values;
        for (auto it: s) {
            cur_unique_values.push_back(it);
        }
        sort(cur_unique_values.begin(), cur_unique_values.end());
        unique_values.push_back(cur_unique_values);
    }
}

void generate_rules() {
    for (int i=0;i<unique_values.size();++i) {
        for (int j=0;j<unique_values[i].size();++j) {
            rules.push_back(make_pair(make_pair(i, unique_values[i][j]), 0));
            if (j > 1) {
                rules.push_back(make_pair(make_pair(i, unique_values[i][j]), 1));
            }
            if (j < unique_values[i].size()-2) {
                rules.push_back(make_pair(make_pair(i, unique_values[i][j]), 2));
            }
        }
    }
}

float get_rule_elements(int rule_index) {
    rule_elements.clear();
    int str_index = rules[rule_index].first.first;
    int val = rules[rule_index].first.second;
    int type = rules[rule_index].second;
    
    for (int k=0; k<values_matrix[str_index].size(); ++k) {
        if ((type == 0) && (values_matrix[str_index][k] == val)) {
            rule_elements.push_back(k);
        }
        if ((type == 1) && (values_matrix[str_index][k] < val)) {
            rule_elements.push_back(k);
        }
        if ((type == 2) && (values_matrix[str_index][k] > val)) {
            rule_elements.push_back(k);
        }
    }
    if (rule_elements.size() == 1) {
        return 0;
    }
    vector<float> elements(files_list.size());
    vector<float> cur_elements(files_list.size(), 0);
    for (int i=0;i<rule_elements.size();++i) {
        for (int k=0;k<hierarchy[i].size();++k) {
            elements[hierarchy[i][k]] += 1;
        }
    }
    for (int i=0;i<files_list.size();++i) {
        elements[i] /= rule_elements.size();
    }
    float total_sq = 0;
    for (int i=0;i<rule_elements.size();++i) {
        float cur_sq = 0;
        for (int k=0;k<hierarchy[i].size();++k) {
            cur_elements[hierarchy[i][k]] += 1;
        }
        for (int k=0;k<files_list.size();++k) {
            float cur_weight = 0.01 + 0.98 * (files_list.size() - k - 1) / (files_list.size() - 1);
            if (cur_elements[k] != elements[k]) {
                cur_sq += cur_weight * (cur_elements[k] - elements[k]) * (cur_elements[k] - elements[k]);
            }
        }
        for (int k=0;k<hierarchy[i].size();++k) {
            cur_elements[hierarchy[i][k]] -= 1;
        }
        total_sq += cur_sq;
    }
    return total_sq / rule_elements.size();
}

int main() {
    read_genes_file("genes_new.txt");
    read_columns_file("columns_new.txt");
    read_values_file("values_new.txt");
    read_hierarchy_file("hierarchy_new.txt");
    get_unique_values();
    generate_rules();
    
    string output_file = "variances_new.csv";
    ofstream outfile(output_file.c_str());
    if (outfile.is_open()) {
        for (int i=0;i<rules.size();++i) {
            if ((i+1) % 10000 == 0) {
                cout << i+1 << endl;
            }
            float res = get_rule_elements(i);
            outfile << genes_list[rules[i].first.first] << " " << rules[i].first.second << " " << rules[i].second << " " << res << '\n';
        }
        outfile.close();
    }
    return 0;
}