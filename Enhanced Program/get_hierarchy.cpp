#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

using namespace std;

vector<string> filenames;
vector<vector<unsigned short> > final_clustering;
vector<pair<string, int> > clustering_str;
string line;

void read_filenames_from_folder(string filenames_file) {
    ifstream myfile(filenames_file);
    if (myfile.is_open()) {
    	while (getline(myfile,line)) {
      	    filenames.push_back(line);
    	}
        myfile.close();
    }
}

void read_clustering_iteration(string iteration_file) {
    ifstream myfile(iteration_file);
    vector<unsigned short> data;
    if (myfile.is_open()) {
        unsigned short cur_number;
    	while (myfile >> cur_number) {
            data.push_back(cur_number);
    	}
        myfile.close();
    }
    final_clustering.push_back(data);
}


int main(int argc, char *argv[]) {
    read_filenames_from_folder(string(argv[1]) + "/filenames.txt");
    for (int i=0;i<filenames.size()-1;++i) {
        read_clustering_iteration(string(argv[1]) + "/iter" + to_string(i+1) + ".txt");
    }
    vector<unsigned short> indices;
    for (int i=0;i<filenames.size();++i) {
        indices.push_back(i);
    }
    final_clustering.push_back(indices);
    for (int i=0;i<final_clustering.size();++i) {
        for (int j=0;j<i;++j) swap(final_clustering[i][j], final_clustering[j][i]);
    }
    sort(final_clustering.begin(), final_clustering.end());
    ofstream myfile(string(argv[1]) + "/total_data.csv");
    if (myfile.is_open()) {
        for (int i=0;i<final_clustering.size();++i) {
            myfile << filenames[final_clustering[i][final_clustering.size()-1]];
            for (int j=0;j<final_clustering.size()-1;++j) {
                myfile << ',';
                if (i == 0 || final_clustering[i][j] != final_clustering[i-1][j]) {
                    myfile << final_clustering[i][j];
                }
            }
            myfile << '\n';
        }
        myfile.close();
    }
    return 0;
}