I. General.
This software package performs hierarchical clustering of samples based on transcriptomic (or any other numerical) data and subsequently extracts signatures from the identified clusters. Clustering can be executed with specified constraints (i.e., must-links and cannot-links). The package is built upon the Constrained–Divisive and Tjala algorithms, written in C++ and Python, and can be operated from the command prompt. The input data can be provided in various file formats such as:
- text formats for tabular data (.tsv, .csv, .dat, .counts)
- electronic spreadsheet formats (.xlsx...)
- data exchange formats (.json...)
- statistical formats (.dta, .sav, .por)
- data storage formats (.h5, .parquet, .orc).

II. Requirements.
Visual Studio, MinGW, and Python are installed.

III. Getting Started.
On the first run, the .cpp files need to be compiled. One possible option is as follows:
1. Open the Visual Studio 2022 folder from the Start Menu and launch the x64 Native Tools Command Prompt Line.
2. Enter the command "cl."
3. Navigate to the folder containing the files to be compiled.
4. Compile the files using the command "cl /EHsc <file_name>."

IV. Clustering (clustering.py, make_clustering.cpp, compute_distances.cpp, get_hierarchy.cpp, build_dendrogram.py).
1. Navigate to the folder containing the clustering code and request information about the order numbers of distance metrics using the command "python clustering.py."
2. Parse the input data using the command "python clustering.py <order_number_of_the_chosen_distance_metrics>". This command requests to specify the paths to folders containing the input data, parses the input data, and creates an empty "clustering_results" folder for storing future results. This folder is created in the parent folder where the input data is located.
The order number of the distance metrics:
1 - Euclidean distance
2 - Pearson metric
3 - Spearman metric
4 - Manhattan distance
NB: The Pearson metric is likely the most suitable standard metric for analyzing transcriptomic Big Data.
3. Perform clustering using the command "make_clustering.exe clustering_results/ <order_number_of_the_chosen_distance_metrics>". This command creates result files for each iteration in the "clustering_results" folder.
4. Display the hierarchy of the created clusters using the command "get_hierarchy.exe clustering_results." This command creates a file in the "clustering_results" folder showing which clusters each sample belongs to.
5. Optionally, visualize the clustering results as a dendrogram using the command "python build_dendrogram.py" without parameters.

V. Signature Extraction (count_variances.py, count_variances.exe, extract_clusters.py, genes.txt, tjala.py, samples_extractor.py)
1. Move the aforementioned files to the "clustering_results" folder.
2. Navigate to the "clustering_results" folder.
3. Specify the list of genes whose transcription level has to be used to generate signatures using one of the following commands:
- "python count_variances.py None" to operate on the full list of genes
- "python count_variances.py genes.txt" when utilizing a specific subgroup of gene.
The "genes.txt" file should be in UTF-8 format.
4. Calculate the variances using the command "count_variances.exe total_data.csv."
5. Generate a list of subclusters and the samples comprising each cluster using the command "python extract_clusters.py. Here, the "exit" command terminates only the "extract_clusters.py" command, not the entire pipeline.
6. Generate signatures using the command "python tjala.py." Here is the option to specify cannot-clusters, which are clusters that should not be covered by the signatures (e.g., non-cancer clusters when generating signatures for cancer cells). The list of cannot-clusters can be generated using the "extract_clusters.py" command. If a cannot-cluster is specified, the command requests to specify a threshold indicating the maximum proportion of samples in the cannot cluster that may be covered by the signatures. The command generates a file named "final_rules.txt" containing the signatures. For each signature, it provides the proportion of samples covered by it in each cluster. For each term of the signature, it provides a step scale indicating the proportion of the cannot-cluster that would be excluded with each change in its value. This scale ranges from 10-fold changes to 100-fold changes with 10 equal steps. The scale is provided to assess the potential for drug development based on the given signature.
7. Obtain a list of samples that comply with the selected signature using the command "python samples_extractor.py." It creates a file with the list of samples in the "clustering_results" folder.
