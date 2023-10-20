import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs
from difflib import SequenceMatcher

# Load your dataset
df = pd.read_csv("SMILES_augmented.csv", usecols=[1])  # Use the second column containing SMILES

# Function to calculate Tanimoto similarity
def calculate_similarity(input_smiles, dataset_smiles):
    mol1 = Chem.MolFromSmiles(input_smiles)
    mol2 = Chem.MolFromSmiles(dataset_smiles)
    if mol1 is not None and mol2 is not None:
        fps1 = Chem.RDKFingerprint(mol1)
        fps2 = Chem.RDKFingerprint(mol2)
        similarity = DataStructs.TanimotoSimilarity(fps1, fps2)
        return similarity
    else:
        return 0.0

# Input SMILES
input_smiles = input('Your_Input_SMILES_Goes_Here: ')

# Threshold for similarity
threshold = 0.75

# Check similarity for each compound in the dataset
found = False
for index, row in df.iterrows():
    smiles = row['Smiles']
    similarity = calculate_similarity(input_smiles, smiles)
    
    if similarity >= threshold:
        print(f"This is a cancer-related drug (from SMILES_augmented.csv): Similarity: {similarity:.2f}")
        found = True
        break

if not found:
    # Check for similar SMILES using SequenceMatcher
    for index, row in df.iterrows():
        smiles = row['Smiles']
        similarity = SequenceMatcher(None, input_smiles, smiles).ratio()
        if similarity >= threshold:
            print(f"This is a cancer-related drug (from SMILES_augmented.csv): Similarity with other cancer drugs: {similarity:.2f}")
            break

if not found:
    print(f"Not a cancer-related drug (Similarity {similarity:.2f})")