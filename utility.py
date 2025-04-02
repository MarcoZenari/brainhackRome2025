import os
import glob
import h5py
import pandas as pd
from functools import reduce

def load_fmri(data_root, file_pattern):
    data = []
    subjects = []
    folders = sorted([
        f for f in os.listdir(data_root)
        if os.path.isdir(os.path.join(data_root, f))
    ])
    
    for folder in folders:
        matched_files = glob.glob(os.path.join(data_root, folder, file_pattern))
        if not matched_files:
            print(f"No files found matching {file_pattern} in {folder}")
            continue
            
        for file_path in matched_files:
            try:
                with h5py.File(file_path, 'r') as hdf:
                    data.append(hdf["dataset"][:])
                    subjects.append(folder)
            except Exception as e:
                print(f"Error loading {file_path}: {str(e)}")
                
    return data, subjects

def load_phenotype(data_root):
    # Get all non-description TSV files
    files = sorted([
        f for f in os.listdir(data_root)
        if f.endswith('.tsv') 
           and '_definitions' not in f
           and 'notes' not in f
           and 'stroop' not in f
           and 'hammer' not in f
           and os.path.isfile(os.path.join(data_root, f))
    ])
    
    # Read all dataframes
    dfs = []
    for file in files:
        file_path = os.path.join(data_root, file)

        # Skip the first line if the file is named demos.tsv
        if file == 'demos.tsv':
            df = pd.read_csv(file_path, sep=',', skiprows=1, encoding='latin1')
        else:
            df = pd.read_csv(file_path, sep='\t')

        if 'subjectkey' not in df.columns:
            raise ValueError(f"File {file} missing 'subjectkey' column")
        dfs.append(df.set_index('subjectkey'))
    
    # Merge all dataframes horizontally
    merged_df = reduce(
        lambda left, right: pd.merge(
            left, right, 
            left_index=True, 
            right_index=True, 
            how='outer', 
            suffixes=('', f'_{right.columns[0]}')  # Add suffix for duplicate columns
        ),
        dfs
    )
    
    return merged_df.reset_index()


def load_definitions(data_root):
    # Get all definition TSV files
    files = sorted([
        f for f in os.listdir(data_root)
        if f.endswith('.tsv') 
        and '_definitions' in f
        and os.path.isfile(os.path.join(data_root, f))
    ])
    
    if not files:
        return pd.DataFrame()  # Return empty dataframe if no files found
    
    dfs = []
    
    for file in files:
        file_path = os.path.join(data_root, file)
        df = pd.read_csv(file_path, sep='\t', header=0)
        dfs.append(df)
    
    # Concatenate all dataframes vertically
    concatenated_df = pd.concat(dfs, axis=0, ignore_index=True)  # Ignore original indices
    
    return concatenated_df


        
