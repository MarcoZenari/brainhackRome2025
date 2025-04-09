import os
import glob
import h5py
import pandas as pd
from functools import reduce
import numpy as np
from scipy.signal import butter, filtfilt
import matplotlib.pyplot as plt

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

def filter(data, fs, cutoff, filter_type, order):

    filtered_data = []
    for bold in data:
        filtered_bold = butterworth_filter(bold, fs, cutoff, filter_type, order)
        filtered_data.append(filtered_bold)

    return filtered_data

def butterworth_filter(matrix, fs, cutoff, filter_type='high', order=3):
    """
    Applies high-pass or band-pass Butterworth filter to a parcels × time matrix.
    
    Args:
        matrix: 2D array (parcels × time steps)
        fs: Sampling frequency (Hz)
        cutoff: Cutoff frequency/frequencies:
            - Single value for high-pass (e.g., 0.01)
            - Tuple/list for band-pass (e.g., [0.01, high_cutoff])
        filter_type: 'high' or 'band'
        order: Filter order (default=3)
    
    Returns:
        Filtered matrix with same shape as input
    """
    # Input validation
    if filter_type == 'band' and not isinstance(cutoff, (list, tuple)):
        raise ValueError("For band-pass, cutoff must be a list/tuple of [low, high]")
    
    # Design filter
    nyq = 0.5 * fs
    if filter_type == 'high':
        normal_cutoff = cutoff / nyq
        b, a = butter(order, normal_cutoff, btype='high', analog=False)
    elif filter_type == 'band':
        normal_cutoff = [c / nyq for c in cutoff]
        b, a = butter(order, normal_cutoff, btype='band', analog=False)
    
    # Apply zero-phase filtering
    filtered = filtfilt(b, a, matrix, axis=1)  # Axis 1 assumes time is columns
    
    return filtered

def remove_movement(data, subjects, fd):

    good_data = []
    for n, subj in enumerate(subjects):
        subj_idx = fd.index[fd.iloc[:, 0] == subj].to_list()
        if not subj_idx:  # Subject not found case
            print(f"Warning: Subject {subj} not found in FD file. Skipping.")
            continue

        movement = fd.iloc[subj_idx, 1:]
        bold = data[n]
        no_move = np.squeeze(movement < 0.5)
        good_points = bold[:, no_move]
        good_data.append(good_points)

    return good_data

def load_events(data_root, subjects, file_pattern):
    data = []
    subj_folder = ["sub-" + sub.replace("_", "") for sub in subjects]
    folders = [os.path.join(data_root, sub, 'func') for sub in subj_folder]
    
    for folder in folders:
        matched_files = glob.glob(os.path.join(folder, file_pattern))
        if not matched_files:
            print(f"No files found matching {file_pattern} in {folder}")
            continue
            
        for file_path in matched_files:
            try:
                data.append(pd.read_csv(file_path, sep='\t'))
            except Exception as e:
                print(f"Error loading {file_path}: {str(e)}")
                
    return data

def calculate_volume_pattern(events, tr, n_volumes):

    volume_pattern = np.zeros((len(events), n_volumes))
    for nsubj, event_df in enumerate(events):
        nevents = len(event_df)
        onsets = event_df["onset"].to_numpy()
        onsets = np.floor(onsets / tr)
        durations = event_df["duration"].to_numpy()
        durations = np.floor(durations / tr)
        mapping = {value: index for index, value in enumerate(event_df['trial_type'].unique())}
        trial_types = (event_df['trial_type'].map(mapping)).to_numpy()
        for trial in range(0,nevents):
            onset = int(onsets[trial])
            duration = int(durations[trial])
            trial_type = trial_types[trial]
            volume_pattern[nsubj, onset:(onset+duration)] = trial_type + 1

    return volume_pattern

def get_info(data, column, subjects):
    infos = []
    for subj in subjects:
        subj_idx = data.index[data.iloc[:, 0] == subj].to_list()
        info = data[column].iloc[subj_idx].item()
        infos.append(info)

    return infos

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


def plot_signals(x, dt=1, n_plot=10, n_parcels=None, title=None, figsize=(10,5)):
    '''
    Plot the trajectory of a signal over time
    ---------------------------------------------- PARAMETERS ----------------------------------------------
    x (ndarray): 2D array (time_steps x n_parcels) representing the trajectory of the signal
    dt (float): Time step used in the simulation
    num_trajectories (int): Number of trajectories to randomly select for plotting
    time_steps (int): Total number of time steps in the simulation
    n_parcels (int): Total number of parcels in the simulation
    --------------------------------------------------------------------------------------------------------
    '''
    time_steps = x.shape[0]
    n_parcels = x.shape[1]

    # Generate the time array
    time = np.arange(0, time_steps * dt, dt)

    # Select random trajectories to plot
    random_indices = np.random.choice(n_parcels, size=n_plot, replace=False)
    selected_trajectories = x[:, random_indices]

    # Create the plot
    plt.figure(figsize=figsize)
    plt.plot(time, selected_trajectories, lw=1.5)
    plt.title(title, fontsize=14)
    plt.xlabel('Time', fontsize=13)
    plt.ylabel('x(t)', fontsize=13)
    plt.yticks([np.round(np.max(selected_trajectories), 0), 0, np.round(np.min(selected_trajectories), 0)])

    # Customize the plot
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout()
