import pandas as pd
import uproot
#import uproot3
import numpy as np
from fnmatch import filter
from os import listdir


br_list_data = ['evt_timestamp', 'evt_number', 'evt_flags','n_hits', 'tofpet_id', 'tofpet_channel', 'timestamp', 't_coarse', 't_fine', 'timestamp', 'v_coarse', 'v_fine', 'value']
#br_list_evt = ]
# evt_tree = 'event_data'
tree = 'event_data'


def load_dataset(file_path):
    """
    Load single .root file containing events and hits (board_x) tree
    ## events tree is global
    ## board_x tree contains hits for each event

    Extract global time from events df and add to the one with hits
    """


    with uproot.open(file_path) as tree:
        # hits_dict = tree[hits_tree].arrays(br_list_data, library="np")
        hits_dict = tree[tree].arrays(br_list_data, library="np")
        #evts_dict = tree[evt_tree].arrays(br_list_evt, library="np")
    
   # df_evts = pd.DataFrame.from_dict(evts_dict)
    df_hits = pd.DataFrame.from_dict(hits_dict)
    #df_hits['timestamp_event'] = df_evts['timestamp']
    df_hits['timestamp_event'] = df_hits['evt_timestamp']
    
    return df_hits


def load_run(run_path, n_hits_min=6, n_hits_max=50):
    """
    Arguments :
        -run_path : path to the directory of a run containing .root files (The separator "/" or "\" should already be put at the end of the path)
    Returns :
        -df_hits_total : pandas dataframe made of the concatenation of all the events in each .root file
        -df_hits : pandas dataframe with the events only containing 6 <= n_hits <= 50
        -og_len : length of df_hits_total
        -new_len : length of df_hits
    """
    files = filter(listdir(run_path),'*.root')
    for i,file in enumerate(files):
        if i == 0:
            df = [load_dataset(run_path+file)]
        else:
            df_i = load_dataset(run_path+file)
            df_i.index += df[-1].index[-1]+1
            df.append(df_i)

    df_hits_total = pd.concat(df)
    df_hits = pd.DataFrame.copy(df_hits_total, deep=True)

    og_len = len(df_hits_total)

    min_condition = "n_hits > "+str(n_hits_min)
    max_condition = "n_hits < "+str(n_hits_max)
    df_hits.query(min_condition, inplace=True)
    df_hits.query(max_condition, inplace=True)

    new_len = len(df_hits)

    print('selected {:.2f}% of all events'.format(new_len/og_len * 100))
    return df_hits_total, df_hits, og_len, new_len

