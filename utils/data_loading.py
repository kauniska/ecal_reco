import pandas as pd
import uproot
#import uproot3
import numpy as np


br_list_data = ['n_hits', 'tofpet_id', 'tofpet_channel', 'timestamp', 't_coarse', 't_fine', 'timestamp', 'v_coarse', 'v_fine', 'value']
br_list_evt = ['timestamp', 'evt_number', 'flags']
evt_tree = 'event'
hits_tree = 'board_57'


def load_dataset(file_path):
    """
    Load single .root file containing events and hits (board_x) tree
    ## events tree is global
    ## board_x tree contains hits for each event

    Extract global time from events df and add to the one with hits
    """


    with uproot.open(file_path) as tree:
        hits_dict = tree[hits_tree].arrays(br_list_data, library="np")
        evts_dict = tree[evt_tree].arrays(br_list_evt, library="np")
    
    df_evts = pd.DataFrame.from_dict(evts_dict)
    df_hits = pd.DataFrame.from_dict(hits_dict)
    df_hits['timestamp_event'] = df_evts['timestamp']

    return df_hits

