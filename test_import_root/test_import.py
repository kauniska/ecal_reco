import pandas as pd
import uproot
#import uproot3
import numpy as np

br_list_data = ['n_hits', 'tofpet_id', 'tofpet_channel', 'timestamp', 't_coarse', 't_fine', 'timestamp', 'v_coarse', 'v_fine', 'value']
br_list_evt = ['timestamp', 'evt_number', 'flags']
evt_tree = 'event'
hits_tree = 'board_57'

with uproot.open('C:\\Users\\Pascal\\Desktop\\TP4a\\git\\test_data_loading\\data_0000.root') as h:
    h['event_data'].plot()