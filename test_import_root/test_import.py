import pandas as pd
import uproot
#import uproot3
import numpy as np
import matplotlib.pyplot as plt

br_list_data = ['n_hits', 'tofpet_id', 'tofpet_channel', 'timestamp', 't_coarse', 't_fine', 'timestamp', 'v_coarse', 'v_fine', 'value']
br_list_evt = ['timestamp', 'evt_number', 'flags']
evt_tree = 'event'
hits_tree = 'board_57'

data = uproot.open(r'C:\Users\kimyk\OneDrive\Bureau\Master 1\Projet_LPHE_I\ecal_reco\test_import_root\run_000002\data_0000.root')
#print(data['event_data'].keys())


x = data['event_data']['evt_number']
y = data['event_data']['n_hits']
print(y)
print(x)
fig, ax = plt.subplots()
ax.plot(x,'+')

ax.set(xlabel='x', ylabel='voltage (mV)',
       title='test')
ax.grid()

fig.savefig("test.png")
plt.show()