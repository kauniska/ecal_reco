import pandas as pd
import uproot
#import uproot3
import numpy as np
import matplotlib.pyplot as plt

br_list_data = ['n_hits', 'tofpet_id', 'tofpet_channel', 'timestamp', 't_coarse', 't_fine', 'timestamp', 'v_coarse', 'v_fine', 'value']
br_list_evt = ['timestamp', 'evt_number', 'flags']
evt_tree = 'event'
hits_tree = 'board_57'

data = uproot.open('C:\\Users\\Pascal\\Desktop\\TP4a\\git\\test_data_loading\\data_0000.root')

x = data['event_data']['evt_number']
y = data['event_data']['n_hits']
print(y)
print(data['event_data'].keys())
print(x)
fig, ax = plt.subplots()
ax.plot(y,'+')

ax.set(xlabel='x', ylabel='voltage (mV)',
       title='test')
ax.grid()

fig.savefig("test.png")
plt.show()