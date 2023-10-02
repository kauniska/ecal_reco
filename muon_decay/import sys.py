import sys
sys.path.insert(1, r"C:\Users\kimyk\Projet LPHE I\ecal_reco\utils")
sys.path.insert(1, r"C:\Users\kimyk\Projet LPHE I\ecal_reco\tracking")
import pandas as pd
import matplotlib.pyplot as plt
import os
#from hit import Hit
raw_data_directory = r"C:\Users\kimyk\OneDrive\Bureau\Master 1\Projet_LPHE_I\ecal_reco\muon_decay\extracted_data\pickle" #path to the ecal data
current_directory = os.getcwd()
print(current_directory, '\n')
data_storage = current_directory+ "\extracted_data\pickle\\"
run = []


# Create a list of all the values (energies) of the hits made by electrons
import functools
import operator
import fnmatch

if len(run) == 0:
    runs = fnmatch.filter(os.listdir(raw_data_directory), '*')

print(runs, '\n\n')

values = []
for run in runs:
    #decay_data = pd.read_pickle(data_storage+"pickle_decay_data_"+run)

    decay_data = pd.read_pickle(data_storage+run)
    hits_electron = decay_data['hits_electron']
    hits_electron = functools.reduce(operator.iconcat, hits_electron, [])
    for hit in hits_electron:
        values.append(hit.value)


# Histogram of values of electron hits
plt.hist(values,40,density=True)
plt.xlabel("value")
plt.grid(True)
plt.show()