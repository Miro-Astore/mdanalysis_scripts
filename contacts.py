#####
#USAGE: python mdanalysis_scripts/contacts.py TOPFILE traj sel1 sel2 out_name
#
#TODO better detection of lysine bridges, trying to fix argenine detection
#####


import MDAnalysis as mda
import os 
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis import contacts

def contacts_within_cutoff(u, group_a, group_b, radius=8.0):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])

    return np.array(timeseries)

print ('TOP ' + str(sys.argv[1]))
print ('traj ' + str(sys.argv[2]))
print ('sel1 ' + str(sys.argv[3]))
print ('sel2 ' + str(sys.argv[4]))
print ('out_name ' + str(sys.argv[5]))

u=mda.Universe(sys.argv[1],sys.argv[2])

file_out_name = str(sys.argv[5])

sel1txt=str(sys.argv[3])
sel2txt=str(sys.argv[4])

sel_1_atoms =  u.select_atoms (sel1txt)
sel_2_atoms =  u.select_atoms (sel2txt)

ca = contacts_within_cutoff (u, sel_1_atoms, sel_2_atoms,radius=15)

np.savetxt(str(sys.argv[5]),ca)


print(ca)

#ca_df = pd.DataFrame(ca, columns=['Frame', '# Contacts']) 

#ca_df.plot(x='Frame')
#plt.ylabel('# salt bridges')

## reference groups (first frame of the trajectory, but you could also use a
## separate PDB, eg crystal structure)
#acidic_sel_1 = sel_1_atoms.select_atoms(sel_acidic)
#basic_sel_1 = sel_1_atoms.select_atoms(sel_basic)
#
#acidic_sel_2 = sel_2_atoms.select_atoms(sel_acidic)
#basic_sel_2 = sel_2_atoms.select_atoms(sel_basic)
#
## set up analysis of native contacts ("salt bridges"); salt bridges have a
## distance <6 A
#ca1 = contacts.Contacts(u, select=(sel_acidic, sel_basic), radius=4.5)
#                        # iterate through trajectory and perform analysis of "native contacts" Q
#ca1.run()
## print number of averave contacts
#average_contacts = np.mean(ca1.timeseries[:, 1])
#print('average contacts = {}'.format(average_contacts))
## plot time series q(t)
#fig, ax = plt.subplots()
#ax.plot(ca1.timeseries[:, 0], ca1.timeseries[:, 1])
#ax.set(xlabel='frame', ylabel='fraction of native contacts',
#title='Native Contacts, average = {:.2f}'.format(average_contacts))
#fig.show()
#
#
#
