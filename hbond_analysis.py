#####
#USAGE: python mdanalysis_scripts/distance.py TOPFILE traj sel1 sel2 out_name
#
#TODO better detection of lysine bridges, trying to fix argenine detection
#####

import MDAnalysis as mda
import os 
import sys
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

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

sel_donor = "(resname ARG and name HH* HE) or (resname ASN and name HD*) or (resname GLN and name HE*) or (resname HSD and name HD*) or (resname LYS and name HZ*) or (resname SER and name HG*) or (resname THR and name HG1) or (resname TRP and name HE1) or (resname TYR and name HH)"

sel_acceptor = "(resname ASN and name OD1) or (resname ASP and name OD*) or (resname GLU and name OE*) or (resname GLN and name OE*) or (resname HSD and name NE*) or (resname SER and name OG*) or (resname THR and name OG1) or (resname TYR and name OH)" 
#sel_acceptor = "name O" 

sel_donor_1 = '( ' + sel_donor + ') and (' + sel1txt + ')' 
sel_acceptor_1 = '( ' + sel_acceptor + ') and (' + sel1txt + ')' 
sel_donor_2 = '( ' + sel_donor + ') and (' + sel2txt + ')' 
sel_acceptor_2 = '( ' + sel_acceptor + ') and (' + sel2txt + ')' 

hbonds_1 = HBA(universe=u, hydrogens_sel=sel_donor_1, acceptors_sel=sel_acceptor_1)
hbonds_2 = HBA(universe=u, hydrogens_sel=sel_donor_2, acceptors_sel=sel_acceptor_2)

hbonds_1.run()
hbonds_2.run()

print (hbonds_1.results)
print (hbonds_2.results)

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
