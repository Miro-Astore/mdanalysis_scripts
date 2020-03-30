import MDAnalysis as mda
import re 
import matplotlib.pyplot as plt
import numpy as np 


data=np.loadtxt('aligned/310K_R352Q_3_snaps_0_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='1')

data=np.loadtxt('aligned/310K_R352Q_3_snaps_1_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='2')


data=np.loadtxt('aligned/310K_R352Q_3_snaps_2_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='3')

data=np.loadtxt('aligned/310K_R352Q_3_snaps_3_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='4')

data=np.loadtxt('aligned/310K_R352Q_3_snaps_4_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='5')

plt.savefig('R352Q.pdf')

plt.figure()


data=np.loadtxt('aligned/310K_wt_2_snaps_0_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='1')

data=np.loadtxt('aligned/310K_wt_2_snaps_1_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='2')

data=np.loadtxt('aligned/310K_wt_2_snaps_2_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='3')

data=np.loadtxt('aligned/310K_wt_2_snaps_3_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='4')

data=np.loadtxt('aligned/310K_wt_2_snaps_4_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='5')

plt.savefig('WT.pdf')

plt.figure()

data=np.loadtxt('aligned/350K_I37R_2_snaps_0_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='1')

data=np.loadtxt('aligned/350K_I37R_2_snaps_1_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad)
plt.plot(z,rad,label='2')

data=np.loadtxt('aligned/350K_I37R_2_snaps_2_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='3')

data=np.loadtxt('aligned/350K_I37R_2_snaps_3_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='4')

data=np.loadtxt('aligned/350K_I37R_2_snaps_4_aligned.sph',dtype=str)
rad=np.array(data[:,-1],dtype=float)
z=np.array(data[:,8],dtype=float)[rad<12]
rad=rad[rad<12]
z=np.array(range(len(rad)))*0.2
plt.plot(z,rad,label='5')
plt.legend()

plt.xlabel('distance into pore (z axis $\AA$)')
plt.ylabel('pore radius $\AA$')

plt.savefig('I37R.pdf')
