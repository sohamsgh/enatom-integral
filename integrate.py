import numpy as np

a = np.loadtxt("dos_free.dat", skiprows=1)

energy = a[:,0]
dos_tetra = a[:,1]
dos_theory = a[:,2]
tot_tetra = np.trapz(dos_tetra[0:51], energy[0:51])
tot_theory = np.trapz(dos_theory[0:51], energy[0:51])
print (tot_tetra, tot_theory)
