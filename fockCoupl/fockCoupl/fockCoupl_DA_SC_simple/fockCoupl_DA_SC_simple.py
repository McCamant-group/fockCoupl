
def main():
	import numpy as np
	from numpy.linalg import matrix_power
	from numpy.linalg import eig
	from scipy.optimize import curve_fit
	from scipy.linalg import sqrtm
	import matplotlib.pyplot as plt
	import sys
	import warnings 
	import os
	from mpl_toolkits.mplot3d import Axes3D

	warnings.filterwarnings('ignore') # Supresses 'casting to real discards complex part' warning

	###################### Calculate a-b electronic coupling ######################

	paramsIn = sys.argv[1] # Input parameters
	f = np.loadtxt(sys.argv[2]) # Fock matrix
	s = np.loadtxt(sys.argv[3]) # Overlap matrix

	with open(paramsIn,'r') as i:
		
		inp = i.readlines()

		j=0
		for l in inp:
			l = l.partition('#')[0]
			l = l.rstrip()
			inp[j] = l
			j+=1

		dim = int(inp[0]) # total number of basis functions
		dima = int(inp[1])  # basis functions of fraction a
		dimb = int(inp[2])  # basis functions of fraction b
		b_orb = int(inp[3])
		hbar = 1.054e-27 # erg s
		T = float(inp[4]) # K
		kBT = 8.6173303e-5*T # eV
		reorg = float(inp[5]) # eV
		Efb = float(inp[6]) # eV
		Eox = float(inp[7]) # eV O2=3.343 O1=3.207
		au2ev = 27.2113966413079
		ev2erg = 1.6e-12 # eV to erg 

	i.close()

	## S^-1/2 ##
	s_sqrt = sqrtm(s)
	s_inv_sqrt = matrix_power(s_sqrt,-1)

	## F_orth = S^-1/2 * F * S^-1/2
	f_orth = np.dot(s_inv_sqrt,np.dot(f,s_inv_sqrt))

	## Diagonalize aa, bb and cc blocks of F_orth ##
	f_aa_eval,f_aa_evec = eig(f_orth[:dima,:dima])
	f_bb_eval,f_bb_evec = eig(f_orth[dima:(dima+dimb),dima:(dima+dimb)])

	## Use eigenvectors of above block diagonalization ##
	## to calculate off diagonal blocks which are the coupling ##
	f_ab = np.dot(f_aa_evec.conj().T,np.dot(f_orth[:dima,dima:(dima+dimb)],f_bb_evec))

	Ea = np.real(f_aa_eval)*au2ev#*au2ev2cm1#*au2ev
	Eb = np.real(f_bb_eval)*au2ev#*au2ev2cm1#*au2ev

	Eb_orb = Eb[b_orb]
	print(Eb_orb)

	V_ab = np.zeros(len(Ea),dtype=complex)
	V_ab = f_ab[:,b_orb]*au2ev#*au2ev2cm1#*au2ev
	V_tot = V_ab 

	j=0
	V_avg=0
	for i in range(len(Ea)):
		
		if Ea[i] < 0.0 and Ea[i] > -4.0:
			V_avg += abs(V_tot[i])
			j+=1
			
	V_avg = ev2erg*V_avg/j
	print V_avg/ev2erg


	def rate():
		
		sigma = 32.5e-24 # cm^3
		me = 15*9.109e-28 # g
		C = (((2*me)**(3./2.))*(V_avg**2)*sigma)/(np.pi*hbar**4)
		ket = C*np.sqrt(-ev2erg*(Eox-Efb)-ev2erg*reorg) #
		print C
		print ket
		print (1/ket)*1e12

	rate()


