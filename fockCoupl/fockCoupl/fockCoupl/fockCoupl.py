def main():

	import numpy as np
	from numpy.linalg import matrix_power
	from numpy.linalg import eig
	from scipy.linalg import sqrtm
	import matplotlib.pyplot as plt
	import warnings 

	warnings.filterwarnings('ignore') # Supresses 'casting to real discards complex part' warning

	###################### Calculate a-b electronic coupling ######################

	dim = 1291 # total number of basis functions
	dima = 986 # basis functions of fraction a
	dimb = 305 # basis functions of fraction b

	aorb = 347 # Orbital of interest fraction a
	borb = 103 # Orbital of interest fraction b 

	aorb = aorb-1 # because Python counts from zero? shouldn't need this
	borb = borb-1

	f = np.loadtxt("fock.dat") # read Fock matrix
	s = np.loadtxt("overlap.dat") # read overlap matrix

	## S^-1/2 ##
	s_sqrt = sqrtm(s)
	s_inv_sqrt = matrix_power(s_sqrt,-1)

	## F_orth = S^-1/2 * F * S^-1/2
	f_orth = np.dot(s_inv_sqrt,np.dot(f,s_inv_sqrt))

	## Diagonalize aa and bb blocks of F_orth ##
	f_aa_eval,f_aa_evec = eig(f_orth[:dima,:dima])
	f_bb_eval,f_bb_evec = eig(f_orth[dima:dim,dima:dim])

	## Use eigenvectors of above block diagonalization ##
	## to calculate off diagonal blocks which are the coupling ##
	f_ab = np.dot(f_aa_evec.conj().T,np.dot(f_orth[:dima,dima:dim],f_bb_evec))

	#############################################################################
	######################### Calculate the Marcus eT rate ######################
	########### for the donor orbital coupled to all acceptor orbitals ##########
	#############################################################################

	hbar = 4.135667662e-3 # eV*ps
	T = 298.0 # K
	kBT = 8.6173303e-5*T # eV
	reorg = 1500./8100. #8100cm-1/eV
	donor_orb = 169 # the orbital numbers get reordered after diagonalization, need to solve this. Just check output for now
	V_ab = f_ab[:,donor_orb]*27.2113966413079
	p_Ea = np.ones(len(np.real(f_aa_eval)))
	Ea = np.real(f_aa_eval)*27.2113966413079
	Eb = np.real(f_bb_eval)*27.2113966413079
	Eb_orb = sorted(Eb)[borb]
	time = np.arange(0,1000,1)
	k = np.zeros(len(p_Ea))

	## individual rates ##
	for i in range(len(V_ab)):
		k[i] = ((V_ab[i]*np.conj(V_ab[i]))/hbar)*np.sqrt(np.pi/(reorg*kBT))*np.exp(-((reorg+(Ea[i]-Eb_orb))**2/(4*reorg*kBT)))

	k_tot = np.sum(k)
	tau = 1./k_tot

	## Dynamics with total rate ##
	eTdynD = np.exp(-k_tot*time)
	eTdynA = 1.0 - eTdynD

	## Print some relavant output data ##
	with open("output.txt",'w') as o:
		o.write("Results: \n")
		o.write("Donor orbital energy = "), o.write(str(Eb_orb)),o.write(" eV \n")
		o.write("Rate Constant = "), o.write(str(k_tot)), o.write(" ps^-1 \n")
		o.write("Time Constant = "), o.write(str(tau)), o.write(" ps \n")
	o.close()

	## Plot the results for preliminary inspection ##
	plt.subplot(2,2,1)
	plt.subplots_adjust(hspace=0.5,wspace=0.5)
	plt.title('D-A Electronic Coupling')
	plt.xlabel('Orbital Energy (eV)')
	plt.ylabel('V_DA (meV)')
	plt.xlim(-9,0)
	plt.ylim(0,15)
	plt.vlines(Ea,0,V_ab*1e3)
	plt.vlines(Eb_orb,0,15,colors='r')
	plt.legend(['V_DA','D energy'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.subplot(2,2,2)
	plt.title('D and A Block Eigenvalues')
	plt.xlabel('Orbital Energy (eV)')
	plt.ylim(0,2)
	plt.xlim(-9,0)
	plt.vlines(Ea,0,1,colors='b')
	plt.vlines(Eb,1,2,colors='r')
	plt.legend(['A','D'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.subplot(2,2,3)
	plt.title('Marcus eT Rate')
	plt.xlabel('Orbital Energy (eV)')
	plt.ylabel('Rate (ps^-1)')
	plt.xlim(-9,0)
	plt.vlines(Ea,0,k)
	plt.subplot(2,2,4)
	plt.title('eT Dynamics and Time Constant')
	plt.xlabel('time (ps)')
	plt.ylabel('Pop.')
	plt.plot(time,eTdynA,color='b')
	plt.plot(time,eTdynD,color='r')
	plt.legend(['A','D'],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.text(500,0.1,'tau = '+str(np.floor(1./np.sum(np.real(k))))+' ps')
	plt.show()

	np.savetxt("V_ab.dat",(Ea,V_ab))
	np.savetxt("eigenValsA.dat",Ea)
	np.savetxt("eigenValsB.dat",Eb)
	np.savetxt("marcus_eTrate.dat",(Ea,k))

