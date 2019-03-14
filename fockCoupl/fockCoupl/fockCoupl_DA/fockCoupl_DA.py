
def main():
	import numpy as np
	from numpy.linalg import matrix_power
	from numpy.linalg import eig
	from scipy.linalg import sqrtm
	import matplotlib.pyplot as plt
	import sys
	import warnings 

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
		hbar = 4.135667662e-3 # eV*ps
		T = float(inp[4]) # K
		kBT = 8.6173303e-5*T # eV
		reorg = float(inp[5])/8100 #8100cm-1/eV
		dt = float(inp[6])
		time = np.arange(0,float(inp[7]),dt)
		theta = float(inp[8])
		Econv_range = np.linspace(-1*float(inp[9]),float(inp[9]),float(inp[10]))
		au2ev = 27.2113966413079

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

	Ea = np.real(f_aa_eval)*au2ev
	Eb = np.real(f_bb_eval)*au2ev
	Eb_orb = Eb[b_orb]

	V_ab = f_ab[:,b_orb]*au2ev
	V_tot = V_ab 

	A = np.zeros(len(time)+1)
	A[0] = 1
	B = np.zeros((len(time)+1,len(Ea)))

	kf = np.zeros(len(Ea))
	kb = np.zeros(len(Ea))

	#############################################################################
	######################### Calculate the Marcus eT rate ######################
	########### for the donor orbital coupled to all acceptor orbitals ##########
	#############################################################################

	def H(a):
		if theta == 0.0:
			H = np.zeros(len(Econv_range))
			H[(len(Econv_range)/2)] = 1
		else:
			H = (1/(theta*np.sqrt(2*np.pi)))*np.exp(-((Econv_range-a)**2)/(2*theta**2))
		return H

	## individual rates ##

	en_f = np.zeros((len(Ea)))
	en_b = np.zeros((len(Ea)))

	for i in range(len(V_tot)):
		## Forward ##
		# Convolution with distribution
		en_f[i] = np.sum(np.convolve(np.exp(-((reorg+(Ea[i]-Eb_orb)+Econv_range)**2/(4*reorg*kBT))),H(0),'valid')/np.sum(H(0)))
		kf[i] = ((V_tot[i]*np.conj(V_tot[i]))/hbar)*np.sqrt(np.pi/(reorg*kBT))*en_f[i]

		## Backward ##
		en_b[i] = np.sum(np.convolve(np.exp(-((reorg+(Eb_orb-Ea[i]+Econv_range))**2/(4*reorg*kBT))),H(0),'valid')/np.sum(H(0)))
		kb[i] = ((V_tot[i]*np.conj(V_tot[i]))/hbar)*np.sqrt(np.pi/(reorg*kBT))*en_b[i]

	Ea_conv = np.zeros(len(Econv_range))
	for i in range(len(Ea)):
		Ea_conv += np.convolve(1,H(Ea[i]),'valid')

	## Population dynamics ##

	for i in range(len(time)):
		A[i+1] = -(dt)*(A[i]*np.sum(kf)-np.dot(kb,B[i,:]))+A[i]
		B[i+1,:] = -(dt)*(np.dot(kb,B[i,:])-np.sum(kf)*A[i])+B[i,:]

	time = np.append(time,[time[(len(time)-1)]+dt])

	## Write the data ##


	np.savetxt("V_ab.dat",np.c_[Ea,np.real(V_tot),kf,kb])
	np.savetxt("eigenValsA_srt.dat",np.c_[Ea,sorted(Ea)])
	np.savetxt("eigenValsB_srt.dat",np.c_[Eb,sorted(Eb)])
	np.savetxt("dynamicsA.dat",np.c_[time,A,1-A])

	## save a summary plot ##

	plt.subplot(2,2,1)
	plt.title('A, B, Block Eigenvalues')
	plt.ylim(0,2)
	plt.xlim(-4,-2)
	plt.vlines(Ea,0,0.5,colors='b')
	plt.plot(Econv_range,Ea_conv/np.sum(H(0)))
	plt.vlines(Eb,0.5,1,colors='r')
	plt.legend(['A','B'])
	plt.subplot(2,2,2)
	plt.title('A, B Coupling')
	plt.ylabel('V_AB (meV)')
	plt.ylim(0,20)
	plt.xlim(-4,-2)
	plt.vlines(Eb_orb,10,100) 
	plt.vlines(Ea,0,abs(V_tot)*1e3,colors='r')
	plt.subplot(2,2,3)
	plt.xlim(-4,-2)
	plt.vlines(Ea,0,kf)
	plt.vlines(Ea,0,-kb,colors='r')
	plt.ylabel("rate (1/ps)")
	plt.xlabel('Orbital Energy (eV)')
	plt.subplot(2,2,4)
	plt.plot(time,A)
	plt.ylabel("B Pop.")
	plt.xlabel("Time (ps)")
	plt.savefig("summary")

