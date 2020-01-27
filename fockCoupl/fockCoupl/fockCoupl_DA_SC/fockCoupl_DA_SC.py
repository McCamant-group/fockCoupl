
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
	delta_d = np.loadtxt(sys.argv[4])
	S_d = delta_d**2/2.
	w_d = np.loadtxt(sys.argv[5])
	d_reorg = np.sum(S_d*w_d)/8100.
	delta_a = np.loadtxt(sys.argv[6])
	S_a = delta_a**2/2
	w_a = np.loadtxt(sys.argv[7])


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
		T = float(inp[4]) # K
		kBT = 0.695*T #
		Efb = float(inp[5]) # eV
		Eox = float(inp[6]) # eV O2=3.343 O1=3.207
		s_reorg = float(inp[7])
		a_reorg = float(inp[8])
		deltaG = (Efb-Eox)
		drivingForce = deltaG+(d_reorg+a_reorg+s_reorg)
		au2ev = 27.2113966413079
		ev2erg = 1.6e-12 # eV to erg 
		Econv_range = np.linspace(-2.,2.,2000)*8100.

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



	E = np.linspace(Efb-5,0,2000)
	V_cont = np.zeros(len(E),dtype=complex)
	theta=0.05
	#for i in range(len(Ea)):
	#	V_cont += ((V_tot[i]*8100.)**2)*np.exp(-((E-Ea[i]+d_tot_reorg+a_tot_reorg)**2)/(2*theta**2))

	j=0
	V_avg = 0
	for i in range(len(Ea)):
		if Ea[i]>=Efb and Ea[i]<0.0:
			V_avg += abs(V_tot[i])
			j+=1
	V_avg = (V_avg*8100./j)**2 
	print np.sqrt(V_avg)

	#plt.vlines(Ea,0,(abs(V_tot*8100)))
	#plt.plot(E,np.sqrt(V_cont))
	#plt.show()


	
	def DOS():
		E = np.linspace(Efb-5,0,2000)
		dos = np.zeros(len(E))
		hbar = 1.054e-27 # erg s
		ev2erg = 1.6e-12 # eV to erg 
		sigma = 136.084e-24 #32.5e-24 # cm^3
		me = 9*9.109e-28 # 0.00054858 amu # m0 = 9.109e-28 # g
		for i in range(len(E)):
			if (Efb) <= E[i]:
				dos[i] = (((2*me)**(3./2.))*sigma)/(2*np.pi**2*hbar**3)*np.sqrt(ev2erg*(E[i]-Efb))*ev2erg*(1./8100.)
			else:
				dos[i] = 0.0
		
		return dos	
	
	

	hbar=5.3088

	eta_d = 1./(np.exp(hbar*w_d/kBT)-1)
	eta_a = 1./(np.exp(hbar*w_a/kBT)-1)
	Eox = Eox*8100.	

	t = np.linspace(-0.25,0.25,2000)
	dt = (t[1]-t[0])
	E = np.linspace(Efb-5,0,2000)*8100.
	dE = (E[1]-E[0])	

	k = 0.1
	beta = 1./kBT
	L = k*np.sqrt(2*s_reorg*8100./beta)
	D = L/k
	gamma = D/((1+0.85*k+0.88*k**2)/(2.355+1.76*k))
	#D =  gamma*(1+0.85*k+0.88*k**2)/(2.355+1.76*k) # D parameter 
	#L =  k*D # LAMBDA parameter
	#s_reorg = beta*(L/k)**2/2
	print gamma
	#print 1./L	
	



	def g(t):	

		g = ((D/L)**2)*(L*t/hbar-1+np.exp(-L*t/hbar))+1j*((beta*D**2)/(2*L))*(1-np.exp(-L*t/hbar)) 
		#g = p.gamma*np.abs(t)#
		#s_reorg = beta*(L/k)**2/2
		#print 1./L
		return np.exp(-g)	

	def tdwp_d(t):
		tdwp_d = np.zeros((len(w_d),2000),dtype=complex)
		for i in range(len(S_d)):
			#tdwp[i,:] = S[i]*((2*eta[i]+1)*(1-np.cos(w[i]*t/hbar))+1j*np.sin(w[i]*t/hbar))
			tdwp_d[i,:] = 0.5*(1+eta_d[i])*S_d[i]*(1-np.exp(-1j*w_d[i]*t/hbar))+eta_d[i]*S_d[i]*(1-np.exp(1j*w_d[i]*t/hbar))
		return np.exp(-np.sum(tdwp_d,axis=0))	

	def tdwp_a(t):
		tdwp_a = np.zeros((len(w_a),2000),dtype=complex)
		for i in range(len(S_a)):
			#tdwp[i,:] = S[i]*((2*eta[i]+1)*(1-np.cos(w[i]*t/hbar))+1j*np.sin(w[i]*t/hbar))
			tdwp_a[i,:] = 0.5*(1+eta_a[i])*S_a[i]*(1-np.exp(-1j*w_a[i]*t/hbar))+eta_a[i]*S_a[i]*(1-np.exp(1j*w_a[i]*t/hbar))
		return np.exp(-np.sum(tdwp_a,axis=0))


	def FC():	

		integrand = np.zeros((len(t),len(E)),dtype=complex)
		integral = np.zeros(len(E),dtype=complex)
		
		
		tt,EE = np.meshgrid(t,E,sparse=True)	

		integrand = np.exp(-1j*(EE-Eox)*tt/hbar)*g(tt)*tdwp_a(tt)*tdwp_d(tt)
		

		integral = 0.5*(1/hbar)*np.real(np.trapz(integrand))*dt
		


		return integral
	

	FC = FC()
	kE = (2./hbar)*np.real(V_avg)*FC

	#plt.plot(FC)
	#plt.show()
	#exit()
	DOS = DOS()
	
	rate = (2./hbar)*np.trapz(np.real(V_avg)*FC*DOS)*dE
	print rate
	print 1./rate
	
	#ovlp = ovlp()	

	plt.plot(E/8100.,kE)
	plt.plot(E/8100.,DOS)	

	plt.show()	
	
	if any([i == 'data' for i in os.listdir('./')]) == True:
		pass
	else:
		os.mkdir('./data')

	np.savetxt("data/FC.dat",FC)
	np.savetxt("data/kE.dat",kE)
	np.savetxt("data/DOS.dat",DOS)
	np.savetxt("data/reorg.dat",S_a*w_a)
	np.savetxt("data/E_axis.dat",E/8100.)
	np.savetxt("data/V_cont.dat",np.real(np.sqrt(V_cont)))

	with open("data/output.txt",'w') as o:
		#o.write("V = "),o.write(str(V_avg)),o.write(" cm-1 \n")
		o.write("gamma = "),o.write(str(gamma)),o.write(" cm-1 \n")
		o.write("V = "),o.write(str(np.sqrt(V_avg))),o.write(" cm-1 \n")
		o.write("d_tot_reorg = "),o.write(str(d_reorg)),o.write(" eV \n")
		o.write("a_tot_reorg = "),o.write(str(a_reorg)),o.write(" eV \n")
		o.write("rate = "),o.write(str(rate)),o.write(" ps-1 \n")
		o.write("tau = "),o.write(str(1./rate)),o.write(" ps \n")
		o.write("deltaG = "),o.write(str(deltaG)), o.write(" eV \n")
		o.write("deltaG + tot_reorg = "),o.write(str(drivingForce)), o.write(" eV \n")

	o.close()

