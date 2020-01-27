
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
	#w = np.loadtxt(sys.argv[4])
	#S = np.loadtxt(sys.argv[5])


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
		hbar = 4.135667662e-3  # plancks constant 5.3088 cm^-1*ps #4.135667662e-3 # eV*ps
		T = float(inp[4]) # K
		#kBT = 8.6173303e-5*T # eV
		kBT = 0.695*T # cm-1
		reorg = float(inp[5])/8065.54 #8100cm-1/eV
		reorg_range = np.arange(200.0,1100.0,50.0)/8100. 
		dt = float(inp[6])
		time = np.arange(0.0,float(inp[7]),dt) #time = np.logspace(1,7,10000,base=2.71)
		time_adj = np.append(time,[time[(len(time)-1)]+time[len(time)-1]-time[len(time)-2]]) #np.append(time,[time[(len(time)-1)]+dt])
		#dt = np.zeros(len(time))

		theta = float(inp[8])*8100. # np.arange(0.0,2.0,0.2) #
		theta_range = np.array([0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.5,1.0,2.0])
		Econv_range = np.linspace(-1*float(inp[9]),float(inp[9]),float(inp[10]))*8100.
		au2ev = 27.2113966413079
		au2ev2cm1 = 8065.54*au2ev

	i.close()

	#kBT = kBT/8100.
	#V = 73/8100.
	#hbar = hbar/8100.
	#print((4*(np.pi**2)/hbar)*(V**2)*(1./np.sqrt(4*np.pi*kBT))*np.exp(-(((1.984-2.75)**2)/(4*1.984*kBT))))
	

	
	#for i in range(len(time)-1):
		#dt[i] = time[i+1]-time[i]

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
			
	V_avg = 1.6e-12*V_avg/j
	print V_avg

	Efb = 3.86 # eV
	Eox = 3.207 # eV O2=3.343 O1=3.207

	reorg = 0.56 # eV

	#kf = np.zeros(len(Ea))
	#kb = np.zeros(len(Ea))

	
	#############################################################################
	######################### Calculate the Marcus eT rate ######################
	########### for the donor orbital coupled to all acceptor orbitals ##########
	#############################################################################



	def rate():
		hbar = 1.054e-27 # erg s
		sigma = 32.5e-24 # cm^3
		me = 15*9.109e-28 # g
		C = (((2*me)**(3./2.))*(V_avg**2)*sigma)/(np.pi*hbar**4)
		ket = C*np.sqrt(-1.6e-12*(Eox-Efb-reorg)) #
		print C
		print ket
		print (1/ket)*1e12

	rate()


	#def H(a,theta):
	#	if theta == 0.0:
	#		H = np.zeros(len(Econv_range))
	#		H[(len(Econv_range)/2)] = 1
	#	else:
	#		H = (1/(theta*np.sqrt(2*np.pi)))*np.exp(-((Econv_range-a)**2)/(2*theta**2))
	#	return H


	#Ea_active = []
	#V_active = []
	#for k in range(len(Ea)):
	#	if Ea[k] < -24000. and Ea[k] > -26000.: 
	#		Ea_active.append(Ea[k])
	#		V_active.append(V_tot[k])
	#Ea_active = np.array(Ea_active)
	#V_active = np.array(np.abs(V_active))

	#print V_active
	#print Eb_orb

	#eta = 1./(np.exp(hbar*w/kBT)-1)
	#t = np.linspace(0.0,0.25,1000)


	#def g(t):
	#	k = 0.1
	#	gamma = 5200.
	#	beta = 1./kBT
	#	D =  gamma*(1+0.85*k+0.88*k**2)/(2.355+1.76*k) # D parameter 
	#	L =  k*D # LAMBDA parameter
	#	g = ((D/L)**2)*(L*t/hbar-1+np.exp(-L*t/hbar))+1j*((beta*D**2)/(2*L))*(1-np.exp(-L*t/hbar)) 
	#	#g = p.gamma*np.abs(t)#
	#	s_reorg = beta*(L/k)**2/2
	#	#print 1./L
	#	return np.exp(-g)
#
	#def L(t):
	#	L = np.zeros((len(w)))
	#	for i in range(len(S)):
	#		L[i] = S[i]*((2*eta[i]+1)*(1-np.cos(w[i]*t/hbar))+1j*np.sin(w[i]*t/hbar))
	#	return np.exp(-np.sum(L,axis=0))
#
	#def rate():
	#	integ_f = np.zeros((len(Ea_active),len(t)))
	#	integ_b = np.zeros((len(Ea_active),len(t)))
	#	rate_f_conv = np.zeros(len(Ea_active))
	#	rate_b_conv = np.zeros(len(Ea_active))
#
	#	for i in range(len(Ea_active)):
	#		for j in range(len(t)):
	#			#integ_f[i,j] = np.sum(np.convolve(2*((V_active[i]/hbar)**2)*np.exp(-1j*(Ea_active[i]-Eb_orb+Econv_range)*t[j]/hbar)*L(t[j])*g(t[j]),H(0,theta),'valid')/np.sum(H(0,theta)))
	#			#integ_b[i,j] = np.sum(np.convolve(2*((V_active[i]/hbar)**2)*np.exp(-1j*(Eb_orb-Ea_active[i]+Econv_range)*t[j]/hbar)*L(t[j])*g(t[j]),H(0,theta),'valid')/np.sum(H(0,theta)))
	#			integ_f[i,j] = np.sum(np.convolve(2*((V_active[i]/hbar)**2)*np.exp(-1j*(G0+Econv_range)*t[j]/hbar)*L(t[j])*g(t[j]),H(0,theta),'valid')/np.sum(H(0,theta)))
	#			integ_b[i,j] = np.sum(np.convolve(2*((V_active[i]/hbar)**2)*np.exp(-1j*(-G0+Econv_range)*t[j]/hbar)*L(t[j])*g(t[j]),H(0,theta),'valid')/np.sum(H(0,theta)))
	#	plt.plot(t,np.transpose(np.real(integ_f)))
	#	plt.show()
	#	rate_f = np.real(np.trapz(integ_f,axis=1))*(t[1]-t[0])
	#	rate_b = np.real(np.trapz(integ_b,axis=1))*(t[1]-t[0])
#
	#	return(rate_f,rate_b)
#
#
	#rates = rate()
	#kf = np.abs(rates[0])
	#kb = np.abs(rates[1])
	#print kf
	#print kb
	##print np.shape(rates)
	##print rates[0]-rates[1]
	#plt.vlines(Ea_active,0,kf)
	#plt.vlines(Ea_active,0,-kb,colors='r')
	#plt.show()
#
#
#
	#	## Population dynamics ##
#
	#D = np.zeros(len(time)+1)
	#D[0] = 1
	##B_deact = np.zeros(len(time)+1)
	##A_ss = 0
	#A = np.zeros((len(time)+1,len(Ea_active)))
	##A_bulk = np.zeros(len(time)+1)
	##A_ss = np.zeros(len(Ea))
#
	#for i in range(len(time)):
	#	D[i+1] = -(dt)*(D[i]*np.sum(kf)-np.dot(kb,A[i,:]))+D[i]
	#	A[i+1,:] = -(dt)*(kb[:]*A[i,:]-kf*D[i])+A[i,:]
#
	#plt.plot(time_adj,A)
	#plt.plot(time_adj,D)
	#plt.show()


####
	#def dynamics(theta,reorg):
	#def dynamics():


	## individual rates ##

		#en_f = np.zeros((len(Ea)))
		#en_b = np.zeros((len(Ea)))

#		for i in range(len(V_tot)):
#			## Forward ##
#			# Convolution with distribution
#			en_f[i] = np.sum(np.convolve(np.exp(-((reorg+(Ea[i]-Eb_orb)+Econv_range)**2/(4*reorg*kBT))),H(0,theta),'valid')/np.sum(H(0,theta)))
#			kf[i] = ((V_tot[i]*np.conj(V_tot[i]))/hbar)*np.sqrt(np.pi/(reorg*kBT))*en_f[i]#

#			## Backward ##
#			en_b[i] = np.sum(np.convolve(np.exp(-((reorg+(Eb_orb-Ea[i]+Econv_range))**2/(4*reorg*kBT))),H(0,theta),'valid')/np.sum(H(0,theta)))
#			kb[i] = ((V_tot[i]*np.conj(V_tot[i]))/hbar)*np.sqrt(np.pi/(reorg*kBT))*en_b[i]#

#		Ea_conv = np.zeros(len(Econv_range))
#		for i in range(len(Ea)):
#			Ea_conv += np.convolve(1,H(Ea[i],theta),'valid')

		## Population dynamics ##

		# Purely Ab Initio dynamics #
		#for i in range(len(time)):
		#	B[i+1] = -(dt[i])*(B[i]*np.sum(kf)-np.dot(kb,A[i,:]))+B[i]
		#	A[i+1,:] = -(dt[i])*(kb[:]*A[i,:]-kf*B[i])+A[i,:]

		# FrankenDynamics #
		#k_deact = 1./52.7
		#kBER = 1./(.160)
		#for i in range(len(time)):
		#	B[i+1] = -(dt[i])*(B[i]*np.sum(kf)-np.dot(kb,A[i,:])+k_deact*B[i])+B[i]
		#	B_deact[i+1] = (dt[i])*(k_deact*B[i])+B_deact[i]
		#	for k in range(len(Ea)):
		#		if Ea[k] < -2.2 and Ea[k] > -2.8:
		#			#if Ea[k]-Ea[k-1] < 1.0:
		#			A[i+1,k] = -(dt[i])*(kb[k]*A[i,k]-kf[k]*B[i]+kBER*A[i,k])+A[i,k]
		#	A_bulk[i+1] = (dt[i])*(kBER*np.sum(A[i,:]))+A_bulk[i]
					#else:
						#A[i+1,k] = -(dt[i])*(kb[k]*A[i,k]-kf[k]*B[i])+A[i,k]

		##	

		#for i in range(len(kb)):
		#	if kb[i] == 0.0:
		#		kb[i] =  np.min(kb[np.nonzero(kb)])

		##

		#B_ss = 1./(1+np.sum(kf/kb))
		#A_ss = B_ss*(kf/kb)


		#return(Ea_conv, A, B, kf, kb, A_ss, B_ss)

	#def fit_func(t,a1,tau1,a2,tau2):
	#	return a1*np.exp(-t/tau1)+a2*np.exp(-t/tau2)



	

	#(Ea_conv, A, B, kf, kb, A_ss, B_ss) = dynamics(0.05,700./8100.)

	#(popt,pcov) = curve_fit(fit_func,time_adj,A,p0=(0.7,0.5,0.5,50.0))


	### Write the main data ##
#	path = "frankenDynamics/"
#	try:
#		os.mkdir(path)
#	except OSError:
#		print("Creation of directory %s failed" % path)
#		exit()
#	else:
#		print("Directory %s created, data being saved" % path)
#	with open(path+"B_ss.dat",'w') as f:
#		f.write("Donor orbital energy = "+str(Eb_orb)+"\n")
#		f.write("Donor steady state population = "+str(B_ss))
#	f.close()
#	np.savetxt(path+"A_ss.dat",np.c_[Ea,A_ss])
#	np.savetxt(path+"V_kf_kb.dat",np.c_[Ea,np.real(V_tot),kf,kb])
#	np.savetxt(path+"eigenValsA_srt.dat",np.c_[Ea,sorted(Ea)])
#	np.savetxt(path+"eigenValsB_srt.dat",np.c_[Eb,sorted(Eb)])
#	np.savetxt(path+"dynamicsA.dat",np.c_[time_adj,B,A])
#	np.savetxt(path+"eigenValsA_conv.dat",np.c_[Econv_range,Ea_conv/np.sum(H(0,theta))])#

#	## save a summary plot ##
#	fig = plt.figure()
#	plt.subplot(2,2,1)
#	plt.title('A, B, Block Eigenvalues')
#	plt.ylim(0,2)
#	plt.xlim(-10,0)
#	plt.vlines(Ea,0,0.5,colors='b')
#	plt.plot(Econv_range,Ea_conv/np.sum(H(0,theta)))
#	plt.vlines(Eb,0.5,1,colors='r')
#	plt.legend(['B','A'])
#	plt.subplot(2,2,2)
#	plt.title('A, B Coupling')
#	plt.ylabel('V_AB (meV)')
#	plt.ylim(0,20)
#	plt.xlim(-10,0)
#	plt.vlines(Eb_orb,10,100) 
#	plt.vlines(Ea,0,abs(V_tot)*1e3,colors='r')
#	plt.subplot(2,2,3)
#	plt.xlim(-10,0)
#	plt.vlines(Ea,0,kf)
#	plt.vlines(Ea,0,-kb,colors='r')
#	plt.ylabel("rate (1/ps)")
#	plt.xlabel('Orbital Energy (eV)')
#	plt.subplot(2,2,4)
#	plt.plot(time_adj,B)
#	plt.plot(time_adj,B_deact)
#	plt.plot(time_adj,A_bulk)
#	plt.ylabel("B Pop.")
#	plt.xlabel("Time (ps)")
#	plt.savefig(path+"summary")#

#	fig = plt.figure()
#	for k in range(len(Ea)):
#		if Ea[k] < -2.3 and Ea[k] > -3.2:
#			plt.plot(time_adj,A[:,k])
#	plt.show()
#	plt.savefig(path+"acceptor_dyn")#

#	exit()
########
#	
#	path = "theta/"
#	try:
#		os.mkdir(path)
#	except OSError:
#		print("Creation of directory %s failed" % path)
#		exit()
#	else:
#		print("Directory %s created, data being saved" % path)#

#	B_equil_theta = np.zeros(len(theta_range))
#	B_dyn_theta = np.zeros((len(theta_range),len(time_adj)))
#	A_dyn_theta = np.zeros((len(theta_range),len(time_adj),len(Ea)))
#	kf_theta = np.zeros((len(theta_range),len(V_tot)))
#	kb_theta = np.zeros((len(theta_range),len(V_tot)))
#	A_ss_theta = np.zeros((len(theta_range),len(Ea)))
#	B_ss_theta = np.zeros((len(theta_range)))#

#	l=0#

#	for i in theta_range:
#		
#		print(i,reorg)
#		(Ea_conv, A_dyn_theta[l], B_dyn_theta[l], kf_theta[l], kb_theta[l],A_ss_theta[l], B_ss_theta[l]) = dynamics(i,reorg)
#		B_equil_theta[l] = B_dyn_theta[l][len(time_adj)-1]#

#		path = "theta/theta_"+str(l)+"/"
#		try:
#			os.mkdir(path)
#		except OSError:
#			print("Creation of directory %s failed" % path)
#			exit()
#		else:
#			print("Directory %s created, data being saved" % path)
#		with open(path+"B_ss.dat",'w') as f:
#			f.write("Donor orbital energy = "+str(Eb_orb)+"\n")
#			f.write("Donor steady state population = "+str(B_ss))
#		f.close()
#		np.savetxt(path+"A_ss.dat",np.c_[Ea,A_ss])
#		np.savetxt(path+"kf_kb.dat",np.c_[Ea,kf_theta[l],kb_theta[l]])
#		np.savetxt(path+"dynamics.dat",np.c_[time_adj,B_dyn_theta[l],A_dyn_theta[l]])
#		np.savetxt(path+"eigenValsA_conv.dat",np.c_[Econv_range,Ea_conv/np.sum(H(0,i))])
#		l+=1#

#	path = "reorg/"
#	try:
#		os.mkdir(path)
#	except OSError:
#		print("Creation of directory %s failed" % path)
#		exit()
#	else:
#		print("Directory %s created, data being saved" % path)#

#	B_equil_reorg = np.zeros(len(reorg_range))
#	B_dyn_reorg = np.zeros((len(reorg_range),len(time_adj)))
#	A_dyn_reorg = np.zeros((len(reorg_range),len(time_adj),len(Ea)))
#	kf_reorg = np.zeros((len(reorg_range),len(Ea)))
#	kb_reorg = np.zeros((len(reorg_range),len(Ea)))
#	A_ss_reorg = np.zeros((len(reorg_range),len(Ea)))
#	B_ss_reorg = np.zeros((len(reorg_range)))
#	l=0#

#	for j in reorg_range:
#		
#		print(theta,j)
#		(Ea_conv, A_dyn_reorg[l], B_dyn_reorg[l], kf_reorg[l], kb_reorg[l],A_ss_reorg[l],B_ss_reorg[l]) = dynamics(theta,j)
#		B_equil_reorg[l] = B_dyn_reorg[l][len(time_adj)-1]#

#		path = "reorg/reorg_"+str(l)+"/"
#		try:
#			os.mkdir(path)
#		except OSError:
#			print("Creation of directory %s failed" % path)
#			exit()
#		else:
#			print("Directory %s created, data being saved" % path)
#		with open(path+"B_ss.dat",'w') as f:
#			f.write("Donor orbital energy = "+str(Eb_orb)+"\n")
#			f.write("Donor steady state population = "+str(B_ss))
#		f.close()
#		np.savetxt(path+"A_ss.dat",np.c_[Ea,A_ss])
#		np.savetxt(path+"kf_kb.dat",np.c_[Ea,kf_reorg[l],kb_reorg[l]])
#		np.savetxt(path+"dynamics.dat",np.c_[time_adj,B_dyn_reorg[l],A_dyn_reorg[l]])
#		np.savetxt(path+"eigenValsA_conv.dat",np.c_[Econv_range,Ea_conv/np.sum(H(0,theta))])
#		l+=1



