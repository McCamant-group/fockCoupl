import numpy as np
from numpy.linalg import matrix_power
from numpy.linalg import eig
from scipy.linalg import sqrtm
import matplotlib.pyplot as plt
import warnings 

warnings.filterwarnings('ignore') # Supresses 'casting to real discards complex part' warning

###################### Calculate a-b electronic coupling ######################

dim = 1399 # total number of basis functions
dima = 988 # basis functions of fraction a
dimb = 411 # basis functions of fraction b

#aorb = 347 # Orbital of interest fraction a
#borb = 161 # Orbital of interest fraction b 

#aorb = aorb-1 # because Python counts from zero? shouldn't need this
#borb = borb-1

f = np.loadtxt("fock.dat") # read Fock matrix
s = np.loadtxt("overlap.dat") # read overlap matrix

## S^-1/2 ##
s_sqrt = sqrtm(s)
s_inv_sqrt = matrix_power(s_sqrt,-1)

## F_orth = S^-1/2 * F * S^-1/2
f_orth = np.dot(s_inv_sqrt,np.dot(f,s_inv_sqrt))

## Diagonalize aa, bb and cc blocks of F_orth ##
f_aa_eval,f_aa_evec = eig(f_orth[:dima,:dima])
f_bb_eval,f_bb_evec = eig(f_orth[dima:(dima+dimb),dima:(dima+dimb)])
#f_cc_eval,f_cc_evec = eig(f_orth[(dima+dimb):dim,(dima+dimb):dim])

## Use eigenvectors of above block diagonalization ##
## to calculate off diagonal blocks which are the coupling ##
f_ab = np.dot(f_aa_evec.conj().T,np.dot(f_orth[:dima,dima:(dima+dimb)],f_bb_evec))
#f_bc = np.dot(f_bb_evec.conj().T,np.dot(f_orth[dima:(dima+dimb),(dima+dimb):dim],f_cc_evec))
#f_ac = np.dot(f_aa_evec.conj().T,np.dot(f_orth[:dima,(dima+dimb):dim],f_cc_evec))

Ea = np.real(f_aa_eval)*27.2113966413079
Eb = np.real(f_bb_eval)*27.2113966413079
#Ec = np.real(f_cc_eval)*27.2113966413079

#c_orb = 157 #150 157
#a_orb = 771 # CT populates acceptor orbital directly as initial condition
b_orb = 192 #217

V_ab = f_ab[:,b_orb]*27.2113966413079
#V_bc = f_bc[b_orb,c_orb]*27.2113966413079
#V_ac = f_ac[:,c_orb]*27.2113966413079
V_tot = V_ab #V_bc*V_ab/((Eb[b_orb]-Ec[c_orb])*27.2113966413079)


hbar = 4.135667662e-3 # eV*ps
T = 298.0 # K
kBT = 8.6173303e-5*T # eV
reorg = 700./8100. #8100cm-1/eV
#Ec_orb = Ec[c_orb]
Eb_orb = Eb[b_orb] #sorted(Eb)[borb]
dt = 0.1
time = np.arange(0,50,dt)
A = np.zeros(len(time)+1)
A[0] = 1
B = np.zeros((len(time)+1,len(Ea)))
#B[0,a_orb] = 0.14 
kf = np.zeros(len(Ea))
kb = np.zeros(len(Ea))

#############################################################################
######################### Calculate the Marcus eT rate ######################
########### for the donor orbital coupled to all acceptor orbitals ##########
#############################################################################

## individual rates ##

for i in range(len(V_tot)):
	## Forward ##
	kf[i] = ((V_tot[i]*np.conj(V_tot[i]))/hbar)*np.sqrt(np.pi/(reorg*kBT))*np.exp(-((reorg+(Ea[i]-Eb_orb))**2/(4*reorg*kBT)))

	## Backward ##
	kb[i] = ((V_tot[i]*np.conj(V_tot[i]))/hbar)*np.sqrt(np.pi/(reorg*kBT))*np.exp(-((reorg+(Eb_orb-Ea[i]))**2/(4*reorg*kBT)))


# print((Ea-Ec_orb)),print(kf),print(kb)

# plt.subplot(2,1,1)
# plt.xlim(-4,0)
# # plt.vlines(Ea,0,(Ea-Ec_orb),colors='r')
# plt.vlines(Ea,0,kf,colors='g')
# plt.subplot(2,1,2)
# plt.xlim(-4,0)
# plt.vlines(Ea,0,kb,colors='b')
# plt.show()



#k_tot = np.sum(k)
#tau = 1./k_tot

for i in range(len(time)):
	A[i+1] = -(dt)*(A[i]*np.sum(kf)-np.dot(kb,B[i,:]))+A[i]
	B[i+1,:] = -(dt)*(np.dot(kb,B[i,:])-np.sum(kf)*A[i])+B[i,:]



## Dynamics with total rate ##
#eTdynD = np.exp(-k_tot*time)
#eTdynA = 1.0 - eTdynD

#np.savetxt("V_ab.dat",(Ea,V_ab))
np.savetxt("eigenValsA_srt.dat",sorted(Ea))
np.savetxt("eigenValsB_srt.dat",sorted(Eb))
#np.savetxt("eigenValsC_srt.dat",sorted(Ec))
np.savetxt("eigenValsA.dat",Ea)
np.savetxt("eigenValsB.dat",Eb)
#np.savetxt("eigenValsC.dat",Ec)
np.savetxt("V_ab.dat",np.real(V_tot))
np.savetxt("kf.dat",kf)
np.savetxt("kb.dat",kb)
np.savetxt("dynamicsA.dat",A)
np.savetxt("dynamicsB.dat",B)
np.savetxt("time.dat",time)

plt.subplot(2,2,1)
plt.title('A, B, Block Eigenvalues')
plt.ylim(0,3)
plt.xlim(-4,-2)
plt.vlines(Ea,0,1,colors='b')
plt.vlines(Eb,1,2,colors='r')
#plt.vlines(Ec,2,3,colors='g')
plt.legend(['A','B'])#,bbox_to_anchor=(1.05, 1), loc=9, borderaxespad=0.)
plt.subplot(2,2,2)
plt.title('A, B Coupling')
plt.ylabel('V_AB (meV)')
plt.ylim(0,20)
plt.xlim(-4,-2)
#plt.vlines(Ec_orb,0,abs(V_bc)*1e3,colors='b')
plt.vlines(Eb_orb,10,100) #plt.vlines(Eb[b_orb],0,abs(V_bc)*1e3,colors='b')
plt.vlines(Ea,0,abs(V_tot)*1e3,colors='r')
#plt.vlines(Ea,0,abs(V_tot)*1e3,colors='g')
plt.subplot(2,2,3)
plt.xlim(-4,-2)
plt.vlines(Ea,0,kf)
plt.vlines(Ea,0,-kb,colors='r')
plt.ylabel("rate (1/ps)")
plt.xlabel('Orbital Energy (eV)')
plt.subplot(2,2,4)
plt.plot(np.append(time,time[len(time)-1]+dt),A)
plt.ylabel("B Pop.")
plt.xlabel("Time (ps)")
plt.show()

'''
## Print some relavant output data ##
with open("output.txt",'w') as o:
	o.write("Results: \n")
	o.write("Donor orbital energy = "), o.write(str(Eb_orb)),o.write(" eV \n")
	#o.write("Rate Constant = "), o.write(str(k_tot)), o.write(" ps^-1 \n")
	#o.write("Time Constant = "), o.write(str(tau)), o.write(" ps \n")
o.close()


np.savetxt("V_ab.dat",(Ea,V_ab))
np.savetxt("eigenValsA.dat",Ea)
np.savetxt("eigenValsB.dat",Eb)
np.savetxt("eigenValsC.dat",Ec)
np.savetxt("eigenValsC_srt.dat",sorted(Eb))
np.savetxt("marcus_eTrate.dat",(Ea,k))

'''
