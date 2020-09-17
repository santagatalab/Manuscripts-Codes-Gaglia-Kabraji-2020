#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
from math import *
from params import *


# functions
def BB(A1,A2,A3,A4):
	return A2-A1+A3*A2+A4*A1

def GK(A1,A2,A3,A4):
	return 2.0*A4*A1/(BB(A1,A2,A3,A4)+ \
		(BB(A1,A2,A3,A4)**2.0-4.0*(A2-A1)*A4*A1)**.5)

def Mass_u5f_(k1):
	return k1

def Mass_0001(k1,S1):
	return k1*S1

def Mass_0002(k1,S1,S2):
	return k1*S1*S2

def Michaelis(M1,J1,k1,S1):
	return k1*S1*M1/(J1+S1)

# Species

def PreMPF(x):
	return x[10]+x[11]

def TriB(x):
	return x[4]+x[11]

def CycBT(x):
	return x[1]+x[10]+x[4]+x[11]

def CycAT(x):
	return x[0]+x[12]

def CycET(x):
	return x[2]+x[13]

def CycD(x):
	return CycD0*x[9]*palbo_mult

def CKIT(x):
	return x[8]+x[4]+x[11]+x[12]+x[13]

def Cdc20T(x):
	return x[5]+x[6]

def Cdc14(x):
	return x[5]

def Wee1(x):
	return GK(kawee_p+kawee_pp*Cdc14(x),kiwee_p+kiwee_pp*x[1],Jawee,Jiwee)

def Vwee(x):
	return kwee_p+kwee_pp*Wee1(x)

def Cdc25(x):
	return GK(ka25_p+ka25_pp*x[1],ki25_p+ki25_pp*Cdc14(x),Ja25,Ji25)

def V25(x):
	return k25_p+k25_pp*Cdc25(x)

def TFB(x):
	return GK(kafb*x[1],kifb,Jafb,Jifb)

def Vatf(x):
	return katf_p+katfa_pp*x[0]+katfe_pp*x[2]+katfd_pp*CycD(x)

def Vitf(x):
	return kitf_p+kitfa_pp*x[0]+kitfb_pp*x[1]

def TFE(x):
	return GK(Vatf(x),Vitf(x),Jatf,Jitf)

def TFI(x):
	return GK(kafi*Cdc14(x),kifi_p+kifib_pp*x[1],Jafi,Jifi)

def Vsb(x):
	return (ksb_p+ksb_pp*TFB(x))*x[9]

def Vsa(x):
	return (ksa_p+ksa_pp*TFE(x))*x[9]

def Vse(x):
	return (kse_p+kse_pp*TFE(x))*x[9]

def Vah1(x):
	return kah1_p+kah1_pp*Cdc14(x)

def Vih1(x):
	return kih1_p+kih1a_pp*x[0]+kih1b_pp*x[1]+kih1e_pp*x[2]+kih1d_pp*CycD(x)

def Vdb(x):
	return kdb_p+kdbh_pp*x[7]+kdbc_pp*x[5]

def Vda(x):
	return kda_p+(kda_pp+kda_ppp)*x[5]+kda_ppp*x[6]

def Vde(x):
	return kde_p+kdee_pp*x[2]+kdea_pp*x[0]+kdeb_pp*x[1]

def Vsi(x):
	return ksi_p+ksi_pp*TFI(x)

def Vdi(x):
	return (kdi_p+kdia_pp*x[0]+kdib_pp*x[1]+ \
		kdie_pp*x[2]+kdid_pp*CycD(x))/(1.0+k14di*Cdc14(x))

#Velocities and Dependent Variables

def APC(x):
	return (-(x[3]-APCT)/1)

def Cdh1i(x):
	return (-(x[7]-Cdh1T)/1)
 
# Independent species

def derivCycA(x):
	return - (kassa*x[8]*x[0]) + (kdissa*x[12]) + (Vdi(x)*x[12]) + \
	(Vsa(x)) - (Vda(x)*x[0])

def derivCycB(x):
	return (Vsb(x)) - (Vdb(x)*x[1]) + (V25(x)*x[10]) - (Vwee(x)*x[1]) - \
	(kassb*x[1]*x[8]) + (kdissb*x[4]) + (Vdi(x)*x[4])

def derivCycE(x):
	return - (kasse*x[8]*x[2]) + (kdisse*x[13]) + (Vdi(x)*x[13]) + \
	(Vse(x)) - (Vde(x)*x[2])

def derivAPCP(x):
	return (Michaelis(x[1],Jaie,kaie,APC(x))) - \
	(Michaelis(1.0,Jiie,kiie,x[3]))

def derivBCKI(x):
	return (kassb*x[1]*x[8]) - (kdissb*x[4]) + (V25(x)*x[11]) - \
	(Vwee(x)*x[4]) - (Vdb(x)*x[4]) - (Vdi(x)*x[4])

def derivCdc20A(x):
	return (Michaelis(x[3],Ja20,ka20,x[6])) - \
	 (Michaelis(1.0,Ji20,ki20,x[5])) - (kd20*x[5])

def derivCdc20i(x):
	return (ks20_p+ks20_pp*x[1]**n20/(J20**n20+x[1]**n20)) - (kd20*x[6]) - \
	 (Michaelis(x[3],Ja20,ka20,x[6])) + (Michaelis(1.0,Ji20,ki20,x[5]))

def derivCdh1(x):
	return (Michaelis(Vah1(x),Jah1,1.0,Cdh1i(x))) - \
	(Michaelis(Vih1(x),Jih1,1.0,x[7]))

def derivCKI(x):
	return (kassb*x[1]*x[8]) + (kdissb*x[4]) - (kassb*x[10]*x[8]) + \
	 (kdissb*x[11]) + (Vdb(x)*x[4]) + (Vdb(x)*x[11]) + (Vsi(x)) - \
	 (Vdi(x)*x[8]) - (kassa*x[8]*x[0]) + (kdissa*x[12]) + (Vda(x)*x[12]) - \
	 (kasse*x[8]*x[2]) + (kdisse*x[13]) + (Vde(x)*x[13])

def derivMass(x):
	return mu*x[9]*(1-x[9]/MaxMass)

def derivPB(x):
	return  - (V25(x)*x[10]) + (Vwee(x)*x[1]) - (Vdb(x)*x[10]) - \
	(kassb*x[10]*x[8]) + (kdissb*x[11]) + (Vdi(x)*x[11])

def derivPBCKI(x):
	return (kassb*x[10]*x[8]) - (kdissb*x[11]) - (V25(x)*x[11]) + \
	 (Vwee(x)*x[4]) - (Vdb(x)*x[11]) - (Vdi(x)*x[11])

def derivTriA(x):
	return (kassa*x[8]*x[0]) - (kdissa*x[12]) - (Vdi(x)*x[12]) - (Vda(x)*x[12])

def derivTriE(x):
	return (kasse*x[8]*x[2]) - (kdisse*x[13]) - (Vdi(x)*x[13]) - (Vde(x)*x[13])



###################################################################################


def deriv_v2(x, t):
	global prev_sign, palbo_mult

	if x[1] < KEZ and prev_sign != -1:
		dMass = -x[9]*0.5/dt
	else:
		dMass = derivMass(x)
	prev_sign = -1 if x[1] < KEZ else 1

	if palbo_on and t > 600:  # after one cycle
		# print(palbo_mult)
		palbo_mult = 0 #palbo_mult * 0.1

	dCycA = derivCycA(x)  
	dCycB = derivCycB(x) 
	dCycE = derivCycE(x)
	dAPCP = derivAPCP(x)
	dBCKI = derivBCKI(x)
	dCdc20A = derivCdc20A(x)
	dCdc20i = derivCdc20i(x)
	dCdh1 = derivCdh1(x)
	dCKI = derivCKI(x) 
	# dMass = derivMass(x)
	dPB = derivPB(x)
	dPBCKI = derivPBCKI(x)
	dTriA = derivTriA(x)
	dTriE = derivTriE(x)

	return np.array([dCycA, dCycB, dCycE, dAPCP, dBCKI,  dCdc20A, dCdc20i, \
	        dCdh1, dCKI,  dMass, dPB,   dPBCKI, dTriA,   dTriE])

# MAIN 

# initial vector
xinit = [CycA, CycB, CycE, APCP, BCKI,  Cdc20A, Cdc20i, \
	     Cdh1, CKI,  Mass, pB,   pBCKI, TriA,   TriE];


# global
prev_sign = 0

# settings
noise_bkg = True
decouple = False
palbo_on = True
palbo_mult = 1
MaxMass = 1.75 # 1.75
Qe = 0.5 # quantum efficiency
t_integration = 100
B = 50  # background photon flux

save = True
label = 'b' + str(B) + '-mm' + str(MaxMass)
if not noise_bkg:
	label = 'no-noise-mm' + str(MaxMass)
if palbo_on:
	label = label + '-palbo'

# Choose timestep
dt=0.01
maxT = 300000

# Perform the integration
t=0
x_orig=[]
x_fluor = []
time = []
xi=xinit

for i in range(maxT):
	if decouple:
		# xi=(xi+deriv_v2(xi,t)*dt) * np.random.normal(1,0.005*sqrt(dt),(len(xi)))
		mult = np.random.uniform(-0.8, .8, len(xi))
		xi=xi+(deriv_v2(xi,t)*dt * 2**mult) #np.random.choice([2/3, 3/4, 4/5, 1, 5/4, 4/3, 3/2], len(xi)))
	else:
		xi=xi+deriv_v2(xi,t)*dt
		
	if noise_bkg:
		# signal = alpha_s * np.random.normal(xi, mult_s*abs(xi), (len(xi)))
		# bkg = alpha_b * np.random.normal(beta, mult_b*beta, (len(xi)))    # BACKGROUND	

		signal = xi*Qe*t_integration
		poisson_noise = np.random.poisson(signal + B)

		xi_fluor = signal + poisson_noise
		x_fluor.append(xi_fluor)
	
	t+=dt
	x_orig.append(xi)
	time.append(t)

if noise_bkg:
	x = np.matrix(x_fluor)
else:
	x = np.matrix(x_orig)

print(len(x), len(x[0,:]), len(time))

# save results
if save:
	np.savetxt("cellcycle-"+label+".csv", x, delimiter=",")

# Plot the solutions
fig, axs = plt.subplots(2)
plots = []
labels = ['CycA', 'CycB', 'CycE', 'APCP', 'BCKI',  'Cdc20A', 'Cdc20i', \
	     'Cdh1', 'CKI',  'Mass', 'pB',   'pBCKI', 'TriA',   'TriE']

selected = [0,1,2,3,5,6,7,8,9] #,10,12,13]
x = x[:,selected]
labels = np.array(labels)[selected]
print(labels)

prev_idx = 0
start_i = 0
for i in range(len(labels)):
	idx = floor(i/5)
	axs[idx].plot(time,x[:,i])
	if idx != prev_idx:
		axs[idx-1].legend(labels[start_i:i])
		start_i = i
	elif i == len(labels)-1:
		axs[idx].legend(labels[start_i:i+1])
	prev_idx = idx

# axs[idx].plot(time, cycD)

plt.xlabel('t')
# plt.ylabel('Population')
plt.show()
