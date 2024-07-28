import numpy as np
from scipy.integrate import quad
import scipy.interpolate
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
import math
import pandas as pd
# import excel
import matplotlib.style as style
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rcParams
import time
plt.rcParams.update({'font.family':'Arial'})
# There is no need to mess with them unless you really know what you're doing.
m = 9.10938356E-31 # mass of electron in Kg
kb = 1.3806505E-23 # boltzmann constant (J/K)
kb_eV = 8.6173E-5# boltzmann constant (eV/K)
e = 1.60217653E-19 # Charge on electron (Coulomb)
r = -0.5 # scattering parameter. -0.5 = APS, 1.5 = IMP, 0.5 = POP, -0.5 = point defect
R = 8.314 # J / mol K (This is the gas constant)
h = 6.62607015e-34 # This is plancks contant in J * s
hbar = h/(2*math.pi)
h_eV = 4.135667696E-15 # Plancks constant in eV * s
integralLowerLimit = 0  # Set the lower limit of the fermi integrals
integralUpperLimit = 100  # Set the upper limit of the fermi integrals

temperature = 323 #Kelvin
E_gap = .18 #eV; bandgap of SnTe
alpha = (kb_eV*temperature)/E_gap #non-parabolicity factor; set to zero if band is parabolic
b = 4 # ratio of mobility between light and heavy bands
K = 4 # ratio of effective mass along parallel vs perpendicular directions; anisotropy of VB
deltav = 0.35/(kb_eV*temperature) #eV; this is band offset between sigma and L bands in SnTe
effmheavy = 1.92
effmlight = 0.168
xi = (effmheavy/effmlight)*(0.25) #ratio of effectivemasses * acoustic deformation potential ratio (of heavy/light)squared
from matplotlib import cm # import the colormap

################### old SPB #################################
def seebeckfunction(eta):
    def FermiIntegral(epsilon, eta, j):
        return epsilon ** j / (1 + np.exp(epsilon - eta))
    # eta = -3
    j = r + 3/2
    Itop = quad(FermiIntegral, 0 , 100, args = (eta , j))
    j = r + 1/2
    Ibot = quad(FermiIntegral, 0, 100, args = (eta,j))
    coefficient = (r + 5/2)/(r + 3/2)
    seebeck = round(1E6 * kb/e * (coefficient * (Itop[0]/Ibot[0])- eta),3)
    return seebeck
    # print(str(round(seebeck,3)))

def carrierfunction(eta, effm):
    def FermiIntegral(epsilon, eta, j):
        return epsilon ** j / (1 + np.exp(epsilon - eta))
    # eta = -3
    j = 0.5 + r
    Itop = quad(FermiIntegral, 0 , 100, args = (eta , j))[0]
    # temperature = 273 + 50
    j= 2 * r + 0.5
    Ibot = quad(FermiIntegral, 0, 100, args = (eta, j))[0]
    rcoefficient = ((r + 3/2) ** 2)/(2*r + 3/2)
    numerator = (2 * effm* m * kb * temperature)**(3/2)
    denominator = (3 * math.pi**2) * hbar ** 3
    carrier = (1E-6 * rcoefficient * numerator * Itop**2) / (denominator * Ibot)
    return carrier


################### SKB #################################
def KaneFermiIntegral(eta,n,m,l,alpha):
    def dervivFDdist(epsilon, eta):
        return np.exp(epsilon-eta)/((np.exp(epsilon-eta)+1)**2)
    def nterm(epsilon, eta, n):
        return -dervivFDdist(epsilon,eta)*(epsilon**n)
    def mterm(epsilon, eta, m):
        return (epsilon + alpha*epsilon**2)**m
    def lterm(epsilon, eta, alpha, l):
        # return ((1 + 2*alpha*epsilon)**2 + 2)**(l/2)
        return (1 + 2*alpha*epsilon)**l # Kamil's edit 1/12/24
    def KaneFermiIntegrand(epsilon, eta, alpha, n, m, l):
        return nterm(epsilon, eta, n) * mterm(epsilon, eta, m) * lterm(epsilon, eta, alpha, l)
    Kaneintegral = quad(KaneFermiIntegrand, 0,100, args=(eta, alpha, n, m, l))[0]
    return Kaneintegral

#Hall factors, from Wiendlocha et al, Residual resistivity as a ....
def simpleFermiIntegral(eta,j):
    def Fermiintegrand(epsilon, eta, j):
        return epsilon ** j / (1 + np.exp(epsilon - eta))
    integratedFermi = quad(Fermiintegrand, 0 , 100, args = (eta , j))[0]
    return integratedFermi

def A_hh(eta):
    heavyHallFactor = 1.5 * simpleFermiIntegral(eta,1/2) * (simpleFermiIntegral(eta, -1/2)/(2*(simpleFermiIntegral(eta,0)**2)))
    return heavyHallFactor

def A_lh(eta):
    top = 3*K*(K+2)* KaneFermiIntegral(eta, 0, 1/2, -4, alpha=float(alpha)) * KaneFermiIntegral(eta, 0, 3/2, 0, alpha=float(alpha))
    bottom = ((2*K + 1)**2) * (KaneFermiIntegral(eta,0,1,-2, alpha=float(alpha))**2)
    lightHallFactor = top/bottom
    return lightHallFactor

def kanelightseebeckfunction(eta):
    top = KaneFermiIntegral(eta, 1, 1, -2, alpha=float(alpha))
    bottom = KaneFermiIntegral(eta, 0, 1, -2, alpha=float(alpha))
    # print(bottom)
    kaneseebeck = round(1E6*(kb/e * ((top/bottom)-eta)), 3)
    return kaneseebeck

def kaneheavyseebeckfunction(potato): #potato should be set equal to eta-deltav
    top = KaneFermiIntegral(potato, 1, 1, -2, 0)
    bottom = KaneFermiIntegral(potato, 0, 1, -2, 0)
    kaneheavyseebeck = round(1E6*(kb/e * ((top/bottom)-(potato))), 3)
    return kaneheavyseebeck


## make sure to change temperature (global variable) as needed
def kanelightcarrierfunction(eta, effm):
    top = (2*effm*m*kb*temperature)**(3/2)
    bottom = (3*(np.pi**2))*(hbar**3)
    # print(KaneFermiIntegral(eta, 1, 1, -2)/ KaneFermiIntegral(eta, 0, 1, -2))
    # print(KaneFermiIntegral(eta, 0,3/2,0))
    kanecarrier = 1E-6*((top/bottom) * KaneFermiIntegral(eta, 0, 1.5, 0, alpha=float(alpha)))
    return -kanecarrier #ugh I hate this hard coding


def kaneheavycarrierfunction(potato, effm): #potato should be set equal to eta-deltav
    top = (2*effm*m*kb*temperature)**(3/2)
    bottom = (3*(np.pi**2))*(hbar**3)
    # print(KaneFermiIntegral(eta, 1, 1, -2)/ KaneFermiIntegral(eta, 0, 1, -2))
    # print(KaneFermiIntegral(eta, 0,3/2,0))
    kanecarrier = 1E-6*((top/bottom) * KaneFermiIntegral(potato, 0, 1.5, 0, alpha= 0))
    return -kanecarrier #ugh I hate this hard coding

def TKBseebeck(eta):
    top = xi*(KaneFermiIntegral(eta, 1, 1, -2, alpha=float(alpha))-eta*KaneFermiIntegral(eta,0,1,-2, alpha=float(alpha))) +\
          (KaneFermiIntegral(eta-deltav, 1,1,-2,0) - (eta-deltav)*KaneFermiIntegral(eta-deltav, 0,1,-2,0))
    bottom = xi*KaneFermiIntegral(eta,0,1,-2, alpha=float(alpha)) + KaneFermiIntegral(eta - deltav,0,1,-2,0)
    kane2bseebeck = 1E6*((kb/e) * (top/bottom))
    return kane2bseebeck


def TKBcarrier(eta):
    pl = kanelightcarrierfunction(eta, effmlight)
    ph = kaneheavycarrierfunction(eta-deltav, effmheavy)
    top = (b*pl + ph)**2
    bottom = (A_lh(eta)*(b**2)*pl + A_hh(eta)*ph)
    totalcarriers = top/bottom
    return totalcarriers




emptylist = []
etalist = []
carrierlist = []
KFIlist = [] #KFI = kane fermi integral
kaneseebecklist = []
kanecarrierlist = []
lightcarriers = []
heavycarriers = []
lightseebeck = []
heavyseebeck = []
lightplusheavycarriers = []
phpluspllist=[]
SPBseebecklist = []
SPBcarrierlist = []
###############p-TYPE CARRIERS####################################################
for i in np.arange(0,25,0.25):
    eta = i
    etalist.append(i)
    emptylist.append(seebeckfunction(i-deltav))
    carrierlist.append(carrierfunction(i-deltav,effmheavy)/1E20) #second input here is your effective mass (band)
    # KFIlist.append(KaneFermiIntegral(i, 1, 1, -2))
    lightcarriers.append(kanelightcarrierfunction(i, effmlight)/1E20)
    heavycarriers.append(kaneheavycarrierfunction(i-deltav, effmheavy)/1E20)
    lightseebeck.append(kanelightseebeckfunction(i))
    heavyseebeck.append((kaneheavyseebeckfunction(i-deltav)))
    kaneseebecklist.append(TKBseebeck(i))
    kanecarrierlist.append(TKBcarrier(i)/1E20)
    SPBcarrierlist.append(carrierfunction(i, 1.28)/1E20)
    SPBseebecklist.append(seebeckfunction(i))
fig, ax = plt.subplots(figsize=(8,8))

# plt.scatter(etalist, kaneseebecklist)


data = pd.read_csv('SGTA_50C_pisarenko.csv')
# data = undoped.fillna(0)
# data = undoped.values

df = pd.DataFrame(data = data, columns= ['Sample number','Carrier', 'Seebeck'])

plt.ylim(0,120)
plt.xlim(0,20)
# if 'In' in df['label']:
#     plt.scatter(df['n']/1E20, df['s'], color = 'magenta')
# else:
#     plt.scatter(df['n'] / 1E20, df['s'], color = 'limegreen')
#
x = df['Carrier']/1E20
y = df['Seebeck']
tag = df['Sample number']
# size = df['size']
colors = cm.viridis(np.linspace(0, 1, len(x)))  # Use the viridis colormap

for i in range(len(x)):
    if 'In' not in tag[i]:
        plt.scatter(x[i],y[i],edgecolors='black',marker='D',label = tag[i], s=100, color=colors[i])
    if 'In' in tag[i]:
        plt.scatter(x[i],y[i],color = 'magenta',edgecolors='black',marker='o', alpha = 0.5, label = tag[i])#, s = size[i]*5)

plt.plot(SPBcarrierlist, SPBseebecklist, color = 'orange', label = 'SPB, $m^*$=1.28$m_e$')
plt.plot(kanecarrierlist, kaneseebecklist, color = 'royalblue', linewidth='2.5',label = f"Two band Kane,\n $m*_h$ = {effmheavy}$m_e$, \n$m*_l$ ={effmlight}$m_e$")#, $\\alpha$=$k_B$T/$E_g$")
# plt.plot(lightcarriers, lightseebeck, linewidth='2',label = f'light Kane band, $m^*$= {effmlight}$m_e$')
# plt.plot(heavycarriers, heavyseebeck, linewidth='2',label = f'heavy parabolic band, $m^*$= {effmheavy}$m_e$')
ax.yaxis.set_minor_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.tick_params(axis='both', which='major', labelsize=16)
# plt.scatter(carrierlist, emptylist, label = f'SPB, $m^*$= {effmheavy}$m_e$')

# plt.scatter(1.05E+20/1E20, 4.25605261, marker='*', s = 100, color= 'yellow', edgecolors='black', label = 'SnTe')
# plt.scatter(1.5,45, marker='D', edgecolors='black', color='magenta',label = '0.25% In doped SnTe')
# plt.scatter(1.6,87, marker='^', edgecolors='black', color='magenta',label = '1.00% In doped SnTe')
plt.annotate('GeTe!', xy=(4.5,49))
# plt.title(f"TBK 300K, $m*_h$ = {effmheavy}, $m*_l$ ={effmlight}", fontsize=20)
# # # plt.scatter(etalist, KFIlist)
plt.legend(frameon= False, fontsize =18, ncol=2)
# # plt.xlabel('$\\eta$ = $E_F$/$k_BT$', fontsize = 14)
plt.ylabel('50$^o$C Seebeck coefficient ($\\mu$V/K)', fontsize = 20)


plt.xlabel('50$^o$C Carrier concentration (10$^{20}$ cm$^{-3}$)', fontsize = 20)
# plt.xscale('log')
plt.tight_layout(pad=.5)
ax.set_aspect(1.0/ax.get_data_ratio())

timestamp = str(time.time()).replace(':', '-')
directory = 'C:\\Users\\Admin\\Documents\\Resonant doping\\transport data\\plots\\Pisarenko plots\\SnTe\\'
# plt.savefig('SGTA_50CPisarenko.png', dpi=1000)
plt.show()
# plt.plot(carrierlist,emptylist, label = 'm*$_{DOS}$ = 0.50, SPB')





