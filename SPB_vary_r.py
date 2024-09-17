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
r = -0.5 # scattering parameter. -0.5 = APS/ADP, 1.5 = IMP, 0.5 = POP, -0.5 = point defect
R = 8.314 # J / mol K (This is the gas constant)
h = 6.62607015e-34 # This is plancks contant in J * s
hbar = h/(2*math.pi)
h_eV = 4.135667696E-15 # Plancks constant in eV * s
integralLowerLimit = 0  # Set the lower limit of the fermi integrals
integralUpperLimit = 100  # Set the upper limit of the fermi integrals

temperature = 323 #Kelvin
effm = 1.28 #i.e. 1.28, or some reasonable float value between 0.001 and 3

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



seeblist = []
etalist = []
carrierlist = []


###############p-TYPE CARRIERS####################################################
for i in np.arange(0,25,0.25):
    eta = i
    etalist.append(i)
    seeblist.append(seebeckfunction(i))
    carrierlist.append(carrierfunction(i, effm)/1E20)

fig, ax = plt.subplots(figsize=(7,7))

#plot calculated SPB data for a given effective mass and eta range (lines 67-71)
plt.plot(carrierlist, seeblist, label = f'SPB, $m^*$= {effm}$m_e$', color = 'r')


data = pd.read_excel('GeTe_data.xlsx')
# data = undoped.fillna(0)
# data = undoped.values

df = pd.DataFrame(data = data, columns= ['sample', 'carrier', 'seebeck', 'percentage indium'])

#### Plot lab data from GeTe_data.xlsx

x = df['carrier']/1E20
y = df['seebeck']
tag = df['sample']
size = df['percentage indium']
for i in range(len(x)):
    if 'In' not in tag[i]:
        plt.scatter(x[i],y[i],edgecolors='black',marker='D',label = tag[i], s=100)
    if 'In' in tag[i]:
        plt.scatter(x[i],y[i],color = 'magenta',edgecolors='black',marker='o',label = tag[i], s = (size[i]*1)*100)


ax.yaxis.set_minor_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.tick_params(axis='both', which='major', labelsize=16)

plt.title(f"SPB {temperature}K $m*$ = {effm}$m_e$  r = {r}", fontsize=20)
# # # plt.scatter(etalist, KFIlist)
plt.legend(frameon= True, fontsize =16, ncol=1)
# # plt.xlabel('$\\eta$ = $E_F$/$k_BT$', fontsize = 14)
plt.ylabel('Seebeck coefficient ($\\mu$V/K)', fontsize = 20)


plt.xlabel('Carrier concentration (10$^{20}$ cm$^{-3}$)', fontsize = 20)
# plt.xscale('log')
plt.tight_layout(pad=.5)
ax.set_aspect(1.0/ax.get_data_ratio())

timestamp = str(time.time()).replace(':', '-')
# directory = 'C:\\Users\\Admin\\Documents\\Resonant doping\\transport data\\plots\\Pisarenko plots\\SnTe\\'
# plt.savefig(directory + 'Kane-TBK_SnTe_simple_' + timestamp + '.png')
plt.savefig(f'GeTe_SBP_effm{effm}_r_{r}_{timestamp}.png')
plt.show()





