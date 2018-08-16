import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

n_H = 130   # cm^-3
B = 14.26   # micro Gauss
Urad_CMB = 0.26 # eVcm^-3
Urad_farIR = 0.415  # eVcm^-3
E_e = np.logspace(-12,3,100)

filename = "rxj1713_leptonic_regional_parameters.txt"
data = np.loadtxt(filename)
reg = data[:,0]
nH = data[:,1]
B = data[:,3]
Ee_1 = 1 #TeV
Ee_10 = 10
Ee_100 = 100

Ee = [1,10,100]
N = len(data[:,0])


t_IC_CMB_10 = []
t_sync_10 = []
t_IC_CMB_100 = []
t_sync_100 = []

for j in range(3):
    t_pp = []
    t_br = []
    t_IC_CMB= []
    t_sync= []
    for i in range(N):
        t_pp.append(5.3e7 / nH[i])
        t_br.append(4.e7 / nH[i])
        t_IC_CMB.append(3.e5 / Urad_CMB / Ee[j])
        t_sync.append(12.e6 / B[i]**2 / Ee[j])

    fig = plt.figure(figsize = (9,9))
    ax = plt.gca()
    ax.axis([0,30,2e2,5e7])
    ax.set_title("$E_e$ = " + str(Ee[j]) + "TeV", size ="16")
    ax.scatter(reg, t_br, c="blue", label = "Bremsstrahlung")
    ax.scatter(reg, t_sync, c="red", label = "Synchrotron")
    ax.scatter(reg, t_IC_CMB, c="green", label = "Inverse Compton")
    # ax.scatter(reg, t_sync_10, c="red", marker = ",", label = "Synchrotron (10 TeV)")
    # ax.scatter(reg, t_IC_CMB_10, c="green", marker = ",", label = "Inverse Compton (10 TeV)")
    # ax.scatter(reg, t_sync_100, c="red",marker = "+", label = "Synchrotron (100 TeV)")
    # ax.scatter(reg, t_IC_CMB_100, c="green", marker = "+", label = "Inverse Compton (100 TeV)")
    ax.legend(loc = 1)
    ax.set_yscale('log')
    plt.xlabel("Region")
    plt.ylabel("Time (yrs)")
    fig.savefig('coolingtime_plots/energy_' + str(Ee[j]) + 'TeV.png')
    plt.close()


# fig.close()

#
# for i in range(len(E_e)):
#     t_pp.append(5.3e7 * (1./n_H))
#     t_br.append(4e7 * (1./n_H))
#     t_IC_CMB.append(3e5 * (1./Urad_CMB) * (1./E_e[i]))
#     t_IC_farIR.append(3e5 * (1./Urad_farIR) * (1./E_e[i]))
#     t_sync.append(12e6 * (1./B)**2 * (1./E_e[i]))
#
# plt.figure(figsize = (9, 9))
# plt.axis([1e-2,200, 1e2,1e8])
# plt.loglog(E_e, t_pp, '-', color = 'y', markersize = 4,  label = r"$\mathrm{PP}$")
# plt.loglog(E_e, t_br, '-', color = 'c', markersize = 4,  label = r"$\mathrm{Brems}$")
# plt.loglog(E_e, t_IC_CMB, '-', color = 'k', markersize = 4,  label = r"$\mathrm{IC_{CMB}}$")
# plt.loglog(E_e, t_IC_farIR, '-', color = 'b', markersize = 4,  label = r"$\mathrm{IC_{farIR}}$")
# plt.loglog(E_e, t_sync, '-', color = 'r', markersize = 4,  label = r"$\mathrm{Sync}$")
# plt.legend()
# plt.xlabel("Energy (TeV)")
# plt.ylabel("Cooling time (years)")
# plt.show()
# plt.close()
