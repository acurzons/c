import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# delta = 0.43
# E_max = 500.     #TeV
# t_sedov = 100*365*24*3600   #seconds
# D_0 = 1e26      #cm^2s^-1
# E_D = 0.01      #TeV
#
# t = np.arange(0,2001,1)
# ts = t*100*365*24*3600
# E = np.arange(100,501,1)
#
# D = []        #diffusion coefficient
# Chi = []      #escape time
# Rd = []       #diffusion Distance
#
# for i in range(len(E)):
#     Rd.append([])
#     D.append(D_0 * np.sqrt(E[i]/E_D))
#     Chi.append(t_sedov * (E[i]/E_max)**(-1./delta))
#     for j in range(len(ts)):
#         if ts[j] <= Chi[i]:
#             Rd[i].append(0.1)
#         if ts[j] > Chi[i]:
#             Rd[i].append(np.sqrt(4 * D[i] * (ts[j] - Chi[i]))/3.1e19)
#
#
#
# plt.figure(figsize = (9, 9))
# plt.axis([0,2000, 0,30])
# plt.plot(t, Rd[0], '-', color = 'c', markersize = 4,  label = "$\mathrm{E = 100 TeV}$")
# plt.plot(t, Rd[50], '-', color = 'k', markersize = 4,  label = "$\mathrm{E = 150 TeV}$")
# plt.plot(t, Rd[100], '-', color = 'b', markersize = 4,  label = "$\mathrm{E = 200 TeV}$")
# plt.plot(t, Rd[150], '-', color = 'r', markersize = 4,  label = "$\mathrm{E = 250 TeV}$")
# plt.plot(t, Rd[200], '-', color = 'g', markersize = 4,  label = "$\mathrm{E = 300 TeV}$")
# plt.plot(t, Rd[400], '-', color = 'g', markersize = 4,  label = "$\mathrm{E = 500 TeV}$")
# plt.legend(loc = 2)
# plt.xlabel("Time (s)")
# plt.ylabel("Diffusion Distance (kpc)")
# plt.show()
# plt.close()

c = 3e10                     #cms^-1
mbtocm = 1.0e-27
cmtopc = 3.24e-19
delta = 0.43
n_ave = 130.0                     #density
t = 1600*365*24*3600        #seconds
t_sedov = 100*365*24*3600   #seconds
D0 = 1.0e26                 #cm^2s^-1
ED = 0.01                   #TeV
eff = 0.3                   #efficiency
ESNR = 6.25e50              #TeV
Ep = np.logspace(-3,2.7,50)    #TeV
Egamma = np.logspace(-1,2,4) #TeV
Emax = np.amax(Ep)
Emin = np.amin(Ep)
x = np.logspace(15,21,30)    #centimeters
y = np.logspace(15,21,30)    #centimeters
z = np.logspace(15,21,30)    #centimeters
xneg = -np.flip(x,0)
yneg = -np.flip(y,0)
zneg = -np.flip(z,0)

x = np.concatenate((xneg, x))
y = np.concatenate((yneg, y))
z = np.concatenate((yneg, z))
#
# n = []
# for i in range(len(x)):
#     n.append([])
#     for j in range(len(y)):
#         if (x[i] < -1.0*3.2 / 3.24e-19) and (x[i] >= -2.0*3.2 / 3.24e-19) and (y[j] < 3.0*3.2 / 3.24e-19) and (y[j] >= 2.0*3.2 / 3.24e-19):
#             n[i].append(156.4)
#             reg = 1
#         elif (x[i] < 0.0*3.2 / 3.24e-19) and (x[i] >= -1.0*3.2 / 3.24e-19) and (y[j] < 3.0*3.2 / 3.24e-19) and (y[j] >= 2.0*3.2 / 3.24e-19):
#             n[i].append(151.5)
#             reg = 2
#         elif (x[i] < 1.0*3.2 / 3.24e-19) and (x[i] >= 0.0*3.2 / 3.24e-19) and (y[j] < 3.0*3.2 / 3.24e-19) and (y[j] >= 2.0*3.2 / 3.24e-19):
#             n[i].append(180.6)
#             reg = 3
#         elif (x[i] < 2.0*3.2 / 3.24e-19) and (x[i] >= 1.0*3.2 / 3.24e-19) and (y[j] < 3.0*3.2 / 3.24e-19) and (y[j] >= 2.0*3.2 / 3.24e-19):
#             n[i].append(112.2)
#             reg = 4
#         elif (x[i] < -2.0*3.2 / 3.24e-19) and (x[i] >= -3.0*3.2 / 3.24e-19) and (y[j] < 2.0*3.2 / 3.24e-19) and (y[j] >= 1.0*3.2 / 3.24e-19):
#             n[i].append(152.3)
#             reg = 5
#         elif (x[i] < -1.0*3.2 / 3.24e-19) and (x[i] >= -2.0*3.2 / 3.24e-19) and (y[j] < 2.0*3.2 / 3.24e-19) and (y[j] >= 1.0*3.2 / 3.24e-19):
#             n[i].append(109.7)
#             reg = 6
#         elif (x[i] < 0.0*3.2 / 3.24e-19) and (x[i] >= -1.0*3.2 / 3.24e-19) and (y[j] < 2.0*3.2 / 3.24e-19) and (y[j] >= 1.0*3.2 / 3.24e-19):
#             n[i].append(159.0)
#             reg = 7
#         elif (x[i] < 1.0*3.2 / 3.24e-19) and (x[i] >= 0.0*3.2 / 3.24e-19) and (y[j] < 2.0*3.2 / 3.24e-19) and (y[j] >= 1.0*3.2 / 3.24e-19):
#             n[i].append(162.3)
#             reg = 8
#         elif (x[i] < 2.0*3.2 / 3.24e-19) and (x[i] >= 1.0*3.2 / 3.24e-19) and (y[j] < 2.0*3.2 / 3.24e-19) and (y[j] >= 1.0*3.2 / 3.24e-19):
#             n[i].append(142.6)
#             reg = 9
#         elif (x[i] < 3.0*3.2 / 3.24e-19) and (x[i] >= 2.0*3.2 / 3.24e-19) and (y[j] < 2.0*3.2 / 3.24e-19) and (y[j] >= 1.0*3.2 / 3.24e-19):
#             n[i].append(101.9)
#             reg = 10
#         elif (x[i] < -2.0*3.2 / 3.24e-19) and (x[i] >= -3.0*3.2 / 3.24e-19) and (y[j] < 1.0*3.2 / 3.24e-19) and (y[j] >= 0.0*3.2 / 3.24e-19):
#             n[i].append(75.8)
#             reg = 11
#         elif (x[i] < -1.0*3.2 / 3.24e-19) and (x[i] >= -2.0*3.2 / 3.24e-19) and (y[j] < 1.0*3.2 / 3.24e-19) and (y[j] >= 0.0*3.2 / 3.24e-19):
#             n[i].append(73.6)
#             reg = 12
#         elif (x[i] < 0.0*3.2 / 3.24e-19) and (x[i] >= -1.0*3.2 / 3.24e-19) and (y[j] < 1.0*3.2 / 3.24e-19) and (y[j] >= 0.0*3.2 / 3.24e-19):
#             n[i].append(100.6)
#             reg = 13
#         elif (x[i] < 1.0* 3.2 / 3.24e-19) and (x[i] >= 0.0) and (y[j] < 1.0 * 3.2 / 3.24e-19) and (y[j] >= 0.0):
#             n[i].append(157.4)
#             reg = 14
#         elif (x[i] < 2.0*3.2 / 3.24e-19) and (x[i] >= 3.2 / 3.24e-19) and (y[j] < 3.2 / 3.24e-19) and (y[j] >= 0.0):
#             n[i].append(106.5)
#             reg = 15
#         elif (x[i] < 3.0*3.2 / 3.24e-19) and (x[i] >= 2.0*3.2 / 3.24e-19) and (y[j] < 3.2 / cmtopc) and (y[j] >= 0.0):
#             n[i].append(70.8)
#             reg = 16
#         elif (x[i] < -2.0*3.2 / cmtopc) and (x[i] >= -3.0*3.2 / cmtopc) and (y[j] < 0.0*3.2 / cmtopc) and (y[j] >= -1.0*3.2 / cmtopc):
#             n[i].append(61.4)
#             reg = 17
#         elif (x[i] < -1.0*3.2 / cmtopc) and (x[i] >= -2.0*3.2 / cmtopc) and (y[j] < 0.0*3.2 / cmtopc) and (y[j] >= -1.0*3.2 / cmtopc):
#             n[i].append(65.3)
#             reg = 18
#         elif (x[i] < 0.0*3.2 / cmtopc) and (x[i] >= -1.0*3.2 / cmtopc) and (y[j] < 0.0*3.2 / cmtopc) and (y[j] >= -1.0*3.2 / cmtopc):
#             n[i].append(78.4)
#             reg = 19
#         elif (x[i] < 1.0*3.2 / cmtopc) and (x[i] >= 0.0*3.2 / cmtopc) and (y[j] < 0.0*3.2 / cmtopc) and (y[j] >= -1.0*3.2 / cmtopc):
#             n[i].append(97.9)
#             reg = 20
#         elif (x[i] < 2.0*3.2 / cmtopc) and (x[i] >= 1.0*3.2 / cmtopc) and (y[j] < 0.0*3.2 / cmtopc) and (y[j] >= -1.0*3.2 / cmtopc):
#             n[i].append(104.9)
#             reg = 21
#         elif (x[i] < 3.0*3.2 / cmtopc) and (x[i] >= 2.0*3.2 / cmtopc) and (y[j] < 0.0*3.2 / cmtopc) and (y[j] >= -1.0*3.2 / cmtopc):
#             n[i].append(97.0)
#             reg = 22
#         elif (x[i] < -2.0*3.2 / cmtopc) and (x[i] >= -3.0*3.2 / cmtopc) and (y[j] < -1.0*3.2 / cmtopc) and (y[j] >= -2.0*3.2 / cmtopc):
#             n[i].append(86.0)
#             reg = 23
#         elif (x[i] < -1.0*3.2 / cmtopc) and (x[i] >= -2.0*3.2 / cmtopc) and (y[j] < -1.0*3.2 / cmtopc) and (y[j] >= -2.0*3.2 / cmtopc):
#             n[i].append(68.3)
#             reg = 24
#         elif (x[i] < 0.0*3.2 / cmtopc) and (x[i] >= -1.0*3.2 / cmtopc) and (y[j] < -1.0*3.2 / cmtopc) and (y[j] >= -2.0*3.2 / cmtopc):
#             n[i].append(71.5)
#             reg = 25
#         elif (x[i] < 1.0*3.2 / cmtopc) and (x[i] >= 0.0*3.2 / cmtopc) and (y[j] < -1.0*3.2 / cmtopc) and (y[j] >= -2.0*3.2 / cmtopc):
#             n[i].append(69.9)
#             reg = 26
#         elif (x[i] < 2.0*3.2 / cmtopc) and (x[i] >= 1.0*3.2 / cmtopc) and (y[j] < -1.0*3.2 / cmtopc) and (y[j] >= -2.0*3.2 / cmtopc):
#             n[i].append(101.9)
#             reg = 27
#         elif (x[i] < 0.0*3.2 / cmtopc) and (x[i] >= -1.0*3.2 / cmtopc) and (y[j] < -2.0*3.2 / cmtopc) and (y[j] >= -3.0*3.2 / cmtopc):
#             n[i].append(78.8)
#             reg = 28
#         elif (x[i] < 1.0*3.2 / cmtopc) and (x[i] >= 0.0*3.2 / cmtopc) and (y[j] < -2.0*3.2 / cmtopc) and (y[j] >= -3.0*3.2 / cmtopc):
#             n[i].append(78.9)
#             reg = 29
#         else:
#             n[i].append(0)
#



rhoCR = []
for i in range(len(Ep)):
    print i
    rhoCR.append([])
    D = D0*np.sqrt(Ep[i]/ED)
    tescape = t_sedov*(Ep[i]/Emax)**(-1.0/delta)
    rhoBG = 7.5e-10*(Ep[i]/0.001)**(-2.75)
    if t >= tescape:
        Rd = np.sqrt(4 * D * (t - tescape))
    else:
        Rd = 0
    for j in range(len(x)):
        rhoCR[i].append([])
        for k in range(len(y)):
            rhoCR[i][j].append([])
            for h in range(len(z)):
                d = np.sqrt(x[j]**2 + y[k]**2 + z[h]**2)
                if Rd != 0:
                    rhoSNR = eff * ESNR / (np.pi**(1.5) * np.log(Emax/Emin)) * np.exp(-(d/Rd)**2) / (Rd**3 * Ep[i]**2)
                else:
                    rhoSNR = 0
                rhoCR[i][j][k].append(rhoSNR + rhoBG)

e = None
e = []
for g in range(len(Egamma)):
    print g
    e.append([])
    for i in range(len(x)):
        print i
        e[g].append([])
        for j in range(len(y)):
            e[g][i].append([])
            for k in range(len(z)):
                for h in range(len(Ep)):
                    L = np.log(Ep[h])
                    sigma = (34.3 + 1.88*L + 0.25*L**2) * mbtocm   #cm^2
                    B = 1.30 + 0.14*L + 0.011*L**2
                    beta = (1.79 + 0.11*L + 0.008*L**2)**(-1)
                    kk = (0.801 + 0.049*L + 0.014*L**2)**(-1)
                    if h == 0:
                        e[g][i][j].append(0)
                    else:
                        if Ep[h] > Egamma[g]:
                            xx = Egamma[g]/Ep[h]
                            Fgamma = (B * np.log(xx)/xx * ((1. - xx**beta)/(1. + kk*xx**beta*(1. - xx**beta)))**4.0 * (1.0/np.log(xx) -
                            4.0*beta*xx**beta/(1. - xx**beta) - 4*kk*beta*xx**beta*(1. - 2.0*xx**beta)/(1. + kk*xx**beta*(1. - xx**beta))))
                            dEp = Ep[h] - Ep[h - 1]
                            int = sigma * Fgamma / Ep[h] * c * n_ave * rhoCR[h][i][j][k] * dEp
                        else:
                            int = 0
                        e[g][i][j][k] = e[g][i][j][k] + int

flux = None
flux = []
for h in range(len(Egamma)):
    print h
    flux.append([])
    for i in range(len(x)):
        flux[h].append([])
        for j in range(len(y)):
            flux[h][i].append(0)
            for k in range(len(z)):
                    sum = e[h][i][j][k] / 4.0 / np.pi
                    flux[h][i][j] = flux[h][i][j] + sum

# np.savetxt("gamma_ray_emissivity.out", e, delimiter= '  ')



xplot = x * 3.24e-19  #pc
yplot = y * 3.24e-19

print "Hello", Ep[15]


X, Y = np.meshgrid(xplot, yplot)
f, ( (ax1,ax2) , (ax3,ax4)) = plt.subplots(2,2, sharex = True, sharey = True)
cont1 = ax1.contourf(X, Y, np.log10(flux[0]), 10)
ax1.set_title("E$_{\gamma}$ = 0.1 TeV")
ax1.set_xlim([-50,50])
ax1.set_ylim([-50,50])
plt.colorbar(cont1, ax=ax1, format = '%i')
cont2 = ax2.contourf(X, Y, np.log10(flux[1]), 1000)
ax2.set_title("E$_{\gamma}$ = 1 TeV")
plt.colorbar(cont2, ax=ax2, format = '%i')
cont3 = ax3.contourf(X, Y, np.log10(flux[2]), 1000)
ax3.set_title("E$_{\gamma}$ = 10 TeV")
plt.colorbar(cont3, ax=ax3, format = '%i')
cont4 = ax4.contourf(X, Y, np.log10(flux[3]), 900)
ax4.set_title("E$_{\gamma}$ = 100 TeV")
plt.colorbar(cont4, ax=ax4, format = '%i')
plt.show()
plt.close()


# plt.figure(figsize = (8,8))
# plt.title("Diffusion")
# plt.contourf(X, Y, np.log10(e[8]), 800)
# plt.colorbar(format = '%i')
# plt.axis([-3, 3, -3, 3])
# plt.xlabel(r'$\mathrm{X\: (kpc)}$', size = 15)
# plt.ylabel(r'$\mathrm{Y\: (kpc)}$', size = 15)
# plt.show()
# plt.close()
