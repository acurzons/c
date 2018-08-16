import numpy as np
import matplotlib.pyplot as plt
from fit_str8_line import str8_line

regions = 29

for i in range(regions):
    filename = "correlation_study_data/gamma_Np_data/reg" + str(i + 1) + ".txt"
    data = np.loadtxt(filename, skiprows = 1)

    gamma_ec = data[:,0] #gamma ray excess Counts
    gamma_ec_err = data[:,2]
    nHI_cd = data[:,1] #nanten + SGPS column Density
    nHI_cd_err = data[:,3]

    N = len(data[:,0])

    x = np.arange(min(nHI_cd) - 5, 30, 0.1)

    #calculate y = a + bx
    beta, variance, var_y, pearsons, pvalue= str8_line(nHI_cd, gamma_ec, nHI_cd_err, gamma_ec_err, x)
    y = beta[0] + beta[1] * x
    print i+1, pearsons, pvalue

    if pvalue < 0.001:
        label = "$\\rho$ = %.2f, p-value = < 0.001" %(pearsons)
    else:
        label = "$\\rho$ = %.2f, p-value = %.3f" %(pearsons, pvalue)

    plt.figure()
    plt.errorbar(nHI_cd, gamma_ec, xerr = nHI_cd_err, yerr = gamma_ec_err, fmt = "o")
    plt.plot(x,y, label = label, color = 'k')
    plt.fill_between(x, y + var_y**0.5, y - var_y**0.5, color = 'k', alpha = 0.01)
    plt.xlim((np.min(nHI_cd) - np.max(nHI_cd_err)*1.1,np.max(nHI_cd) + np.max(nHI_cd_err)*1.1))
    plt.ylim((0,np.max(gamma_ec) + np.max(gamma_ec_err)*1.6))
    plt.legend(loc=2)
    plt.xlabel("$N_p$ Column Density (x10$^{21}$ cm$^{-2}$)")
    plt.ylabel("$\gamma$-ray Excess Counts")
    plt.title("Region %d" %(i+1))
    plt.savefig("correlation_study_data/plots2/gamma_nHI_reg" + str(i + 1) + ".ps")
    plt.clf()
    plt.close()
