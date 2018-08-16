import numpy as np
import matplotlib.pyplot as plt
from fit_str8_line import str8_line
import os.path
from scipy.stats import spearmanr

regions = 29

for i in range(regions):
    #filename = "correlation_study_data/suzaku_mopHI_data/reg" + str(i + 1) + ".txt"
    filename = "correlation_study_data/suzaku_mopHI_data/NRreg" + str(i + 1) + ".txt"
    if os.path.isfile(filename) == False:
        continue
    data = np.loadtxt(filename, skiprows = 1)

    xray_ec = data[:,1] #gamma ray excess Counts
    xray_ec_err = data[:,3]
    mHI_cd = data[:,0]/1.0e21 #mopra + SGPS column Density
    mHI_cd_err = data[:,2]/1.0e21

    N = len(data[:,0])

    x = np.arange(min(mHI_cd) - 5, max(mHI_cd) + 5, 0.1)

    #calculate y = a + bx and Pearsons correlation coef
    beta, variance, var_y, pearsons, pvalue = str8_line(x = mHI_cd, y = xray_ec, xerr = mHI_cd_err, yerr = xray_ec_err, line_x = x)
    y = beta[0] + beta[1] * x

    spearman, pvalues = spearmanr(mHI_cd,xray_ec)
    print i+1, pearsons, pvalue, spearman, pvalues

    if pvalue < 0.001:
        label = "$\\rho$ = %.2f, p-value = < 0.001" %(pearsons)
    else:
        label = "$\\rho$ = %.2f, p-value = %.3f" %(pearsons, pvalue)

    plt.figure()
    plt.errorbar(mHI_cd, xray_ec, xerr = mHI_cd_err, yerr = xray_ec_err, fmt = "o")
    plt.plot(x,y, label = label, color = 'k')
    plt.fill_between(x, y + var_y**0.5, y - var_y**0.5, color = 'k', alpha = 0.1)
    plt.xlim((np.min(mHI_cd) - np.max(mHI_cd_err),np.max(mHI_cd) + np.max(mHI_cd_err)))
    plt.ylim((0,np.max(xray_ec)*1.3))
    plt.legend(loc=2)
    plt.xlabel("$N_p$ Column Density (x10$^{21}$ cm$^{-2}$)")
    plt.ylabel("X-ray Excess Counts")
    plt.title("Region %d" %(i+1))
    #plt.savefig("correlation_study_data/plots2/suzaku_mHI_reg" + str(i + 1) + ".ps")
    plt.savefig("correlation_study_data/plots2/suzakuNR_mHI_reg" + str(i + 1) + ".ps")
    plt.clf()
    plt.close()
