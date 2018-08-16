import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.optimize import curve_fit
from scipy import stats
import os.path

def korr(x,y,xerr,yerr,plot=False):

    xx = np.random.normal(x,xerr)
    yy = np.random.normal(y,yerr)
    if plot==True:
        plt.errorbar(xx, yy, xerr=xerr, yerr=yerr, fmt='o', color='k')
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', color='g')
        plt.show()
        print(pearsonr(xx,yy))

    return pearsonr(xx,yy)[0]

def korr2(x_mean, x_sigm, y_mean, y_sigm, number):

    x = np.random.normal(x_mean,x_sigm,number)
    y = np.random.normal(x_mean,x_sigm,number)

    return pearsonr(x,y)[0]

#read data
for i in range(29):
    filename = "correlation_study_data/suzaku_mopHI_data/NRreg" + str(i + 1) + ".txt"
    if os.path.isfile(filename) == False:
        continue
    data = np.loadtxt(filename, skiprows = 1)

    gas = data[:,0]
    xray = data[:,1]
    gas_err = data[:,2]
    xray_err = data[:,3]

    #korr(gas,xray,gas_err,xray_err, plot=True)

    k = []
    for i in range(100000):
        k.append(korr(gas,xray,gas_err,xray_err))

    plt.hist(k, bins=81, range=[0,0.25])
    # 3 sigma
    plt.axvline(np.percentile(k, (100.0 - 99.73)/2.0), color='black', ls="--",)
    plt.axvline(np.percentile(k, (100.0 - (100.0 - 99.73)/2.0)), color='black', ls="--")

    # 5 sigma
    plt.axvline(np.percentile(k, (100.0 - 99.99994)/2.0), color='black', ls=":")
    plt.axvline(np.percentile(k, (100.0 - (100.0 - 99.99994)/2.0)), color='black', ls=":")

    # correlation of original data set
    plt.axvline(pearsonr(xray, gas)[0], color='black')

    #plt.show()


    x_mean = np.mean(gas)
    x_sigm = np.std(gas)
    y_mean = np.mean(xray)
    y_sigm = np.std(xray)
    number = len(gas)

    #print(x_mean, x_sigm, y_mean, y_sigm)


    kk = []
    for i in range(100000):
        kk.append(korr2(x_mean, x_sigm, y_mean, y_sigm, number))

    plt.hist(np.abs(kk), bins=81, range=[0,1])

    print('2 sigma', np.percentile(np.abs(kk),95.449))
    print('3 sigma', np.percentile(np.abs(kk),99.73))
    print('5 sigma', np.percentile(np.abs(kk),99.99994))

    # plt.show()

    print(np.abs(pearsonr(xray,gas)[0]))
    print(np.percentile(np.abs(kk),100-pearsonr(xray,gas)[1]*100.))

    print(i+1, np.percentile(np.abs(kk),95.449, np.percentile(np.abs(kk),100-pearsonr(xray,gas)[1]*100.))
