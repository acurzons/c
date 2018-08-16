import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr

#Define
# x_model = np.logspace(20,23,1000)
x_model = np.arange(0,100,1)


#User inputs
# ifile = input("What is the file name?\n")
correlation_type = input("What data are we comparing?\n 1.Mopra+HI Hess\n 2.Nanten+HI Hess\n 3.Mopra vs Naten\n 4.Mopra+HI vs Suzaku\n")
column = input("What column is x data?\n 1.Left\n 2.Right\n")
noise_removed_data = 0
if correlation_type == 1:
    ifile = "gamma_mopraHI_data"
elif correlation_type == 2:
    ifile = "gamma_Np_data"
elif correlation_type == 3:
    ifile = "nanten_mopra_data"
    noise_removed_data = input("Do you want to use the dataset with noise removed?\n 1.Yes\n 2.No\n")
elif correlation_type == 4:
    ifile = "suzaku_mopHI_data"
    x_model = np.logspace(-8,-3,500)
    # x_model = np.logspace(20,24,500)

print "\n\nReg. || rho p rs p sig? sigb siga covab"
print "==============================================================="

# Load the line of best fit data computed in c
# cfilename = "correlation_study_data/gamma_mopra.txt"
# cdata = np.loadtxt("correlation_study_data/gamma_mopra.txt")
# a = cdata[:,0]
# b = cdata[:,1]
# a_err = cdata[:,2]
# b_err = cdata[:,3]
# covab = cdata[:,7]

if correlation_type == 1 or correlation_type == 2:
    no_regions = 29
elif correlation_type == 3 or correlation_type == 4:
    no_regions = 1
if noise_removed_data == 1:
    ifile_n = "_NR"
else:
    ifile_n = ""

for i in range(no_regions):
    #Load the data into arrays
    if correlation_type == 1 or correlation_type == 2:
        filename = "correlation_study_data/" + str(ifile) + "/reg" + str(i + 1) + ".txt"
    elif correlation_type == 3:
        filename = "correlation_study_data/" + str(ifile) + "/all_data" + str(ifile_n) + ".txt"
    elif correlation_type == 4:
        filename = "correlation_study_data/" + str(ifile) + "/reg3.txt"
    data = np.loadtxt(filename, skiprows = 1)

    N = len(data[:,0])

    #Determine which column to use for x/y data
    if column == 1:
        x = data[:,0]
        y = data[:,1]
        x_err = data[:,2]
        y_err = data[:,3]
    elif column == 2:
        y = data[:,0]
        x = data[:,1]
        y_err = data[:,2]
        x_err = data[:,3]


    #Create the lines of best fit (and confidence bands) for data
    p,var = np.polyfit(x, y, 1, w=y_err**0.5, cov=True)
    sigma = np.sqrt(abs(var))
    line = p[0]*x_model + p[1]
    sig_line = np.sqrt(abs((sigma[0,0] * x_model)**2 + sigma[1,1]**2 + 2 * x_model * sigma[0,1]))
    upper_band = line + sig_line
    lower_band = line - sig_line


    #Calculate the pearson correlation coefficient and determine whether it is significant or not
    rho = pearsonr(x,y)
    if rho[1] <= 0.05:
        significant = "yes"
    else:
        significant = "no"

    #calculate the spearman correlaiton coefficient and p-value
    rs  = spearmanr(x,y)
    if rs[1] <= 0.05:
        significants = "yes"
    else:
        significants = "no"

    #print '{} || {:03.2f} {:04.3f} {} {:03.2f} {:04.3f} {} {:04.3f} {:04.3f} {:04.2f}'.format(i+1, rho[0], rho[1], significant, rs[0], rs[1], significants, sigma[0,0], sigma[1,1], sigma[0,1])
    print '{:03.2f} {:04.3f}'.format(p[0], p[1])

    #Plot each region
    plt.figure()
    plt.errorbar(x, y, xerr = x_err, yerr = y_err, fmt = "o")
    if rho[1] < 0.001:
        label = "$\\rho$ = %.2f, p-value = < 0.001" %(rho[0])
    else:
        label = "$\\rho$ = %.2f, p-value = %.3f" %(rho[0], rho[1])
    if correlation_type == 3:
        label = "Slope = %.2f $\pm$ %.2f\nIntercept = %.3f" %(p[0], sigma[0][0], p[1])
    plt.plot(x_model, line, label = label, color = 'k')
    # plt.plot(x_model, upper_band, color = 'w')
    # plt.plot(x_model, lower_band, color = 'w')
    plt.fill_between(x_model, upper_band, lower_band, color = 'k', alpha = 0.1)
    plt.xlim((np.min(x) - 0.5*np.min(x),np.max(x) + 0.5*np.max(x)))
    plt.ylim((0,np.max(y) + 0.5*np.max(y)))
    plt.legend(loc =2)

    # change the following 3 according to data input
    # plt.xlabel("Brightness Temperature (K km s$^{-1}$)")
    # plt.ylabel("$\gamma$-ray Excess Counts")
    # plt.savefig("correlation_study_data/plots/g_m_reg" + str(i+1)+ ".pdf")
    if correlation_type == 1:
        plt.xlabel("$N_p$ Column Density (Mopra + SGPS) (x10$^{21}$ cm$^{-2}$)")
        plt.ylabel("$\gamma$-ray Excess Counts")
        plt.title("Region %d" %(i+1))
        plt.savefig("correlation_study_data/plots/g_mHI_reg" + str(i+1)+ ".pdf")
    elif correlation_type == 2:
        plt.xlabel("$N_p$ Column Density (Nanten) (x10$^{21}$ cm$^{-2}$)")
        plt.ylabel("$\gamma$-ray Excess Counts")
        plt.title("Region %d" %(i+1))
        plt.savefig("correlation_study_data/plots/g_Np_reg" + str(i+1)+ ".pdf")
    elif correlation_type == 3:
        plt.plot(x_model,x_model, color = 'r')
        plt.xlabel("Nanten Brightness Temperature (K km s$^{-1}$)")
        plt.ylabel("Mopra Brightness Temperature (K km s$^{-1}$)")
        plt.xlim((0,2500))
        plt.ylim((0,2500))
        plt.savefig("correlation_study_data/plots/mop_nan_reg" + str(ifile_n) + ".pdf")
    elif correlation_type == 4:
        plt.xlabel("Suzaku X-ray Excess Counts")
        plt.ylabel("$N_p$ Column Density (Mopra) (x10$^{21}$ cm$^{-2}$)")
        plt.savefig("correlation_study_data/plots/mopHI_suz_reg3.pdf")

    # plt.xlabel("HI Column Density (x10$^{21}$)(cm$^{-2}$)")
    # plt.ylabel("$\gamma$-ray Excess Counts")
    # plt.savefig("correlation_study_data/plots/g_HI_reg" + str(i+1)+ ".pdf")


    plt.clf()

    plt.close()
