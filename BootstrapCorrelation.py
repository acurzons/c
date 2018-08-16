
# coding: utf-8

# In[235]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.optimize import curve_fit
from scipy import stats


# # Uncertainty of Correlation Coefficient
#
# ... due to uncertainty of data points.

# In[184]:


def korr(x,y,xerr,yerr,plot=False):

    xx = np.random.normal(x,xerr)
    yy = np.random.normal(y,yerr)

    if plot==True:
        plt.errorbar(xx, yy, xerr=xerr, yerr=yerr, fmt='o')
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color='y')
        plt.show()
        print(pearsonr(xx,yy))

    return pearsonr(xx,yy)[0]


# In[191]:


# Put your data here
xray = np.random.normal(5,2,20)
gas = np.random.normal(15,5,20)

xray_err = np.random.normal(0.5,0.1,20)
gas_err = np.random.normal(2,0.5,20)


# In[192]:


korr(xray, gas, xray_err, gas_err, True)


# In[198]:


k = []
for i in range(100000):
    k.append(korr(xray,gas,xray_err,gas_err))

plt.hist(k, bins=81, range=[-1,1])
# 3 sigma
plt.axvline(np.percentile(k, (100.0 - 99.73)/2.0), color='black', ls="--",)
plt.axvline(np.percentile(k, (100.0 - (100.0 - 99.73)/2.0)), color='black', ls="--")

# 5 sigma
plt.axvline(np.percentile(k, (100.0 - 99.99994)/2.0), color='black', ls=":")
plt.axvline(np.percentile(k, (100.0 - (100.0 - 99.99994)/2.0)), color='black', ls=":")

# correlation of original data set
plt.axvline(pearsonr(xray, gas)[0], color='black')

plt.show()


# --> One could use (correlation coefficient + uncertainty) as some extreme value like "in the most optimistic case (uncertainty of Z sigma),
#       the most extreme correlation coefficient is X." And then use this extreme value in the following to determine the corresponding p-value.

# # p-value

# In[200]:


def korr2(x_mean, x_sigm, y_mean, y_sigm, number):

    x = np.random.normal(x_mean,x_sigm,number)
    y = np.random.normal(y_mean,y_sigm,number)

    return pearsonr(x,y)[0]


# In[202]:


# Make some assumptions about your data
# For example that they are normal distributed with specific values
# Specifiy the values

x_mean = np.mean(xray)
x_sigm = np.std(xray)
y_mean = np.mean(gas)
y_sigm = np.std(gas)
number = len(xray)

print(x_mean, x_sigm, y_mean, y_sigm)


# --> Decide a confidence level for your test if the data is correlated and determine the corresponding correlation coefficient.
#     Then check if your optimistic value is larger (then it is correlated wrt to the CL) or smaller. You can also quote the p-value then.

# In[232]:


kk = []
for i in range(100000):
    kk.append(korr2(x_mean, x_sigm, y_mean, y_sigm, number))

plt.hist(np.abs(kk), bins=81, range=[0,1])

print('2 sigma', np.percentile(np.abs(kk),95.449))
print('3 sigma', np.percentile(np.abs(kk),99.73))
print('5 sigma', np.percentile(np.abs(kk),99.99994))

plt.show()


# In[234]:


# A check of the method
# Compare if bootstrap result is consistent with Pearson result
print(np.abs(pearsonr(xray,gas)[0]))
print(np.percentile(np.abs(kk),100-pearsonr(xray,gas)[1]*100.))
