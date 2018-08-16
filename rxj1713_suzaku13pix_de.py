
# coding: utf-8

# In[94]:


import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import coordinates as coord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
from matplotlib.patches import Rectangle


# In[170]:


filename1 = '../DS9/suzaku_rebin.fits'
hdu1 = fits.open(filename1)
suzaku_data = hdu1[0].data

filename2 = '../DS9/rxj1713_Np_mopra_13pixels.fits'
hdu2 = fits.open(filename2)
Np_mop_data =hdu2[0].data


# In[165]:


suzaku_data.shape
Np_mop_data.shape


# In[51]:


reg3 = coord.SkyCoord(ra=258.4457*u.deg, dec=-39.2666*u.deg, frame='fk5')
reg3
reg3.dec.hms
reg3.ra.hms


# In[186]:


suzaku_data2 = suzaku_data[~np.isnan(suzaku_data)]
suzaku_data3 = [i for i in suzaku_data2 if i >= 2.0e-6]
Np_mop_data2 = Np_mop_data[~np.isnan(Np_mop_data)]
Np_mop_data3 = [i for i in Np_mop_data2 if i >= 1.0e21]

#use range when you have nan in your data
# plt.hist(suzaku_data2, bins=81)
# plt.show()
# plt.hist(suzaku_data3, bins=81)
# plt.show()
# plt.hist(Np_mop_data2, bins=81)
# plt.show()
# plt.hist(Np_mop_data3, bins=81)
# plt.show()
#
#
# # In[99]:
#
#
# fig = plt.figure(figsize=(9,9))
#
# wcs = WCS(suzaku_header)
# ax = WCSAxes(fig, [0,0,1,1], wcs=wcs)
# fig.add_axes(ax)
#
# plt.imshow(suzaku_data, cmap=cm.afmhot) #vlim is limiting the colourbar range
#
# cb = plt.colorbar()
# cb.set_label('Counts')
#
# # ax.scatter(258.4457, -39.2666, transform=ax.get_transform('fk5'), s=300, edgecolor='white', facecolor='none')
# r = Rectangle((258.4457, -39.2666), -0.18, -0.18, edgecolor='white', facecolor='none', transform=ax.get_transform('fk5'))
# ax.add_patch(r)
#
# plt.xlabel('RA')
# plt.ylabel('Dec')

#plt.tight_layout()
# plt.savefig('test.eps', bbox_inches='tight')


# In[188]:


s_data = []
Np_data = []
for i in range(29):
    s_data.append([])
    Np_data.append([])

s_data[0].append(suzaku_data[81:94,41:54].flatten())
s_data[1].append(suzaku_data[81:94,54:67].flatten())
s_data[2].append(suzaku_data[81:94,67:80].flatten())
s_data[3].append(suzaku_data[81:94,80:93].flatten())

s_data[4].append(suzaku_data[68:81,28:41].flatten())
s_data[5].append(suzaku_data[68:81,41:54].flatten())
s_data[6].append(suzaku_data[68:81,54:67].flatten())
s_data[7].append(suzaku_data[68:81,67:80].flatten())
s_data[8].append(suzaku_data[68:81,80:93].flatten())
s_data[9].append(suzaku_data[68:81,93:106].flatten())

s_data[10].append(suzaku_data[55:68,28:41].flatten())
s_data[11].append(suzaku_data[55:68,41:54].flatten())
s_data[12].append(suzaku_data[55:68,54:67].flatten())
s_data[13].append(suzaku_data[55:68,67:80].flatten())
s_data[14].append(suzaku_data[55:68,80:93].flatten())
s_data[15].append(suzaku_data[55:68,93:106].flatten())

s_data[16].append(suzaku_data[42:55,28:41].flatten())
s_data[17].append(suzaku_data[42:55,41:54].flatten())
s_data[18].append(suzaku_data[42:55,54:67].flatten())
s_data[19].append(suzaku_data[42:55,67:80].flatten())
s_data[20].append(suzaku_data[42:55,80:93].flatten())
s_data[21].append(suzaku_data[42:55,93:106].flatten())

s_data[22].append(suzaku_data[29:42,28:41].flatten())
s_data[23].append(suzaku_data[29:42,41:54].flatten())
s_data[24].append(suzaku_data[29:42,54:67].flatten())
s_data[25].append(suzaku_data[29:42,67:80].flatten())
s_data[26].append(suzaku_data[29:42,80:93].flatten())

s_data[27].append(suzaku_data[16:29,54:67].flatten())
s_data[28].append(suzaku_data[16:29,67:80].flatten())

Np_data[0].append(Np_mop_data[81:94,41:54].flatten())
Np_data[1].append(Np_mop_data[81:94,54:67].flatten())
Np_data[2].append(Np_mop_data[81:94,67:80].flatten())
Np_data[3].append(Np_mop_data[81:94,80:93].flatten())

Np_data[4].append(Np_mop_data[68:81,28:41].flatten())
Np_data[5].append(Np_mop_data[68:81,41:54].flatten())
Np_data[6].append(Np_mop_data[68:81,54:67].flatten())
Np_data[7].append(Np_mop_data[68:81,67:80].flatten())
Np_data[8].append(Np_mop_data[68:81,80:93].flatten())
Np_data[9].append(Np_mop_data[68:81,93:106].flatten())

Np_data[10].append(Np_mop_data[55:68,28:41].flatten())
Np_data[11].append(Np_mop_data[55:68,41:54].flatten())
Np_data[12].append(Np_mop_data[55:68,54:67].flatten())
Np_data[13].append(Np_mop_data[55:68,67:80].flatten())
Np_data[14].append(Np_mop_data[55:68,80:93].flatten())
Np_data[15].append(Np_mop_data[55:68,93:106].flatten())

Np_data[16].append(Np_mop_data[42:55,28:41].flatten())
Np_data[17].append(Np_mop_data[42:55,41:54].flatten())
Np_data[18].append(Np_mop_data[42:55,54:67].flatten())
Np_data[19].append(Np_mop_data[42:55,67:80].flatten())
Np_data[20].append(Np_mop_data[42:55,80:93].flatten())
Np_data[21].append(Np_mop_data[42:55,93:106].flatten())

Np_data[22].append(Np_mop_data[29:42,28:41].flatten())
Np_data[23].append(Np_mop_data[29:42,41:54].flatten())
Np_data[24].append(Np_mop_data[29:42,54:67].flatten())
Np_data[25].append(Np_mop_data[29:42,67:80].flatten())
Np_data[26].append(Np_mop_data[29:42,80:93].flatten())

Np_data[27].append(Np_mop_data[16:29,54:67].flatten())
Np_data[28].append(Np_mop_data[16:29,67:80].flatten())


# # Remove x-ray points with nan and corresponding Np data points

# In[203]:


# for i in range(29):
#     for j in range(len(s_data[0])):
#         if


# # Remove noisey X-ray data and corresponding gas data

# In[154]:


# s_data_reg1NR = [i for i in s_data_reg1 if i >= 2.0e-6]
# s_data_reg2NR = [i for i in s_data_reg2 if i >= 2.0e-6]
# s_data_reg3NR = [i for i in s_data_reg3 if i >= 2.0e-6]
# s_data_reg4NR = [i for i in s_data_reg4 if i >= 2.0e-6]
# s_data_reg5NR = [i for i in s_data_reg5 if i >= 2.0e-6]
# s_data_reg6NR = [i for i in s_data_reg6 if i >= 2.0e-6]
# s_data_reg7NR = [i for i in s_data_reg7 if i >= 2.0e-6]
# s_data_reg8NR = [i for i in s_data_reg8 if i >= 2.0e-6]
# s_data_reg9NR = [i for i in s_data_reg9 if i >= 2.0e-6]
# s_data_reg10NR = [i for i in s_data_reg10 if i >= 2.0e-6]
# s_data_reg11NR = [i for i in s_data_reg11 if i >= 2.0e-6]
# s_data_reg12NR = [i for i in s_data_reg12 if i >= 2.0e-6]
# s_data_reg13NR = [i for i in s_data_reg13 if i >= 2.0e-6]
# s_data_reg14NR = [i for i in s_data_reg14 if i >= 2.0e-6]
# s_data_reg15NR = [i for i in s_data_reg15 if i >= 2.0e-6]
# s_data_reg16NR = [i for i in s_data_reg16 if i >= 2.0e-6]
# s_data_reg17NR = [i for i in s_data_reg17 if i >= 2.0e-6]
# s_data_reg18NR = [i for i in s_data_reg18 if i >= 2.0e-6]
# s_data_reg19NR = [i for i in s_data_reg19 if i >= 2.0e-6]
# s_data_reg20NR = [i for i in s_data_reg20 if i >= 2.0e-6]
# s_data_reg21NR = [i for i in s_data_reg21 if i >= 2.0e-6]
# s_data_reg22NR = [i for i in s_data_reg22 if i >= 2.0e-6]
# s_data_reg23NR = [i for i in s_data_reg23 if i >= 2.0e-6]
# s_data_reg24NR = [i for i in s_data_reg24 if i >= 2.0e-6]
# s_data_reg25NR = [i for i in s_data_reg25 if i >= 2.0e-6]
# s_data_reg26NR = [i for i in s_data_reg26 if i >= 2.0e-6]
# s_data_reg27NR = [i for i in s_data_reg27 if i >= 2.0e-6]
# s_data_reg28NR = [i for i in s_data_reg28 if i >= 2.0e-6]
# s_data_reg29NR = [i for i in s_data_reg29 if i >= 2.0e-6]

#Np_mopra_HINR = [i for i in da]
