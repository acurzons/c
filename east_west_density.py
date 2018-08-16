import numpy as np

filename = "rxj1713_regional_densities.txt"
data = np.loadtxt(filename)

nH = data[:,2]

nH_east = (nH[0]+nH[1]+nH[4]+nH[5]+nH[6]+nH[10]+nH[11]+nH[12]+nH[16]+nH[17]+nH[18]+nH[22]+nH[23]+nH[24]+nH[27]) / 15.0

nH_west = (nH[2]+nH[3]+nH[7]+nH[8]+nH[9]+nH[13]+nH[14]+nH[15]+nH[19]+nH[20]+nH[21]+nH[25]+nH[26]+nH[28]) / 14.0

print "East:         West:"
print nH_east, nH_west
