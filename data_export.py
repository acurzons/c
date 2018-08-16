import numpy as np
from astropy.io import ascii
import os.path

N = 9
N2 = 169
no_reg = 29
xfactor = 2.0e20
odata1 = np.zeros((9,4))
odata2 = np.zeros((9,4))
odata3 = np.zeros((9,4))
odata4 = np.zeros((9,4))
odata5 = np.zeros((9,4))
odata6 = np.zeros((9*no_reg,4))
odata7_n = np.empty(0)
odata7_m = np.empty(0)
odata7_ne = np.empty(0)
odata7_me = np.empty(0)

for i in range(no_reg):
    odata8 = np.zeros((169,4))

    ifile = "correlation_study_data/rawdata/reg" + str(i + 1) + ".txt"
    data = np.loadtxt(ifile, skiprows = 1)


    gamma = data[:,0]   #excess counts
    gamma_err = gamma**0.5

    mopra = data[:,1]    #Brightness temp
    mopra_err = np.ones(9) * 1.5 / np.sqrt(51.0) * np.sqrt(20.0/0.11) #Trms / myregrid_pixels * integrate_velocity_channels

    nanten = data[:,4] * 1.0e21 #column density
    nanten_err = np.ones(9) * 0.3 * np.sqrt(20.0/0.65) * xfactor * 2.0 / np.sqrt(2.0)    #Trms * sqrt(velocity_window/delta_v) * X_factor * NH2_to_Np_conversion * regrid_factor_sum
    nanten_BT = nanten / xfactor / 2.     #brightness temp K km s^-1
    nanten_BT_err = nanten_err  / xfactor / 2.

    HI = data[:,3] * 1.0e21    #column density
    HI_err = np.ones(9) * 1.6 * np.sqrt(20.0/0.82) * 1.823e18 / np.sqrt(3.0)    #Trms * sqrt(velocity_window/delta_v) *  temp_to_Np_conversion * regrid_factor

    nanten_HI = data[:,2]     #column density
    nanten_HI_err = (nanten_err + HI_err) /1.e21

    # Check that the individual HI + nanten column density maps agree with the total column density maps from Fukui et al
    # for j in range(len(data[:,0])):
    #     if abs(nanten_HI[j] - (nanten[j] + HI[j])) > 0.02:
    #         print i+1, j+1, nanten_HI[j] - (nanten[j] + HI[j])

    mopra_CD = 2. * xfactor * mopra
    mopra_CD_err = 2. * xfactor * mopra_err
    mopra_HI = (mopra_CD + HI ) / 1.e21           #column density
    mopra_HI_err = (mopra_CD_err + HI_err) /1.e21



    odata1[:,0] = gamma
    odata1[:,1] = mopra
    odata1[:,2] = gamma_err
    odata1[:,3] = mopra_err
    odata2[:,0] = gamma
    odata2[:,1] = nanten_HI
    odata2[:,2] = gamma_err
    odata2[:,3] = nanten_HI_err
    odata3[:,0] = gamma
    odata3[:,1] = HI
    odata3[:,2] = gamma_err
    odata3[:,3] = HI_err
    odata4[:,0] = gamma
    odata4[:,1] = mopra_HI
    odata4[:,2] = gamma_err
    odata4[:,3] = mopra_HI_err
    odata5[:,0] = nanten_HI
    odata5[:,1] = mopra_HI
    odata5[:,2] = nanten_HI_err
    odata5[:,3] = mopra_HI_err
    odata6[i*9:9*(i+1),0] = nanten_BT
    odata6[i*9:9*(i+1),1] = mopra
    odata6[i*9:9*(i+1),2] = nanten_BT_err
    odata6[i*9:9*(i+1),3] = mopra_err

    odata7_n = np.hstack((odata7_n,nanten_BT))
    odata7_m = np.hstack((odata7_m,mopra))
    odata7_ne = np.hstack((odata7_ne,nanten_BT_err))
    odata7_me = np.hstack((odata7_me,mopra_err))

    if i == 28:
        n = 0
        for j in range(N*no_reg):
            if odata7_n[j] > 10:
                n += 1

    if i == 28:
        odata7 = np.zeros((n,4))
        k = 0
        for j in range(N*no_reg):
            if odata7_n[j] > 10:
                odata7[k,0] =  odata7_n[j]
                odata7[k,1] =  odata7_m[j]
                odata7[k,2] =  odata7_ne[j]
                odata7[k,3] =  odata7_me[j]
                k += 1


    ofile1 = "correlation_study_data/gamma_mopra_data/reg" + str(i + 1) + ".txt"
    ascii.write(odata1, ofile1, names = ["Gamma_Counts", "Mopra_Gas_Density(K km/s)", "Gamma error", "Mopra Error"], overwrite=True)
    ofile2 = "correlation_study_data/gamma_Np_data/reg" + str(i + 1) + ".txt"
    ascii.write(odata2, ofile2, names = ["Gamma_Counts", "Np Column Density (cm^-2)", "Gamma error", "Np Error"], overwrite=True)
    ofile3 = "correlation_study_data/gamma_HI_data/reg" + str(i + 1) + ".txt"
    ascii.write(odata3, ofile3, names = ["Gamma_Counts", "HI Column Density (cm^-2)", "Gamma error", "HI Error"], overwrite=True)
    ofile4 = "correlation_study_data/gamma_mopraHI_data/reg" + str(i + 1) + ".txt"
    ascii.write(odata4, ofile4, names = ["Gamma_Counts", "Np Column Density (cm^-2)", "Gamma error", "Np Error"], overwrite=True)
    ofile5 = "correlation_study_data/nantenHI_mopraHI_data/reg" + str(i + 1) + ".txt"
    ascii.write(odata5, ofile5, names = ["Np Column Density (Nanten)", "Np Column Density (Mopra)", "Nanten error", "Mopra Error"], overwrite=True)
    if i == 28:
        ofile6 = "correlation_study_data/nanten_mopra_data/all_data.txt"
        ascii.write(odata6, ofile6, names = ["Nanten Brightness Temp", "Mopra Brightness Temp", "Nanten error", "Mopra Error"], overwrite=True)
        ofile7 = "correlation_study_data/nanten_mopra_data/all_data_NR.txt"
        ascii.write(odata7, ofile7, names = ["Nanten Brightness Temp", "Mopra Brightness Temp", "Nanten error", "Mopra Error"], overwrite=True)



    ifile2 = "correlation_study_data/rawdata/reg" + str(i + 1) + "_13pixels.txt"
    if os.path.isfile(ifile2) == False:
        continue
    data2 = np.loadtxt(ifile2, skiprows = 1)
    print i+1

    suzaku = data2[:,0] #units of (1x10^-4)
    suzaku_err = np.abs(0.2*suzaku / np.sqrt(9.0))

    Np_mop_13 = data2[:,1] * 2.8
    Np_mop_13_err = np.ones(N2) * 2. * xfactor * np.sqrt(2.8) / np.sqrt(100) * 0.7 * np.sqrt(20.0/0.11) + np.ones(N2) * 1.9 * np.sqrt(20.0/0.82) * 1.823e18 * np.sqrt(2.8)

    odata8[:,0] = Np_mop_13
    odata8[:,1] = suzaku
    odata8[:,2] = Np_mop_13_err
    odata8[:,3] = suzaku_err

    #remove rows with nan
    ii = 0
    for k in range(len(odata8[:,0])):
        for j in range(len(odata8[0,:])):
            if np.isfinite(odata8[ii,j]) == False:
                odata8 = np.delete(odata8, ii, axis=0)
                ii = ii - 1
                break
        ii += 1

    odata9 = odata8

    #remove rows with negative values (noise)
    ii = 0
    for k in range(len(odata9[:,0])):
        for j in range(len(odata9[0,:])):
            if odata9[ii,j] < 0.0:
                odata9 = np.delete(odata9, ii, axis=0)
                ii = ii - 1
                break
        ii += 1

    ofile8 = "correlation_study_data/suzaku_mopHI_data/reg" + str(i + 1) + ".txt"
    ascii.write(odata8, ofile8, names = ["Np_Col_Dens(cm^-2)", "Suzaku_Xrays(excess_counts)", "Np error", "Suzaku Error"], overwrite=True)
    ofile9 = "correlation_study_data/suzaku_mopHI_data/NRreg" + str(i + 1) + ".txt"
    ascii.write(odata9, ofile9, names = ["Np_Col_Dens(cm^-2)", "Suzaku_Xrays(excess_counts)", "Np error", "Suzaku Error"], overwrite=True)
