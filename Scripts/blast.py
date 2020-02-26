"""
Conwep Eq's

"""
import scipy.io as sio #1.1.0
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt #3.0.2

conwep_import = sio.loadmat(r"C:\Users\jorda\Google Drive\Generic Matlab Scripts\Blast.m\ConWepValues")
air = conwep_import['air']
surface = conwep_import['surface']


def positive_parameters(R,W,bursttype,air,surface):
    conwep_import = sio.loadmat(r"C:\Users\jorda\Google Drive\Generic Matlab Scripts\Blast.m\ConWepValues")
    air = conwep_import['air']
    surface = conwep_import['surface']
    """
    W is Charge mass (kg TNT)
    R is stand-off (m)
    burstype: 1 for hemi-spherical surface burst, else spherical free-air burst
    """
    Wroot = W ** (1/3)
    Z = R / Wroot
    
    if bursttype == 1:
        psomax = np.power(10, griddata(np.log10(surface[:,0]), np.log10(surface[:,1]), np.log10(Z), method='cubic'))
        prmax = np.power(10, griddata(np.log10(surface[:,0]), np.log10(surface[:,2]), np.log10(Z), method='cubic'))
        ta = np.power(10, griddata(np.log10(surface[:,0]), np.log10(surface[:,3]), np.log10(Z), method='cubic'))
        td = np.power(10, griddata(np.log10(surface[:,0]), np.log10(surface[:,4]), np.log10(Z), method='cubic'))
        iso = np.power(10, griddata(np.log10(surface[:,0]), np.log10(surface[:,5]), np.log10(Z), method='cubic'))
        ir = np.power(10, griddata(np.log10(surface[:,0]), np.log10(surface[:,6]), np.log10(Z), method='cubic'))    
    else:
        psomax = np.power(10, griddata(np.log10(air[:,0]), np.log10(air[:,1]), np.log10(Z), method='cubic'))
        prmax = np.power(10, griddata(np.log10(air[:,0]), np.log10(air[:,2]), np.log10(Z), method='cubic'))
        ta = np.power(10, griddata(np.log10(air[:,0]), np.log10(air[:,3]), np.log10(Z), method='cubic'))
        td = np.power(10, griddata(np.log10(air[:,0]), np.log10(air[:,4]), np.log10(Z), method='cubic'))
        iso = np.power(10, griddata(np.log10(air[:,0]), np.log10(air[:,5]), np.log10(Z), method='cubic'))
        ir = np.power(10, griddata(np.log10(air[:,0]), np.log10(air[:,6]), np.log10(Z), method='cubic'))
    return psomax, prmax, ta, td, iso, ir


# ##Generic positive phase parameter plots --------------------------------------------------------------------------
# fig = plt.figure()
# fig.set_size_inches(6, 4)  
# ax = plt.subplot(1,2,1)
# ax.plot(surface[:,0], surface[:,1], label = 'pso, max')
# ax.plot(surface[:,0], surface[:,2], label = 'pr, max')
# ax.plot(surface[:,0], surface[:,3], label = 'ta')
# ax.plot(surface[:,0], surface[:,4], label = 'td')
# ax.plot(surface[:,0], surface[:,5], label = 'iso')
# ax.plot(surface[:,0], surface[:,6], label = 'ir')
# ax.set_title('Surface')
# ax.set_xlabel('scaled distance, Z')
# ax.set_yscale('log')
# ax.set_xscale('log')
# plt.gca().set_xlim(0.6, 50)
# plt.gca().set_ylim(0.01, 100000)
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels, loc='center', bbox_to_anchor=(0.7, 0.8), prop={'size':8})
# plt.tight_layout()

# ax = plt.subplot(1,2,2)
# ax.plot(air[:,0], air[:,1], label = 'pso, max')
# ax.plot(air[:,0], air[:,2], label = 'pr, max')
# ax.plot(air[:,0], air[:,3], label = 'ta')
# ax.plot(air[:,0], air[:,4], label = 'td')
# ax.plot(air[:,0], air[:,5], label = 'iso')
# ax.plot(air[:,0], air[:,6], label = 'ir')
# ax.set_title('Air')
# ax.set_xlabel('scaled distance, Z')
# ax.set_yscale('log')
# ax.set_xscale('log')
# plt.gca().set_xlim(0.6, 50)
# plt.gca().set_ylim(0.01, 100000)
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels, loc='center', bbox_to_anchor=(0.7, 0.8), prop={'size':8})
# plt.tight_layout()
# #--------------------------------------------------------------------------------------------------------------------


# #Input
# psomax, prmax, ta, td, iso, ir = positive_parameters(7,0.5,1,air,surface)