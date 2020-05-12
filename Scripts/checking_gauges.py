"""
Quick checking script to plot gauge history vs time
"""
import preamble_functions as pre
import matplotlib.pyplot as plt #3.0.2
from matplotlib.lines import Line2D
import os
import numpy as np
from scipy import ndimage
from impulse_models import *

params = {'font.family':'serif',
        'axes.labelsize':'small',
        'xtick.labelsize':'x-small',
        'ytick.labelsize':'x-small', 
        'legend.fontsize':'small',
        'legend.title_fontsize':'small',
        'grid.linestyle':'--',
        'grid.linewidth':'0.5',
        'lines.linewidth':'0.5'}
plt.rcParams.update(params)


# main_dataset_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_dataset\*.txt")
# original_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*.txt")
# original_file = original_file[4::]
# latest_1500_r3_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res3\*.txt")
# latest_1500_r4_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res4\*.txt")
# latest_1500_r5_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res5\*.txt")
# latest_var_r4_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\VAR_ZL40mm_res4\*.txt")



# main_dataset = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_dataset\*gtable",1)
# original = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*gtable",1)
# original = original[4::]
# latest_1500_r3 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res3\*gtable",1)
# latest_1500_r4 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res4\*gtable",1)
# latest_1500_r5 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res5\*gtable",1)
# latest_var_r4 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\VAR_ZL40mm_res4\*gtable",1)


# main_dataset = np.asarray([main_dataset[i][:,7] for i in range(len(main_dataset))]).T
# original = np.asarray([original[i][:,7] for i in range(len(original))]).T
# latest_1500_r3 = np.asarray([latest_1500_r3[i][:,7] for i in range(len(latest_1500_r3))]).T
# latest_1500_r4 = np.asarray([latest_1500_r4[i][:,7] for i in range(len(latest_1500_r4))]).T
# latest_1500_r5 = np.asarray([latest_1500_r5[i][:,7] for i in range(len(latest_1500_r5))]).T
# latest_var_r4 = np.asarray([latest_var_r4[i][:,7] for i in range(len(latest_var_r4))]).T




# fig, [ax0, ax1] = plt.subplots(2,1)
# ax0.plot(np.linspace(0,80,300), original, 'k')
# ax0.plot(np.linspace(0,80,200), main_dataset, 'g')
# ax0.plot(np.linspace(0,80,200), latest_1500_r3, 'b')
# ax0.plot(np.linspace(0,80,200), latest_1500_r4, 'y')
# ax0.plot(np.linspace(0,80,200), latest_1500_r5,  'm')
# ax0.plot(np.linspace(0,80,200), latest_var_r4, 'c')
# ax0.set_xlabel('Theta')
# ax1.set_ylabel('Peak specific impulse')

# #ax1.plot(np.linspace(0,80,200), main_dataset/main_dataset.max(0), 'g')
# ax1.plot(np.linspace(0,80,200), latest_1500_r3/latest_1500_r3.max(0), 'b')
# ax1.plot(np.linspace(0,80,200), latest_1500_r4/latest_1500_r4.max(0), 'y')
# ax1.plot(np.linspace(0,80,200), latest_1500_r5/latest_1500_r5.max(0),  'm')
# ax1.plot(np.linspace(0,80,200), latest_var_r4/latest_var_r4.max(0), 'c')
# ax1.set_xlabel('Theta')
# ax1.set_ylabel('Peak specific impulse ratio')
# labels = ['1500mm^3 domain, 100mm ZL, res3', '1500mm^3 domain, 100mm ZL, res4', '1500mm^3 domain, 100mm ZL, res5',
#           'variable domain, 40mm ZL, res4',]
# colors = ['b', 'y', 'm', 'c']
# lines = [Line2D([0], [0], color=c, linewidth=0.5) for c in colors]
# ax1.legend(lines,labels)


# main_dataset_I = [Impulse_CFD(main_dataset[:,i], 1, 80, np.linspace(0,80,200)) for i in range(len(main_dataset_file))]
# original_I = [Impulse_CFD(original[:,i], 1, 80, np.linspace(0,80,300)) for i in range(len(original_file))]
# latest_1500_r3_I = [Impulse_CFD(latest_1500_r3[:,i], 1, 80, np.linspace(0,80,200)) for i in range(len(latest_1500_r3_file))]
# latest_1500_r4_I = [Impulse_CFD(latest_1500_r4[:,i], 1, 80, np.linspace(0,80,200)) for i in range(len(latest_1500_r4_file))]
# latest_1500_r5_I = [Impulse_CFD(latest_1500_r5[:,i], 1, 80, np.linspace(0,80,200)) for i in range(len(latest_1500_r5_file))]
# latest_var_r4_I = [Impulse_CFD(latest_var_r4[:,i], 1, 80, np.linspace(0,80,200)) for i in range(len(latest_var_r4_file))]



#--------------------------------------------------------------------------------
main_dataset_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_055\*.txt")
main_dataset = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_055\*_gauges",1)
Apollo_gauges = [main_dataset[i] for i in [0,1,4,5]]

#Quick test to see if enough time in simulation

for i in range(len(Apollo_gauges)):
    fig00, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
    ax1.set_xlabel(i)
    #theta = 0
    ax1.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,1]) #OP
    ax2.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,201])#imp
    
    #theta = 80
    ax3.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,200])#OP
    ax4.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,400])#imp

# lim = 0.6
# fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
# ax1.plot(Apollo_gauges[0][:,0]*1000,Apollo_gauges[0][:,1]) #OP
# ax2.plot(Apollo_gauges[0][:,0]*1000,Apollo_gauges[0][:,201])#imp
# ax1.set_xlim(0,lim)
# ax2.set_xlim(0,lim)   
# #theta = 80
# ax3.plot(Apollo_gauges[0][:,0]*1000,Apollo_gauges[0][:,200])#OP
# ax4.plot(Apollo_gauges[0][:,0]*1000,Apollo_gauges[0][:,385])#imp
# ax3.set_xlim(0,lim)
# ax4.set_xlim(0,lim)




