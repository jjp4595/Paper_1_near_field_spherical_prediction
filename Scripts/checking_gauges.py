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
        'scatter.marker': 's',
        'lines.linewidth':'0.5'}
plt.rcParams.update(params)


#Some charge properties
charge_rad = 0.0246
charge_mass = 0.1



main_dataset_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_dataset\*.txt")
original_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*.txt")
original_file = original_file[4::]
latest_1500_r3_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res3\*.txt")
latest_1500_r4_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res4\*.txt")
latest_1500_r5_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res5\*.txt")
latest_var_r4_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\VAR_ZL40mm_res4\*.txt")
testing_DMA_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\testing_DMA\*.txt")

#clear Z
main_dataset_z = [(pre.standoff_func(main_dataset_file[i]) - charge_rad)/(charge_mass**(1/3)) for i in range(len(main_dataset_file))]
original_z = [(pre.standoff_func(original_file[i]) - charge_rad)/(charge_mass**(1/3)) for i in range(len(original_file))]

latest_1500_r3_z = [(pre.standoff_func(latest_1500_r3_file[i]) - charge_rad)/(charge_mass**(1/3)) for i in range(len(latest_1500_r3_file))]
latest_1500_r4_z = [(pre.standoff_func(latest_1500_r4_file[i]) - charge_rad)/(charge_mass**(1/3)) for i in range(len(latest_1500_r4_file))]
latest_1500_r5_z = [(pre.standoff_func(latest_1500_r5_file[i]) - charge_rad)/(charge_mass**(1/3)) for i in range(len(latest_1500_r5_file))]
latest_var_r4_z = [(pre.standoff_func(latest_var_r4_file[i]) - charge_rad)/(charge_mass**(1/3)) for i in range(len(latest_var_r4_file))]
testing_DMA_z = [(pre.standoff_func(testing_DMA_file[i]) - charge_rad)/(charge_mass**(1/3)) for i in range(len(testing_DMA_file))]

#Gtables
main_dataset = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_dataset\*gtable",1)
original = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*gtable",1)
original = original[4::]

latest_1500_r3 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res3\*gtable",1)
latest_1500_r4 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res4\*gtable",1)
latest_1500_r5 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res5\*gtable",1)
latest_var_r4 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\VAR_ZL40mm_res4\*gtable",1)
testing_DMA = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\testing_DMA\*gtable", 1)
main_dataset = np.asarray([main_dataset[i][:,7] for i in range(len(main_dataset))]).T
original = np.asarray([original[i][:,7] for i in range(len(original))]).T
latest_1500_r3 = np.asarray([latest_1500_r3[i][:,7] for i in range(len(latest_1500_r3))]).T
latest_1500_r4 = np.asarray([latest_1500_r4[i][:,7] for i in range(len(latest_1500_r4))]).T
latest_1500_r5 = np.asarray([latest_1500_r5[i][:,7] for i in range(len(latest_1500_r5))]).T
latest_var_r4 = np.asarray([latest_var_r4[i][:,7] for i in range(len(latest_var_r4))]).T
testing_DMA = np.asarray([testing_DMA[i][:,7] for i in range(len(testing_DMA))]).T

#Gauges
latest_1500_r3_gauges = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res3\*gauges",1)
latest_1500_r4_gauges = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res4\*gauges",1)
latest_1500_r5_gauges = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res5\*gauges",1)
latest_var_r4_gauges = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\VAR_ZL40mm_res4\*gauges",1)
testing_DMA_gauges = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\testing_DMA\*gauges", 1)






def term_time_diff():
    #peakIs from gauges for a specified term time
    time = 2e-3
    latest_1500_r3_i =  latest_1500_r3_gauges[3][latest_1500_r3_gauges[3][:,0]<time,201::].max(axis=0)
    latest_1500_r4_i =  latest_1500_r4_gauges[3][latest_1500_r4_gauges[3][:,0]<time,201::].max(axis=0)
    latest_1500_r5_i =  latest_1500_r5_gauges[3][latest_1500_r5_gauges[3][:,0]<time,201::].max(axis=0)
    latest_var_r4_i = latest_var_r4_gauges[3][latest_var_r4_gauges[3][:,0]<time,201::].max(axis=0)
    #plot testing if term time of 2ms or 3ms makes a difference, result - not really.
    fig, ax0 = plt.subplots(1,1)
    ax0.set_title("Testing termination time sensitivity")
    ax0.plot(np.linspace(0,80,200),latest_var_r4[:,3], 'c')
    ax0.plot(np.linspace(0,80,200),latest_1500_r3[:,3], 'g')
    ax0.plot(np.linspace(0,80,200),latest_1500_r4[:,3], 'r')
    ax0.plot(np.linspace(0,80,200),latest_1500_r5[:,3], 'm')
    ax0.plot(np.linspace(0,80,200),latest_var_r4_i, 'c:')
    ax0.plot(np.linspace(0,80,200),latest_1500_r3_i, 'g:')
    ax0.plot(np.linspace(0,80,200),latest_1500_r4_i, 'r:')
    ax0.plot(np.linspace(0,80,200),latest_1500_r5_i, 'm:')
term_time_diff()


def dataset_overview():
#Overview of datasets
    fig, [[ax0, ax1, ax2],[ax3, ax4, ax5]] = plt.subplots(2,3)
    ax0.plot(np.linspace(0,80,200),main_dataset)
    ax0.set_title("main dataset - 0.02m, res3 - 2.5mm")
    ax1.plot(np.linspace(0,80,300),original)
    ax1.set_title("original - 0.1m, res4 - 6.25mm")
    ax2.plot(np.linspace(0,80,200), latest_var_r4)
    ax2.set_title("latest domain = 1.2xgaugelength - 0.04m, res4 - 2.5mm")
    ax3.plot(np.linspace(0,80,200),latest_1500_r3)
    ax3.set_title("latest_1500mm - 0.1m, res3 - 12.5mm")
    ax4.plot(np.linspace(0,80,200),latest_1500_r4)
    ax4.set_title("latest_1500mm - 0.1m, res4 - 6.25mm")
    ax5.plot(np.linspace(0,80,200),latest_1500_r5)
    ax5.set_title("latest_1500mm - 0.1m, re5 - 3.125mm")
dataset_overview()


#Looking at 
fig, [ax, ax0, ax1] = plt.subplots(1,3)
fig.set_size_inches(7,2.5)
#ax.set_title('Dataset, 5 samples:    0.055 < Z(clear s/o) < 0.160')
ax.plot(np.linspace(0,80,200),latest_var_r4/1e3, 'c')
ax.plot(np.linspace(0,80,200),latest_1500_r3/1e3, 'g')
ax.plot(np.linspace(0,80,200),latest_1500_r4/1e3, 'r')
ax.plot(np.linspace(0,80,200),latest_1500_r5/1e3, 'm')
ax.plot(np.linspace(0,80,200),testing_DMA[:,0]/1e3, 'k')
ax.plot(np.linspace(0,80,200),testing_DMA[:,1]/1e3, 'b')
ax.set_ylabel('peak specific impulse (MPa.ms)')
ax.set_xlabel('angle of incidence')
labels = ['A', 'B', 'C', 'D', 'E', 'F']
colors = ['c', 'g', 'r', 'm', 'k', 'b']
lines = [Line2D([0], [0], color=c, linewidth=0.5) for c in colors]
ax.legend(lines,labels)
#ax0.set_title('Sample of dataset, z(clear) = 0.133')
# ax0.plot(np.linspace(0,80,200),latest_var_r4[:,3]/1e3, 'c', label = 'A')
# ax0.plot(np.linspace(0,80,200),latest_1500_r3[:,3]/1e3, 'g', label = 'B')
# ax0.plot(np.linspace(0,80,200),latest_1500_r4[:,3]/1e3, 'r', label = 'C')
# ax0.plot(np.linspace(0,80,200),latest_1500_r5[:,3]/1e3, 'm', label = 'D')
ax0.plot(np.linspace(0,80,200),testing_DMA[:,0]/1e3, 'k', label = 'E')
ax0.plot(np.linspace(0,80,200),testing_DMA[:,1]/1e3, 'b',label = 'F')
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='upper right', prop={'size':6})
ax0.set_xlabel('angle of incidence')

ax1.scatter(latest_var_r4_z, np.asarray(latest_var_r4_I)/1e3, c='c', s=10., label = 'A')
ax1.scatter(latest_1500_r3_z, np.asarray(latest_1500_r3_I)/1e3, c='g', s=10., label = 'B')
ax1.scatter(latest_1500_r4_z, np.asarray(latest_1500_r4_I)/1e3, c='r', s=10., label = 'C')
ax1.scatter(latest_1500_r5_z, np.asarray(latest_1500_r5_I)/1e3, c='m', s=10., label = 'D')
ax1.scatter(testing_DMA_z[0], np.asarray(testing_DMA_I[0])/1e3, c='k', s=10., label = 'E')
ax1.scatter(testing_DMA_z[1], np.asarray(testing_DMA_I[1])/1e3, c='b', s=10., label = 'F')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='upper right', prop={'size':6})
ax1.set_xlabel('Z (clear standoff)')
ax1.set_ylabel('Total Impulse (MPa.ms)')
#ax1.set_title(r'$1m^2$ area integrated impulse')
plt.tight_layout()
fig.savefig('dma.pdf', format='pdf')



fig, [ax0, ax1] = plt.subplots(2,1)
ax0.plot(np.linspace(0,80,300), original, 'k')
ax0.plot(np.linspace(0,80,200), main_dataset, 'g')
ax0.plot(np.linspace(0,80,200), latest_1500_r3, 'b') 
ax0.plot(np.linspace(0,80,200), latest_1500_r4, 'y')
ax0.plot(np.linspace(0,80,200), latest_1500_r5,  'm')
ax0.plot(np.linspace(0,80,200), latest_var_r4, 'c')
ax0.set_xlabel('Theta')
ax1.set_ylabel('Peak specific impulse')
ax1.plot(np.linspace(0,80,200), main_dataset/main_dataset.max(0), 'g')
ax1.plot(np.linspace(0,80,200), latest_1500_r3/latest_1500_r3.max(0), 'b')
ax1.plot(np.linspace(0,80,200), latest_1500_r4/latest_1500_r4.max(0), 'y')
ax1.plot(np.linspace(0,80,200), latest_1500_r5/latest_1500_r5.max(0),  'm')
ax1.plot(np.linspace(0,80,200), latest_var_r4/latest_var_r4.max(0), 'c')
ax1.set_xlabel('Theta')
ax1.set_ylabel('Peak specific impulse ratio')
labels = ['1500mm^3 domain, 100mm ZL, res3', '1500mm^3 domain, 100mm ZL, res4', '1500mm^3 domain, 100mm ZL, res5',
          'variable domain, 40mm ZL, res4',]
colors = ['b', 'y', 'm', 'c']
lines = [Line2D([0], [0], color=c, linewidth=0.5) for c in colors]
ax1.legend(lines,labels)


main_dataset_I = [Impulse_CFD(main_dataset[:,i], 1, 80, np.linspace(0,80,200)) for i in range(len(main_dataset_file))]
original_I = [Impulse_CFD(original[:,i], 1, 80, np.linspace(0,80,300)) for i in range(len(original_file))]
latest_1500_r3_I = [Impulse_CFD(latest_1500_r3[:,i], 1, 80, np.linspace(0,80,200)) for i in range(len(latest_1500_r3_file))]
latest_1500_r4_I = [Impulse_CFD(latest_1500_r4[:,i], 1, 80, np.linspace(0,80,200)) for i in range(len(latest_1500_r4_file))]
latest_1500_r5_I = [Impulse_CFD(latest_1500_r5[:,i], 1, 80, np.linspace(0,80,200)) for i in range(len(latest_1500_r5_file))]
latest_var_r4_I = [Impulse_CFD(latest_var_r4[:,i], 1, 80, np.linspace(0,80,200)) for i in range(len(latest_var_r4_file))]
testing_DMA_I = [Impulse_CFD(testing_DMA[:,i], 1, 80, np.linspace(0,80,200)) for i in range(len(testing_DMA_file))]






# #Quick test to see if enough time in simulation
# for i in range(len(Apollo_gauges)):
#     fig00, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
#     ax1.set_xlabel(i)
#     #theta = 0
#     ax1.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,1]) #OP
#     ax2.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,201])#imp
    
#     #theta = 80
#     ax3.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,200])#OP
#     ax4.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,400])#imp



# #-----------------------
# val_80mm_chosenmesh_gauges = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\80mm_validation\*_gauges",1)
# val_80mm_chosenmesh = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\80mm_validation\*_gtable",1)
# Apollo_gauges = val_80mm_chosenmesh_gauges

# Apollo_gtable_z80mm_first = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Near Field Sims\Sims\Latest\80mm_with_afterburn\*gtable",1)
# z80mm_first_theta = np.rad2deg(np.arctan2(Apollo_gtable_z80mm_first[0][:,2], 0.08))

# #Quick test to see if enough time in simulation
# for i in range(len(Apollo_gauges)):
#     fig00, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
#     ax1.set_xlabel(i)
#     #theta = 0
#     ax1.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,1]) #OP
#     ax2.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,201])#imp
    
#     #theta = 80
#     ax3.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,200])#OP
#     ax4.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,400])#imp
    
    
# #peak impulse distribution for first 0.3ms
# term = 0.3e-3
# z80mm_chosenmesh_fin = int(np.argwhere(Apollo_gauges[0][:,0]>term)[0][0])
# fig8, [ax, ax0] = plt.subplots(1,2)
# fig8.set_size_inches(5, 3)
# ax.set_xlabel('theta (degrees)')
# ax.set_ylabel('peak specific scaled impulse')
# ax.plot(np.linspace(0,80,200), np.max(Apollo_gauges[0][0:z80mm_chosenmesh_fin,201:], axis = 0), 'k', label = '3.125mm')

# ax0.set_xlabel('theta (degrees)')
# ax0.set_ylabel('peak impulse ratio of maximum')
# ax0.plot(theta, Apollo_gtable_z80mm_chosenmesh[0][:,7]/max(Apollo_gtable_z80mm_chosenmesh[0][:,7]), 'k', label = '3.125mm')
# ax0.plot(z80mm_first_theta, Apollo_gtable_z80mm_first[0][:,7]/max(Apollo_gtable_z80mm_first[0][:,7]), 'k', label = '3.125mm')


# handles, labels = ax0.get_legend_handles_labels()
# ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.7, 0.80), prop={'size':6})
# plt.tight_layout()




