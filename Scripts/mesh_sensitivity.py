"""
Mesh Sensitivity Analysis
"""

import numpy as np
import preamble_functions as pre
import matplotlib.pyplot as plt


plt.rcParams["font.family"] = "cmr10" #Set Graph fonts to cmr10

#Import Apollo data for ZL 0.2m
Apollo_FileList_2d = pre.FileAddressList(r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\meshes\*.txt")
Apollo_gtable_2d = pre.FileAddressList(r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\meshes\*gtable",1)
Apollo_gauges_2d = pre.FileAddressList(r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\meshes\*gauges",1)
Apollo_log_2d = pre.FileAddressList(r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\meshes\*_log")
ult_cell_sizes_2d=[]
peak_impulse_2d = []
CPU_times_2d=[]


total_mesh_vol = 0.18* 0.08 * 0.08

for i in enumerate(Apollo_FileList_2d):
    ult_cell_sizes_2d.append(pre.ElementSizes(Apollo_FileList_2d[i[0]]))
    CPU_times_2d.append(pre.CPUFinder(Apollo_log_2d[i[0]]))
    peak_impulse_2d.append(max(Apollo_gtable_2d[i[0]][:,7]))
    

ult_cell_sizes_2d = np.asarray(ult_cell_sizes_2d)
CPU_times_2d = np.asarray(CPU_times_2d)    
peak_impulse_2d = np.asarray(peak_impulse_2d)
no_cells = np.divide(total_mesh_vol, ult_cell_sizes_2d**3)



#No of cells vs impulse & comp time
fig, [ax0,ax1] = plt.subplots(1,2)
ax0.scatter(no_cells[0:5], peak_impulse_2d[0:5], marker="o", s=15, label = 'Zone Length 0.02m')
ax0.scatter(no_cells[5::], peak_impulse_2d[5::], marker="o", s=15, label = 'Zone Length 0.01m')
ax0.set_xscale('log')
ax0.set_xlabel('Number of elements')
ax0.set_ylabel('peak specific impulse')
ax0.set_ylim(0,40000)
ax0.grid(which='minor', alpha=0.2)
ax0.grid(which='major', alpha=0.5)
ax1.scatter(no_cells[0:5], CPU_times_2d[0:5], marker="o", s=15, label = 'Zone Length 0.02m')
ax1.scatter(no_cells[5::], CPU_times_2d[5::], marker="o", s=15, label = 'Zone Length 0.01m')
ax1.set_xscale('log')
#ax1.set_ylim(0,10000)
ax1.set_xlabel('Number of elements')
ax1.set_ylabel('Wall time (s)')
ax1.grid(which='minor', alpha=0.2)
ax1.grid(which='major', alpha=0.5)
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='center', bbox_to_anchor=(0.25, 0.80), prop={'size':6})


#Impulse & OP analysis  @ 0 degrees
fig1, [[ax0,ax1],[ax2,ax3], [ax4, ax5]] = plt.subplots(3,2)
ax0.plot(Apollo_gauges_2d[2][:,0]*1000, Apollo_gauges_2d[2][:,1], label = '5mm')
ax0.plot(Apollo_gauges_2d[3][:,0]*1000, Apollo_gauges_2d[3][:,1], label = '2.5mm')
ax0.plot(Apollo_gauges_2d[4][:,0]*1000, Apollo_gauges_2d[4][:,1], label = '1.25mm')
ax0.plot(Apollo_gauges_2d[9][:,0]*1000, Apollo_gauges_2d[9][:,1], label = '0.625mm')
ax0.set_xlim(0,0.08)
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.75, 0.75), prop={'size':6})
ax0.set_title('0 degrees')
ax1.plot(Apollo_gauges_2d[2][:,0]*1000, Apollo_gauges_2d[2][:,201])
ax1.plot(Apollo_gauges_2d[3][:,0]*1000, Apollo_gauges_2d[3][:,201])
ax1.plot(Apollo_gauges_2d[4][:,0]*1000, Apollo_gauges_2d[4][:,201])
ax1.plot(Apollo_gauges_2d[9][:,0]*1000, Apollo_gauges_2d[9][:,201])
ax1.set_xlim(0,0.08)


#Impulse & OP analysis  @ 40 degrees
ax2.plot(Apollo_gauges_2d[2][:,0]*1000, Apollo_gauges_2d[2][:,101])
ax2.plot(Apollo_gauges_2d[3][:,0]*1000, Apollo_gauges_2d[3][:,101])
ax2.plot(Apollo_gauges_2d[4][:,0]*1000, Apollo_gauges_2d[4][:,101])
ax2.plot(Apollo_gauges_2d[9][:,0]*1000, Apollo_gauges_2d[9][:,101])
ax2.set_xlim(0,0.08)
ax2.set_title('40 degrees')

ax3.plot(Apollo_gauges_2d[2][:,0]*1000, Apollo_gauges_2d[2][:,300])
ax3.plot(Apollo_gauges_2d[3][:,0]*1000, Apollo_gauges_2d[3][:,300])
ax3.plot(Apollo_gauges_2d[4][:,0]*1000, Apollo_gauges_2d[4][:,300])
ax3.plot(Apollo_gauges_2d[9][:,0]*1000, Apollo_gauges_2d[9][:,300])
ax3.set_xlim(0,0.08)

#Impulse & OP analysis  @ 80 degrees
ax4.plot(Apollo_gauges_2d[2][:,0]*1000, Apollo_gauges_2d[2][:,200])
ax4.plot(Apollo_gauges_2d[3][:,0]*1000, Apollo_gauges_2d[3][:,200])
ax4.plot(Apollo_gauges_2d[4][:,0]*1000, Apollo_gauges_2d[4][:,200])
ax4.plot(Apollo_gauges_2d[9][:,0]*1000, Apollo_gauges_2d[9][:,200])
ax4.set_xlim(0,0.08)
ax4.set_title('80 degrees')

ax5.plot(Apollo_gauges_2d[2][:,0]*1000, Apollo_gauges_2d[2][:,400])
ax5.plot(Apollo_gauges_2d[3][:,0]*1000, Apollo_gauges_2d[3][:,400])
ax5.plot(Apollo_gauges_2d[4][:,0]*1000, Apollo_gauges_2d[4][:,400])
ax5.plot(Apollo_gauges_2d[9][:,0]*1000, Apollo_gauges_2d[9][:,400])
ax5.set_xlim(0,0.08)
