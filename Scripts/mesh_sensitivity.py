"""
Mesh Sensitivity Analysis
"""

import numpy as np
import preamble_functions as pre
import matplotlib.pyplot as plt


plt.rcParams["font.family"] = "cmr10" #Set Graph fonts to cmr10

#Import Apollo data for ZL 0.2m
Apollo_FileList_no_autostages = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Mesh_Sensitivity_no_autostages\*.txt")
Apollo_gtable_no_autostages = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Mesh_Sensitivity_no_autostages\*gtable",1)
Apollo_log_no_autostages = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Mesh_Sensitivity_no_autostages\*_log")
ult_cell_sizes_no_autostages=[]
peak_impulse_no_autostages = []
CPU_times_no_autostages=[]
Apollo_FileList_autostages = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Mesh_Sensitivity_autostages\*.txt")
Apollo_gtable_autostages = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Mesh_Sensitivity_autostages\*gtable",1)
Apollo_log_autostages = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Mesh_Sensitivity_autostages\*_log")
ult_cell_sizes_autostages=[]
peak_impulse_autostages = []
CPU_times_autostages=[]

for i in enumerate(Apollo_FileList_no_autostages):
    ult_cell_sizes_no_autostages.append(pre.ElementSizes(i[1]))
    CPU_times_no_autostages.append(pre.CPUFinder(Apollo_log_no_autostages[i[0]]))
    peak_impulse_no_autostages.append(max(Apollo_gtable_no_autostages[i[0]][:,7]))
    ult_cell_sizes_autostages.append(pre.ElementSizes(Apollo_FileList_autostages[i[1]]))
    CPU_times_autostages.append(pre.CPUFinder(Apollo_log_autostages[i[0]]))
    peak_impulse_autostages.append(max(Apollo_gtable_autostages[i[0]][:,7]))
    
ult_cell_sizes_no_autostages = np.asarray(ult_cell_sizes_no_autostages)
CPU_times_no_autostages = np.asarray(CPU_times_no_autostages)    
peak_impulse_no_autostages = np.asarray(peak_impulse_no_autostages)
ult_cell_sizes_autostages = np.asarray(ult_cell_sizes_autostages)
CPU_times_autostages = np.asarray(CPU_times_autostages)    
peak_impulse_autostages = np.asarray(peak_impulse_autostages)
    

#Cell size vs peak impulse
fig, [ax0,ax1] = plt.subplots(1,2)
ax0.scatter(ult_cell_sizes_no_autostages, peak_impulse_no_autostages, marker="x", s=15, label = '0.2m ZL - No Autostage')
ax0.scatter(ult_cell_sizes_autostages, peak_impulse_autostages, marker="o", s=15, label = '0.2m ZL - Autostage')
ax0.set_xlabel('Ultimate cell size (m)')
ax0.set_ylabel('peak specific impulse')
ax1.scatter(ult_cell_sizes_no_autostages, CPU_times_no_autostages)
ax1.set_xlabel('Ultimate cell size (m)')
ax1.set_ylabel('Wall time (s)')