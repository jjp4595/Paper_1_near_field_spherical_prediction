"""
Mesh Sensitivity Analysis
"""
#Global info ------------------------------------------------------------------
import numpy as np
import preamble_functions as pre
import matplotlib.pyplot as plt
import scipy.io as sio

cr = 0.0246

#plt.rcParams["font.family"] = "cmr10" #Set Graph fonts to cmr10
plt.rc('font', family = 'serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')


#Load Data for 80mm and 380mm Apollo Experimental
fileID_NF_80mm_gtable = r"C:\Users\jorda\Google Drive\Apollo Sims\Near Field Sims\Sims\Latest\80mm_with_afterburn\*gtable"
gtable_80mm = pre.FileAddressList(fileID_NF_80mm_gtable,1)
Irmax_Ii_80mm = np.divide(gtable_80mm[0][:,7], gtable_80mm[0][:,7].max())
NF_80mm_exp = sio.loadmat(r"C:\Users\jorda\Google Drive\Apollo Sims\Near Field Sims\100gPE4Sphere_80mm") 
coords = np.append(np.arange(-100,101,25), np.arange(-100,-24,25)) 
coords = np.append(coords, np.arange(25,101,25))
coords = np.divide(coords,1000)
coords_mean = np.arange(0,0.101,0.025)
theta_exp_80mm = np.rad2deg(np.arctan(np.abs(coords)/0.08))
theta_exp_80mm_mean = np.rad2deg(np.arctan(np.abs(coords_mean)/0.08))
MxI_1_80mm = np.transpose(NF_80mm_exp['MxI'][:,:,0])
MxI_2_80mm = np.transpose(NF_80mm_exp['MxI'][:,:,1])
MxI_3_80mm = np.transpose(NF_80mm_exp['MxI'][:,:,2])
Mx_mean_80mm = np.transpose(NF_80mm_exp['MEANI'])
#------------------------------------------------------------------------------

"""



"""
#Lower bound mesh sensitivity z = 0.055 ---------------------------------------
#Import Apollo data for mesh strategy
Apollo_FileList_z0_055 = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_055\*.txt")
Apollo_gtable_z0_055 = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_055\*gtable",1)
Apollo_gauges_z0_055 = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_055\*gauges",1)
Apollo_log_z0_055 = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_055\*_log")
ult_cell_sizes_z0_055=[]
peak_impulse_z0_055 = []
CPU_times_z0_055=[]

total_mesh_vol = 0.35* 0.35 * 0.35

for i in enumerate(Apollo_FileList_z0_055):
    ult_cell_sizes_z0_055.append(pre.ElementSizes(Apollo_FileList_z0_055[i[0]]))
    CPU_times_z0_055.append(pre.CPUFinder(Apollo_log_z0_055[i[0]]))
    peak_impulse_z0_055.append(max(Apollo_gtable_z0_055[i[0]][:,7]))
    

ult_cell_sizes_z0_055 = np.asarray(ult_cell_sizes_z0_055)
R_z0_055 = np.divide(0.0255, ult_cell_sizes_z0_055)
CPU_times_z0_055 = np.asarray(CPU_times_z0_055)    
peak_impulse_z0_055 = np.asarray(peak_impulse_z0_055)
no_cells = np.divide(total_mesh_vol, ult_cell_sizes_z0_055**3)

#Indexes
final_100s = 4 + 1
final_200s = 8 + 1

#Graph 1 ---------------------------------------------------------------------
#No of cells vs impulse & comp time
fig, [ax2,ax0,ax1] = plt.subplots(1,3)
fig.set_size_inches(7, 3)
ax2.scatter(R_z0_055[0:final_100s], peak_impulse_z0_055[0:final_100s], c = 'k', marker=".", s=15, label = 'Zone Length 0.05m')
ax2.scatter(R_z0_055[final_100s:final_200s], peak_impulse_z0_055[final_100s:final_200s], c = 'k', marker="*", s=15, label = 'Zone Length 0.02m')
ax2.scatter(R_z0_055[final_200s::], peak_impulse_z0_055[final_200s::], c = 'k', marker="+", s=15, label = 'Zone Length 0.01m')
ax2.set_xlabel('R / cell length')
ax2.set_ylabel('peak specific impulse')
ax2.set_ylim(0,12000)
ax2.grid(which='minor', alpha=0.2)
ax2.grid(which='major', alpha=0.5)
ax0.scatter(no_cells[0:final_100s], peak_impulse_z0_055[0:final_100s], c = 'k', marker=".", s=15, label = 'Zone Length 0.05m')
ax0.scatter(no_cells[final_100s:final_200s], peak_impulse_z0_055[final_100s:final_200s], c = 'k', marker="*", s=15, label = 'Zone Length 0.02m')
ax0.scatter(no_cells[final_200s::], peak_impulse_z0_055[final_200s::], c = 'k', marker="+", s=15, label = 'Zone Length 0.01m')
ax0.set_xscale('log')
ax0.set_xlabel('Number of elements')
ax0.set_ylabel('peak specific impulse')
ax0.set_ylim(0,12000)
ax0.grid(which='minor', alpha=0.2)
ax0.grid(which='major', alpha=0.5)
ax1.scatter(no_cells[0:final_100s], CPU_times_z0_055[0:final_100s], c = 'k', marker=".", s=15, label = 'Zone Length 0.05m')
ax1.scatter(no_cells[final_100s:final_200s], CPU_times_z0_055[final_100s:final_200s], c = 'k', marker="*", s=15, label = 'Zone Length 0.02m')
ax1.scatter(no_cells[final_200s::], CPU_times_z0_055[final_200s::], c = 'k', marker="+", s=15, label = 'Zone Length 0.01m')
ax1.set_xscale('log')
#ax1.set_ylim(0,10000)
ax1.set_xlabel('Number of elements')
ax1.set_ylabel('Wall time (s)')
ax1.grid(which='minor', alpha=0.2)
ax1.grid(which='major', alpha=0.5)
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='center', bbox_to_anchor=(0.40, 0.80), prop={'size':6})
plt.tight_layout()
#fig.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_055_1.pdf', format = 'pdf')



# Graph 2 --------------------------------------------------------------------
#Impulse & OP analysis  @ 0 degrees
fig1, [[ax0,ax1],[ax2,ax3], [ax4, ax5]] = plt.subplots(3,2)
fig1.set_size_inches(7, 3)
ax0.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,1], c = 'k', label = '3.125mm')
ax0.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,1], 'k--' , label = '2.5mm')
ax0.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,1], 'k:', label = '2.5mm')
ax0.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,1], 'k-.', label = '1.25mm')
ax0.set_xlim(0,0.12)
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.75, 0.75), prop={'size':6})
ax0.set_title('0 degrees')
ax1.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,201], c = 'k')
ax1.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,201], 'k--')
ax1.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,201], 'k:')
ax1.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,201], 'k-.')
ax1.set_xlim(0,0.12)

#Impulse & OP analysis  @ 40 degrees
ax2.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,101], c = 'k')
ax2.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,101], 'k--')
ax2.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,101], 'k:')
ax2.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,101], 'k-.')
ax2.set_xlim(0,0.12)
ax2.set_title('40 degrees')

ax3.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,300], c = 'k')
ax3.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,300], 'k--')
ax3.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,300], 'k:')
ax3.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,300], 'k-.')
ax3.set_xlim(0,0.12)

#Impulse & OP analysis  @ 80 degrees
ax4.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,200], c = 'k')
ax4.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,200], 'k--')
ax4.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,200], 'k:')
ax4.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,200], 'k-.')
ax4.set_xlim(0,0.12)
ax4.set_title('80 degrees')

ax5.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,400], c = 'k')
ax5.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,400], 'k--')
ax5.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,400], 'k:')
ax5.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,400], 'k-.')
ax5.set_xlim(0,0.12)
plt.tight_layout()
#fig1.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_055_2.pdf', format = 'pdf')
#-----------------------------------------------------------------------------

#Graph 3 ---------------------------------------------------------------------
#theta vs scaled I plots 
theta = np.linspace(0,80,200)
fig2, [ax, ax0] = plt.subplots(1,2)
fig2.set_size_inches(5, 3)
ax.set_xlabel('theta (degrees)')
ax.set_ylabel('peak specific scaled impulse')
ax.plot(theta, Apollo_gtable_z0_055[4][:,7]/(0.1**(1/3)), 'k', label = '3.125mm')
ax.plot(theta, Apollo_gtable_z0_055[8][:,7]/(0.1**(1/3)), 'k--', label = '2.5mm')
ax.plot(theta, Apollo_gtable_z0_055[11][:,7]/(0.1**(1/3)), 'k:', label = '2.5mm')
ax.plot(theta, Apollo_gtable_z0_055[12][:,7]/(0.1**(1/3)), 'k-.', label = '1.25mm')
ax0.set_xlabel('theta (degrees)')
ax0.set_ylabel('peak impulse ratio of maximum')
ax0.plot(theta, (Apollo_gtable_z0_055[4][:,7]/(0.1**(1/3)))/max(Apollo_gtable_z0_055[4][:,7]/(0.1**(1/3))), 'k', label = '3.125mm')
ax0.plot(theta, (Apollo_gtable_z0_055[8][:,7]/(0.1**(1/3)))/max(Apollo_gtable_z0_055[8][:,7]/(0.1**(1/3))), 'k--', label = '2.5mm')
ax0.plot(theta, (Apollo_gtable_z0_055[11][:,7]/(0.1**(1/3)))/max(Apollo_gtable_z0_055[11][:,7]/(0.1**(1/3))), 'k:', label = '2.5mm')
ax0.plot(theta, (Apollo_gtable_z0_055[12][:,7]/(0.1**(1/3)))/max(Apollo_gtable_z0_055[12][:,7]/(0.1**(1/3))), 'k-.', label = '1.25mm')
ax0.scatter(theta_exp_80mm, np.divide(MxI_1_80mm, max(MxI_1_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = '80mm Exp')
ax0.scatter(theta_exp_80mm, np.divide(MxI_2_80mm, max(MxI_2_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax0.scatter(theta_exp_80mm, np.divide(MxI_3_80mm, max(MxI_3_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax0.scatter(theta_exp_80mm_mean, np.divide(Mx_mean_80mm, max(Mx_mean_80mm)), marker="o", s=15., label = '80mm Exp Mean')
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.7, 0.80), prop={'size':6})
plt.tight_layout()
#fig2.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_05_3.pdf', format = 'pdf')
#-----------------------------------------------------------------------------







#Upper bound mesh sensitivity z = 0.5 -----------------------------------------

#Import Apollo data for mesh strategy
Apollo_FileList_z0_5 = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_5\*.txt")
Apollo_gtable_z0_5 = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_5\*gtable",1)
Apollo_gauges_z0_5 = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_5\*gauges",1)
Apollo_log_z0_5 = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_5\*_log")
ult_cell_sizes_z0_5=[]
peak_impulse_z0_5 = []
CPU_times_z0_5=[]

total_mesh_vol = 1.75* 1.75 * 1.75

for i in enumerate(Apollo_FileList_z0_5):
    ult_cell_sizes_z0_5.append(pre.ElementSizes(Apollo_FileList_z0_5[i[0]]))
    CPU_times_z0_5.append(pre.CPUFinder(Apollo_log_z0_5[i[0]]))
    peak_impulse_z0_5.append(max(Apollo_gtable_z0_5[i[0]][:,7]))
    

ult_cell_sizes_z0_5 = np.asarray(ult_cell_sizes_z0_5)
R_z0_5 = np.divide(0.2321, ult_cell_sizes_z0_5)
CPU_times_z0_5 = np.asarray(CPU_times_z0_5)    
peak_impulse_z0_5 = np.asarray(peak_impulse_z0_5)
no_cells_z0_5 = np.divide(total_mesh_vol, ult_cell_sizes_z0_5**3)

#Indexes
final_100s = 4 + 1
final_200s = 8 + 1


#Graph 1 ---------------------------------------------------------------------
#No of cells vs impulse & comp time
fig3, [ax2,ax0,ax1] = plt.subplots(1,3)
fig3.set_size_inches(7, 3)
ax2.scatter(R_z0_5[0:final_100s], peak_impulse_z0_5[0:final_100s], c = 'k', marker=".", s=15, label = 'Zone Length 0.05m')
ax2.scatter(R_z0_5[final_100s:final_200s], peak_impulse_z0_5[final_100s:final_200s], c = 'k', marker="*", s=15, label = 'Zone Length 0.02m')
ax2.scatter(R_z0_5[final_200s::], peak_impulse_z0_5[final_200s::], c = 'k', marker="*", s=15, label = 'Zone Length 0.02m')
ax2.set_xlabel('R / cell length')
ax2.set_ylabel('peak specific impulse')
ax2.set_ylim(0,1000)
ax2.grid(which='minor', alpha=0.2)
ax2.grid(which='major', alpha=0.5)
ax0.scatter(no_cells_z0_5[0:final_100s], peak_impulse_z0_5[0:final_100s], c = 'k', marker=".", s=15, label = 'Zone Length 0.05m')
ax0.scatter(no_cells_z0_5[final_100s:final_200s], peak_impulse_z0_5[final_100s:final_200s], c = 'k', marker="*", s=15, label = 'Zone Length 0.02m')
ax0.scatter(no_cells_z0_5[final_200s::], peak_impulse_z0_5[final_200s::], c = 'k', marker="+", s=15, label = 'Zone Length 0.01m')
ax0.set_xscale('log')
ax0.set_xlabel('Number of elements')
ax0.set_ylabel('peak specific impulse')
ax0.set_ylim(0,1000)
ax0.grid(which='minor', alpha=0.2)
ax0.grid(which='major', alpha=0.5)
ax1.scatter(no_cells_z0_5[0:final_100s], CPU_times_z0_5[0:final_100s], c = 'k', marker=".", s=15, label = 'Zone Length 0.05m')
ax1.scatter(no_cells_z0_5[final_100s:final_200s], CPU_times_z0_5[final_100s:final_200s], c = 'k', marker="*", s=15, label = 'Zone Length 0.02m')
ax1.scatter(no_cells_z0_5[final_200s::], CPU_times_z0_5[final_200s::], c = 'k', marker="+", s=15, label = 'Zone Length 0.01m')
ax1.set_xscale('log')
#ax1.set_ylim(0,10000)
ax1.set_xlabel('Number of elements')
ax1.set_ylabel('Wall time (s)')
ax1.grid(which='minor', alpha=0.2)
ax1.grid(which='major', alpha=0.5)
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='center', bbox_to_anchor=(0.40, 0.80), prop={'size':6})
plt.tight_layout()



# # Graph 2 --------------------------------------------------------------------
# #Impulse & OP analysis  @ 0 degrees
# fig4, [[ax0,ax1],[ax2,ax3], [ax4, ax5]] = plt.subplots(3,2)
# fig4.set_size_inches(7, 3)
# ax0.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,1], c = 'k', label = '3.125mm')
# ax0.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,1], 'k--' , label = '2.5mm')
# ax0.plot(Apollo_gauges_z0_5[11][:,0]*1000, Apollo_gauges_z0_5[11][:,1], 'k:', label = '2.5mm')
# ax0.plot(Apollo_gauges_z0_5[12][:,0]*1000, Apollo_gauges_z0_5[12][:,1], 'k-.', label = '1.25mm')
# ax0.set_xlim(0,0.12)
# handles, labels = ax0.get_legend_handles_labels()
# ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.75, 0.75), prop={'size':6})
# ax0.set_title('0 degrees')
# ax1.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,201], c = 'k')
# ax1.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,201], 'k--')
# ax1.plot(Apollo_gauges_z0_5[11][:,0]*1000, Apollo_gauges_z0_5[11][:,201], 'k:')
# ax1.plot(Apollo_gauges_z0_5[12][:,0]*1000, Apollo_gauges_z0_5[12][:,201], 'k-.')
# ax1.set_xlim(0,0.12)

# #Impulse & OP analysis  @ 40 degrees
# ax2.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,101], c = 'k')
# ax2.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,101], 'k--')
# ax2.plot(Apollo_gauges_z0_5[11][:,0]*1000, Apollo_gauges_z0_5[11][:,101], 'k:')
# ax2.plot(Apollo_gauges_z0_5[12][:,0]*1000, Apollo_gauges_z0_5[12][:,101], 'k-.')
# ax2.set_xlim(0,0.12)
# ax2.set_title('40 degrees')

# ax3.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,300], c = 'k')
# ax3.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,300], 'k--')
# ax3.plot(Apollo_gauges_z0_5[11][:,0]*1000, Apollo_gauges_z0_5[11][:,300], 'k:')
# ax3.plot(Apollo_gauges_z0_5[12][:,0]*1000, Apollo_gauges_z0_5[12][:,300], 'k-.')
# ax3.set_xlim(0,0.12)

# #Impulse & OP analysis  @ 80 degrees
# ax4.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,200], c = 'k')
# ax4.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,200], 'k--')
# ax4.plot(Apollo_gauges_z0_5[11][:,0]*1000, Apollo_gauges_z0_5[11][:,200], 'k:')
# ax4.plot(Apollo_gauges_z0_5[12][:,0]*1000, Apollo_gauges_z0_5[12][:,200], 'k-.')
# ax4.set_xlim(0,0.12)
# ax4.set_title('80 degrees')

# ax5.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,400], c = 'k')
# ax5.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,400], 'k--')
# ax5.plot(Apollo_gauges_z0_5[11][:,0]*1000, Apollo_gauges_z0_5[11][:,400], 'k:')
# ax5.plot(Apollo_gauges_z0_5[12][:,0]*1000, Apollo_gauges_z0_5[12][:,400], 'k-.')
# ax5.set_xlim(0,0.12)
# plt.tight_layout()
# #-----------------------------------------------------------------------------

# #Graph 3 ---------------------------------------------------------------------
# #theta vs scaled I plots 
fig5, [ax, ax0] = plt.subplots(1,2)
fig5.set_size_inches(5, 3)
ax.set_xlabel('theta (degrees)')
ax.set_ylabel('peak specific scaled impulse')
ax.plot(theta, Apollo_gtable_z0_5[4][:,7]/(0.1**(1/3)), 'k', label = 'mm')
ax.plot(theta, Apollo_gtable_z0_5[8][:,7]/(0.1**(1/3)), 'k--', label = '2.5mm')
ax0.set_xlabel('theta (degrees)')
ax0.set_ylabel('peak impulse ratio of maximum')
ax0.plot(theta, (Apollo_gtable_z0_5[4][:,7]/(0.1**(1/3)))/max(Apollo_gtable_z0_5[4][:,7]/(0.1**(1/3))), 'k', label = '3.125mm')
ax0.plot(theta, (Apollo_gtable_z0_5[8][:,7]/(0.1**(1/3)))/max(Apollo_gtable_z0_5[8][:,7]/(0.1**(1/3))), 'k--', label = '2.5mm')
ax0.scatter(theta_exp_80mm, np.divide(MxI_1_80mm, max(MxI_1_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = '80mm Exp')
ax0.scatter(theta_exp_80mm, np.divide(MxI_2_80mm, max(MxI_2_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax0.scatter(theta_exp_80mm, np.divide(MxI_3_80mm, max(MxI_3_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax0.scatter(theta_exp_80mm_mean, np.divide(Mx_mean_80mm, max(Mx_mean_80mm)), marker="o", s=15., label = '80mm Exp Mean')
ax2.set_ylim(0,1)
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.7, 0.80), prop={'size':6})
plt.tight_layout()
# #-----------------------------------------------------------------------------

