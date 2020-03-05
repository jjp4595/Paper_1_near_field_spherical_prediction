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
params = {'font.family':'serif',
        'axes.labelsize':'small',
        'xtick.labelsize':'x-small',
        'ytick.labelsize':'x-small', 
        'legend.fontsize':'small',
        'grid.linestyle':'--',
        'grid.linewidth':'0.5',
        'lines.linewidth':'0.5'}
plt.rcParams.update(params)


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
R_z0_055 = np.divide(0.0255+cr, ult_cell_sizes_z0_055)
CPU_times_z0_055 = np.asarray(CPU_times_z0_055)    
peak_impulse_z0_055 = np.asarray(peak_impulse_z0_055)
no_cells = np.divide(total_mesh_vol, ult_cell_sizes_z0_055**3)

#Indexes
final_100s = 4 + 1
final_200s = 8 + 1

temp = max(peak_impulse_z0_055[final_200s::]/1e3)*1.03
#Graph 1 ---------------------------------------------------------------------
#No of cells vs impulse & comp time
fig, [ax0,ax1] = plt.subplots(1,2)
fig.set_size_inches(5, 2.5)
ax0.scatter(R_z0_055[0:final_100s], peak_impulse_z0_055[0:final_100s]/1e3, c = 'k', marker=".", s=10, label = 'Zone Length 0.05m')
ax0.scatter(R_z0_055[final_100s:final_200s], peak_impulse_z0_055[final_100s:final_200s]/1e3, c = 'r', marker="^", s=10, label = 'Zone Length 0.02m')
ax0.scatter(R_z0_055[final_200s::], peak_impulse_z0_055[final_200s::]/1e3, c = 'g', marker="D", s=10, label = 'Zone Length 0.01m')
ax0.plot([0,50], [temp, temp], linewidth = 0.5, c = 'k')
ax0.plot([0,50], [temp*0.9, temp*0.9], linewidth = 0.5, linestyle = '--', c = 'k',label = '$10\%$ convergence')
ax0.set_xlabel('R / cell length')
ax0.set_ylabel('peak specific impulse (MPa.ms)')
ax0.set_ylim(0,15)
ax0.grid(which='minor', alpha=0.2)
ax0.grid(which='major', alpha=0.5)
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.60, 0.40), prop={'size':6})
ax1.scatter(R_z0_055[0:final_100s], CPU_times_z0_055[0:final_100s], c = 'k', marker=".", s=10, label = 'Zone Length 0.05m')
ax1.scatter(R_z0_055[final_100s:final_200s], CPU_times_z0_055[final_100s:final_200s], c = 'r', marker="^", s=10, label = 'Zone Length 0.02m')
ax1.scatter(R_z0_055[final_200s::], CPU_times_z0_055[final_200s::], c = 'g', marker="D", s=10, label = 'Zone Length 0.01m')
ax1.set_xlabel('R / cell length')
ax1.set_ylabel('Wall time (s)')
ax1.grid(which='minor', alpha=0.2)
ax1.grid(which='major', alpha=0.5)
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=6)
ax1.locator_params(axis = 'both',tight=True, nbins=6)
fig.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_055_1.pdf', format = 'pdf')



# Graph 2 --------------------------------------------------------------------
#Impulse & OP analysis  @ 0 degrees
fig1, [ax0,ax1] = plt.subplots(1,2)
fig1.set_size_inches(5, 1.8)
ax0.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,1]/1e6, c = 'k', label = '3.125mm')
ax0.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,1]/1e6, 'k--' , label = '2.5mm')
ax0.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,1]/1e6, 'k:', label = '2.5mm')
ax0.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,1]/1e6, 'k-.', label = '1.25mm')
ax0.set_xlim(0,0.12)
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.7, 0.7), prop={'size':6})
ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('Overpressure (MPa)')
ax1.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,201]/1e3, c = 'k')
ax1.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,201]/1e3, 'k--')
ax1.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,201]/1e3, 'k:')
ax1.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,201]/1e3, 'k-.')
ax1.set_xlim(0,0.12)
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=6)
ax1.locator_params(axis = 'both',tight=True, nbins=6)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Impulse (MPa.ms)')
fig1.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_055_2.pdf', format = 'pdf')

#Impulse & OP analysis  @ 20 degrees
fig1a, [ax0,ax1] = plt.subplots(1,2)
fig1a.set_size_inches(5, 1.8)
ax0.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,51]/1e6, c = 'k')
ax0.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,51]/1e6, 'k--')
ax0.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,51]/1e6, 'k:')
ax0.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,51]/1e6, 'k-.')
ax0.set_xlim(0,0.12)
ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('Overpressure (MPa)')
ax1.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,251]/1e3, c = 'k')
ax1.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,251]/1e3, 'k--')
ax1.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,251]/1e3, 'k:')
ax1.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,251]/1e3, 'k-.')
ax1.set_xlim(0,0.12)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Impulse (MPa.ms)')
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=6)
ax1.locator_params(axis = 'both',tight=True, nbins=6)
fig1a.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_055_2a.pdf', format = 'pdf')

#Impulse & OP analysis  @ 40 degrees
fig1b, [ax0,ax1] = plt.subplots(1,2)
fig1b.set_size_inches(5, 1.8)
ax0.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,101]/1e6, c = 'k')
ax0.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,101]/1e6, 'k--')
ax0.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,101]/1e6, 'k:')
ax0.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,101]/1e6, 'k-.')
ax0.set_xlim(0,0.12)
ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('Overpressure (MPa)')
plt.tight_layout()
ax1.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,301]/1e3, c = 'k')
ax1.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,301]/1e3, 'k--')
ax1.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,301]/1e3, 'k:')
ax1.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,301]/1e3, 'k-.')
ax1.set_xlim(0,0.12)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Impulse (MPa.ms)')
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=6)
ax1.locator_params(axis = 'both',tight=True, nbins=6)
fig1b.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_055_2b.pdf', format = 'pdf')

#Impulse & OP analysis  @ 60 degrees
fig1c, [ax0,ax1] = plt.subplots(1,2)
fig1c.set_size_inches(5, 1.8)
ax0.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,151]/1e6, c = 'k')
ax0.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,151]/1e6, 'k--')
ax0.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,151]/1e6, 'k:')
ax0.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,151]/1e6, 'k-.')
ax0.set_xlim(0,0.12)
ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('Overpressure (MPa)')
ax1.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,351]/1e3, c = 'k')
ax1.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,351]/1e3, 'k--')
ax1.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,351]/1e3, 'k:')
ax1.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,351]/1e3, 'k-.')
ax1.set_xlim(0,0.12)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Impulse (MPa.ms)')
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=6)
ax1.locator_params(axis = 'both',tight=True, nbins=6)
fig1c.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_055_2c.pdf', format = 'pdf')

#Impulse & OP analysis  @ 80 degrees
fig1d, [ax0,ax1] = plt.subplots(1,2)
fig1d.set_size_inches(5, 1.8)
ax0.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,200]/1e6, c = 'k')
ax0.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,200]/1e6, 'k--')
ax0.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,200]/1e6, 'k:')
ax0.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,200]/1e6, 'k-.')
ax0.set_xlim(0,0.12)
ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('Overpressure (MPa)')
ax1.plot(Apollo_gauges_z0_055[4][:,0]*1000, Apollo_gauges_z0_055[4][:,400]/1e3, c = 'k')
ax1.plot(Apollo_gauges_z0_055[8][:,0]*1000, Apollo_gauges_z0_055[8][:,400]/1e3, 'k--')
ax1.plot(Apollo_gauges_z0_055[11][:,0]*1000, Apollo_gauges_z0_055[11][:,400]/1e3, 'k:')
ax1.plot(Apollo_gauges_z0_055[12][:,0]*1000, Apollo_gauges_z0_055[12][:,400]/1e3, 'k-.')
ax1.set_xlim(0,0.12)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Impulse (MPa.ms)')
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=6)
ax1.locator_params(axis = 'both',tight=True, nbins=6)
fig1d.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_055_2d.pdf', format = 'pdf')

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
fig2.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_055_3.pdf', format = 'pdf')
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
R_z0_5 = np.divide(0.2321+cr, ult_cell_sizes_z0_5)
CPU_times_z0_5 = np.asarray(CPU_times_z0_5)    
peak_impulse_z0_5 = np.asarray(peak_impulse_z0_5)
no_cells_z0_5 = np.divide(total_mesh_vol, ult_cell_sizes_z0_5**3)

#Indexes
final_100s = 4 + 1
final_200s = 8 + 1

temp = max(peak_impulse_z0_5[final_200s::]/1e3)*1.05
#Graph 1 ---------------------------------------------------------------------
#No of cells vs impulse & comp time
fig3, [ax0,ax1] = plt.subplots(1,2)
fig3.set_size_inches(5, 2.5)
ax0.scatter(R_z0_5[0:final_100s], peak_impulse_z0_5[0:final_100s]/1e3, c = 'k', marker=".", s=10, label = 'Zone Length 0.05m')
ax0.scatter(R_z0_5[final_100s:final_200s], peak_impulse_z0_5[final_100s:final_200s]/1e3, c = 'r', marker="^", s=10, label = 'Zone Length 0.02m')
ax0.scatter(R_z0_5[final_200s::], peak_impulse_z0_5[final_200s::]/1e3, c = 'g', marker="D", s=10, label = 'Zone Length 0.02m')
ax0.plot([0,50], [temp, temp], linewidth = 0.5, c = 'k')
ax0.plot([0,50], [temp*0.9, temp*0.9], linewidth = 0.5, linestyle = '--', c = 'k',label = '$10\%$ convergence')
ax0.set_xlabel('R / cell length')
ax0.set_ylabel('peak specific impulse (MPa.ms)')
ax0.grid(which='minor', alpha=0.2)
ax0.grid(which='major', alpha=0.5)
ax0.set_ylim(0,1)
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.60, 0.30), prop={'size':6})
ax1.scatter(R_z0_5[0:final_100s], CPU_times_z0_5[0:final_100s], c = 'k', marker=".", s=10, label = 'Zone Length 0.05m')
ax1.scatter(R_z0_5[final_100s:final_200s], CPU_times_z0_5[final_100s:final_200s], c = 'r', marker="^", s=10, label = 'Zone Length 0.02m')
ax1.scatter(R_z0_5[final_200s::], CPU_times_z0_5[final_200s::], c = 'g', marker="D", s=10, label = 'Zone Length 0.01m')
ax1.set_xlabel('R / cell length')
ax1.set_ylabel('Wall time (s)')
ax1.grid(which='minor', alpha=0.2)
ax1.grid(which='major', alpha=0.5)
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=6)
ax1.locator_params(axis = 'both',tight=True, nbins=6)
fig3.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_1.pdf', format = 'pdf')


# Graph 2 --------------------------------------------------------------------
#Impulse & OP analysis  @ 0 degrees
fig4, [ax0,ax1] = plt.subplots(1,2)
fig4.set_size_inches(5, 1.8)
ax0.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,1]/1e6, c = 'k', label = '6.25mm (Res. Lev. 4)')
ax0.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,1]/1e6, 'k--' , label = '6.25mm (Res. Lev. 3)')
ax0.plot(Apollo_gauges_z0_5[10][:,0]*1000, Apollo_gauges_z0_5[10][:,1]/1e6, 'k:', label = '10mm')
ax0.plot(Apollo_gauges_z0_5[7][:,0]*1000, Apollo_gauges_z0_5[7][:,1]/1e6, 'k-.', label = '12.5mm (Res. Lev. 2)')
ax0.set_xlim(0,1)
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.5, 0.7), prop={'size':6})
ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('Overpressure (MPa)')
ax1.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,201]/1e3, c = 'k')
ax1.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,201]/1e3, 'k--')
ax1.plot(Apollo_gauges_z0_5[10][:,0]*1000, Apollo_gauges_z0_5[10][:,201]/1e3, 'k:')
ax1.plot(Apollo_gauges_z0_5[7][:,0]*1000, Apollo_gauges_z0_5[7][:,201]/1e3, 'k-.')
ax1.set_xlim(0,1)
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=4)
ax1.locator_params(axis = 'both',tight=True, nbins=4)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Impulse (MPa.ms)')
fig4.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_2.pdf', format = 'pdf')

#Impulse & OP analysis  @ 20 degrees
fig4a, [ax0,ax1] = plt.subplots(1,2)
fig4a.set_size_inches(5, 1.8)
ax0.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,51]/1e6, c = 'k')
ax0.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,51]/1e6, 'k--')
ax0.plot(Apollo_gauges_z0_5[10][:,0]*1000, Apollo_gauges_z0_5[10][:,51]/1e6, 'k:')
ax0.plot(Apollo_gauges_z0_5[7][:,0]*1000, Apollo_gauges_z0_5[7][:,51]/1e6, 'k-.')
ax0.set_xlim(0,1)
ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('Overpressure (MPa)')
ax1.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,251]/1e3, c = 'k')
ax1.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,251]/1e3, 'k--')
ax1.plot(Apollo_gauges_z0_5[10][:,0]*1000, Apollo_gauges_z0_5[10][:,251]/1e3, 'k:')
ax1.plot(Apollo_gauges_z0_5[7][:,0]*1000, Apollo_gauges_z0_5[7][:,251]/1e3, 'k-.')
ax1.set_xlim(0,1)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Impulse (MPa.ms)')
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=4)
ax1.locator_params(axis = 'both',tight=True, nbins=4)
fig4a.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_2a.pdf', format = 'pdf')

#Impulse & OP analysis  @ 40 degrees
fig4b, [ax0,ax1] = plt.subplots(1,2)
fig4b.set_size_inches(5, 1.8)
ax0.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,101]/1e6, c = 'k')
ax0.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,101]/1e6, 'k--')
ax0.plot(Apollo_gauges_z0_5[10][:,0]*1000, Apollo_gauges_z0_5[10][:,101]/1e6, 'k:')
ax0.plot(Apollo_gauges_z0_5[7][:,0]*1000, Apollo_gauges_z0_5[7][:,101]/1e6, 'k-.')
ax0.set_xlim(0,1)
ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('Overpressure (MPa)')
plt.tight_layout()
ax1.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,301]/1e3, c = 'k')
ax1.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,301]/1e3, 'k--')
ax1.plot(Apollo_gauges_z0_5[10][:,0]*1000, Apollo_gauges_z0_5[10][:,301]/1e3, 'k:')
ax1.plot(Apollo_gauges_z0_5[7][:,0]*1000, Apollo_gauges_z0_5[7][:,301]/1e3, 'k-.')
ax1.set_xlim(0,1)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Impulse (MPa.ms)')
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=4)
ax1.locator_params(axis = 'both',tight=True, nbins=4)
fig4b.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_2b.pdf', format = 'pdf')

#Impulse & OP analysis  @ 60 degrees
fig4c, [ax0,ax1] = plt.subplots(1,2)
fig4c.set_size_inches(5, 1.8)
ax0.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,151]/1e6, c = 'k')
ax0.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,151]/1e6, 'k--')
ax0.plot(Apollo_gauges_z0_5[10][:,0]*1000, Apollo_gauges_z0_5[10][:,151]/1e6, 'k:')
ax0.plot(Apollo_gauges_z0_5[7][:,0]*1000, Apollo_gauges_z0_5[7][:,151]/1e6, 'k-.')
ax0.set_xlim(0,1)
ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('Overpressure (MPa)')
ax1.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,351]/1e3, c = 'k')
ax1.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,351]/1e3, 'k--')
ax1.plot(Apollo_gauges_z0_5[10][:,0]*1000, Apollo_gauges_z0_5[10][:,351]/1e3, 'k:')
ax1.plot(Apollo_gauges_z0_5[7][:,0]*1000, Apollo_gauges_z0_5[7][:,351]/1e3, 'k-.')
ax1.set_xlim(0,1)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Impulse (MPa.ms)')
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=4)
ax1.locator_params(axis = 'both',tight=True, nbins=4)
fig4c.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_2c.pdf', format = 'pdf')

#Impulse & OP analysis  @ 80 degrees
fig4d, [ax0,ax1] = plt.subplots(1,2)
fig4d.set_size_inches(5, 1.8)
ax0.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,200]/1e6, c = 'k')
ax0.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,200]/1e6, 'k--')
ax0.plot(Apollo_gauges_z0_5[10][:,0]*1000, Apollo_gauges_z0_5[10][:,200]/1e6, 'k:')
ax0.plot(Apollo_gauges_z0_5[7][:,0]*1000, Apollo_gauges_z0_5[7][:,200]/1e6, 'k-.')
ax0.set_xlim(0,3)
ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('Overpressure (MPa)')
ax1.plot(Apollo_gauges_z0_5[4][:,0]*1000, Apollo_gauges_z0_5[4][:,400]/1e3, c = 'k')
ax1.plot(Apollo_gauges_z0_5[8][:,0]*1000, Apollo_gauges_z0_5[8][:,400]/1e3, 'k--')
ax1.plot(Apollo_gauges_z0_5[10][:,0]*1000, Apollo_gauges_z0_5[10][:,400]/1e3, 'k:')
ax1.plot(Apollo_gauges_z0_5[7][:,0]*1000, Apollo_gauges_z0_5[7][:,400]/1e3, 'k-.')
ax1.set_xlim(0,3)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Impulse (MPa.ms)')
plt.tight_layout()
ax0.locator_params(axis = 'both',tight=True, nbins=4)
ax1.locator_params(axis = 'both',tight=True, nbins=4)
fig4d.savefig(r'C:\Users\jorda\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_2d.pdf', format = 'pdf')


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

handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.7, 0.80), prop={'size':6})
plt.tight_layout()
# #-----------------------------------------------------------------------------

