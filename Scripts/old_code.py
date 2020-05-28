# -*- coding: utf-8 -*-
"""
Created on Tue May 12 16:06:18 2020

@author: cip18jjp
"""


#Upper bound mesh sensitivity z = 0.5 -----------------------------------------

#Import Apollo data for mesh strategy
Apollo_FileList_z0_5 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_5\*.txt")
Apollo_gtable_z0_5 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_5\*gtable",1)
Apollo_gauges_z0_5 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_5\*gauges",1)
Apollo_log_z0_5 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_5\*_log")
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
fig3.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_1.pdf', format = 'pdf')


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
fig4.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_2.pdf', format = 'pdf')

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
fig4a.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_2a.pdf', format = 'pdf')

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
fig4b.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_2b.pdf', format = 'pdf')

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
fig4c.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_2c.pdf', format = 'pdf')

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
fig4d.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\mesh_convergence_z0_5_2d.pdf', format = 'pdf')



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
ax0.scatter(theta_exp_80mm, np.divide(MxI_1_80mm, max(MxI_1_80mm))/(0.1**(1/3)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = '80mm Exp')
ax0.scatter(theta_exp_80mm, np.divide(MxI_2_80mm, max(MxI_2_80mm))/(0.1**(1/3)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax0.scatter(theta_exp_80mm, np.divide(MxI_3_80mm, max(MxI_3_80mm))/(0.1**(1/3)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax0.scatter(theta_exp_80mm_mean, np.divide(Mx_mean_80mm, max(Mx_mean_80mm))/(0.1**(1/3)), marker="o", s=15., label = '80mm Exp Mean')

handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.7, 0.80), prop={'size':6})
plt.tight_layout()











#Compare 
def graph_impulse_comparisons():      
    
    RPB_MCEER_exp = RPB_MCEER_i(0.1*1.2, 0.08, 80)
    i_jang = jang(np.linspace(0,80, num=80), 1, 0.8)
    i_dharmasena = dharmasena(np.linspace(0,80,num=80), 1)
    
    #Graph of RPB & MCEER contour
    fig2, ax = plt.subplots(1,1)
    fig2.set_size_inches(3, 2.5)
    CS = ax.contourf(RPB_MCEER_exp[0], RPB_MCEER_exp[1], RPB_MCEER_exp[3], levels = [0,0.25,0.5,0.75,1, 1.5, 2, 2.5, 3, 3.5, 4, 5], cmap = plt.cm.magma_r)
    cbar = fig2.colorbar(CS)
    cbar.ax.set_ylabel('peak specific impulse (MPa.ms)')
    ax.set_ylabel('x-position')
    ax.set_xlabel('y-position')
    plt.tight_layout()
    fig2.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\theta_peak_impulse_RPB_MCEER_plate.pdf"), format = 'pdf')



    #Graph comparing different models       
    fig3, ax1 = plt.subplots(1,1)
    fig3.set_size_inches(2.5,2.5)    
    ax1.set_xlabel('angle of incidence')
    ax1.set_ylabel(r'$I_r/I_{r,max}$')
    
    ax1.plot(RPB_MCEER_exp[2][int(len(RPB_MCEER_exp[2])/2), 0:int(len(RPB_MCEER_exp[2])/2)], RPB_MCEER_exp[3][int(len(RPB_MCEER_exp[2])/2), 0:int(len(RPB_MCEER_exp[2])/2)]/ max(RPB_MCEER_exp[3][int(len(RPB_MCEER_exp[2])/2), 0:int(len(RPB_MCEER_exp[2])/2)]), 'k', label = 'RPB-MCEER')    
    ax1.plot(np.linspace(0,80,num=200), Apollo_gtable[0][:,7]/max(Apollo_gtable[0][:,7]), 'k-.', dashes=[12,6,12,6,3,6], label='CFD - Z = 0.05')
    ax1.plot(np.linspace(0,80,num=400), Henrych_i_fit(np.linspace(0,80,num=400), hen_mod['x'][0]), 'k:', label = 'Henrych')    
    ax1.plot(np.linspace(0,80, num=80), i_jang, 'k--', marker='o', markevery=12, ms=3., mfc = 'white', label = 'Jang')    
    ax1.plot(np.linspace(0,80, num=80), i_dharmasena, 'k-.', marker='D', markevery=16, ms=3., label = 'Dharmasena')    
    ax1.plot(theta.mean(1), gaussmod[int(len(gaussmod)/2)::], 'k:', label = 'Gaussian-single')    
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels, loc='upper left', prop={'size':6})
    plt.tight_layout()  
    fig3.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\theta_peak_impulse_i_theta_comparisons.pdf"), format = 'pdf')
#graph_impulse_comparisons()