"""
Mesh Sensitivity Analysis
"""
#Global info ------------------------------------------------------------------
import numpy as np
import preamble_functions as pre
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker
import scipy.io as sio
import os 
from impulse_models import *
from scipy.signal import savgol_filter
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter

cr = 0.0246
charge_mass = 0.1
params = {'font.family':'serif',
        'axes.labelsize':'small',
        'xtick.labelsize':'x-small',
        'ytick.labelsize':'x-small', 
        'axes.linewidth':0.5,
        
        'xtick.major.width':0.5,
        'xtick.minor.width':0.4,
        'ytick.major.width':0.5,
        'ytick.minor.width':0.4,
        'xtick.major.size':3.0,
        'xtick.minor.size':1.5,
        'ytick.major.size':3.0,
        'ytick.minor.size':1.5,
        
        'legend.fontsize':'small',
        'legend.title_fontsize':'small',
        'legend.fancybox': False,
        'legend.framealpha': 1,
        'legend.shadow': False,
        'legend.frameon': True,
        'legend.edgecolor':'black',
        'patch.linewidth':0.5,
        
        'scatter.marker': 's',
        
        'grid.linewidth':'0.5',
        
        'lines.linewidth':'0.5'}
plt.rcParams.update(params)



#Load Data for 80mm and 380mm Apollo Experimental
fileID_NF_80mm_gtable = os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Near Field Sims\Sims\Latest\80mm_with_afterburn\*gtable"
gtable_80mm = pre.FileAddressList(fileID_NF_80mm_gtable,1)
Irmax_Ii_80mm = np.divide(gtable_80mm[0][:,7], gtable_80mm[0][:,7].max())
NF_80mm_exp = sio.loadmat(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Near Field Sims\100gPE4Sphere_80mm") 
theta_80mm_mesh = np.rad2deg(np.arctan(gtable_80mm[0][:,2]/0.08))
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

NF_380mm_exp = sio.loadmat(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Near Field Sims\100gPE4Sphere_380mm")
MxI_1_380mm = np.transpose(NF_380mm_exp['MxI'][:,:,0])
MxI_2_380mm = np.transpose(NF_380mm_exp['MxI'][:,:,1])
MxI_3_380mm = np.transpose(NF_380mm_exp['MxI'][:,:,2])
Mx_mean_380mm = np.transpose(NF_380mm_exp['MEANI'])
fileID_NF_380mm_gtable = r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Near Field Sims\Sims\Latest\380mm\*gtable" 
gtable380mm = pre.FileAddressList(fileID_NF_380mm_gtable,1)
#------------------------------------------------------------------------------


#Lower bound mesh sensitivity z = 0.055 ---------------------------------------
#Import Apollo data for mesh strategy
Apollo_FileList_z0_055 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_055\*.txt")
Apollo_gtable_z0_055 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_055\*gtable",1)
Apollo_gauges_z0_055 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_055\*gauges",1)
Apollo_log_z0_055 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\z0_055\*_log")
ult_cell_sizes_z0_055=[]
peak_impulse_z0_055 = []
CPU_times_z0_055=[]

total_mesh_vol = 0.35* 0.35 * 0.35

for i in enumerate(Apollo_FileList_z0_055):
    ult_cell_sizes_z0_055.append(pre.ElementSizes(Apollo_FileList_z0_055[i[0]]))
    CPU_times_z0_055.append(pre.CPUFinder(Apollo_log_z0_055[i[0]]))
    peak_impulse_z0_055.append(max(Apollo_gtable_z0_055[i[0]][:,7]))
    
ult_cell_sizes_z0_055 = np.asarray(ult_cell_sizes_z0_055)
R_z0_055 = np.divide(0.0255+cr, ult_cell_sizes_z0_055) #No of cells from gauge to centre of charge.
CPU_times_z0_055 = np.asarray(CPU_times_z0_055)    
peak_impulse_z0_055 = np.asarray(peak_impulse_z0_055)
no_cells = np.divide(total_mesh_vol, ult_cell_sizes_z0_055**3)
dist = 0.5/np.tan(np.deg2rad(80))
mesh_sensitivity_I = [Impulse_CFD(Apollo_gtable_z0_055[i][:,7], dist, 80, np.linspace(0,80,200)) for i in range(len(Apollo_gtable_z0_055))]
mesh_sensitivity_I = np.asarray(mesh_sensitivity_I)
#Indexes
start_100s = 0
final_100s = 4 + 1
final_200s = 8 + 1


# #Graph 1 ---------------------------------------------------------------------
# fig, [ax0, ax0a, ax1] = plt.subplots(1,3)
# fig.set_size_inches(7, 2.5)


# ax0.scatter(R_z0_055[start_100s:final_100s], peak_impulse_z0_055[start_100s:final_100s]/1e3, c = 'b', marker="s", s=10, label = 'Zone length 0.05m')
# ax0.scatter(R_z0_055[final_100s:final_200s], peak_impulse_z0_055[final_100s:final_200s]/1e3, c = 'gray',  marker="D", s=10, label = 'Zone length 0.02m')
# ax0.scatter(R_z0_055[final_200s::], peak_impulse_z0_055[final_200s::]/1e3, c = 'r', marker="o",  s=10, label = 'Zone length 0.01m')

# temp = max(peak_impulse_z0_055[final_200s::]/1e3)
# ax0.plot([0,50], [temp, temp], linewidth = 0.5, c = 'k')
# ax0.plot([0,50], [temp*0.9, temp*0.9], linewidth = 0.5, linestyle = '--', c = 'k',label = '$10\%$ convergence')

# ax0.set_xlabel('S / cell length')
# ax0.set_ylabel('Peak specific impulse (MPa.ms)', fontsize = 'x-small')
# ax0.set_ylim(4,12)
# ax0.set_xlim(0,41)
# ax0.minorticks_on()
# ax0.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
# ax0.grid(which='minor', ls=':', dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
# handles, labels = ax0.get_legend_handles_labels()
# ax0.legend(handles, labels, bbox_to_anchor=(0., 1.08, 3.85, .102), loc='lower left', ncol = 4, mode = "expand", borderaxespad=0., prop={'size':6})



# ax0a.scatter(R_z0_055[start_100s:final_100s], mesh_sensitivity_I[start_100s:final_100s]/1e3, c = 'b',  marker="s", s=10, label = 'Zone length 0.05m')
# ax0a.scatter(R_z0_055[final_100s:final_200s], mesh_sensitivity_I[final_100s:final_200s]/1e3, c = 'gray',  marker="D", s=10, label = 'Zone length 0.02m')
# ax0a.scatter(R_z0_055[final_200s::], mesh_sensitivity_I[final_200s::]/1e3, c = 'r', marker="o",  s=10, label = 'Zone length 0.01m')
# temp = max(mesh_sensitivity_I[final_200s::]/1e3)
# ax0a.plot([0,50], [temp, temp], linewidth = 0.5, c = 'k')
# ax0a.plot([0,50], [temp*0.9, temp*0.9], linewidth = 0.5, linestyle = '--', c = 'k',label = '$10\%$ convergence')
# ax0a.set_xlabel('S / cell length')
# ax0a.set_ylabel('Area integrated impulse (MN.ms)', fontsize = 'x-small')
# ax0a.set_ylim(250/1e3,450/1e3)
# ax0a.set_xlim(0,41)
# ax0a.minorticks_on()
# ax0a.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
# ax0a.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)


# ax1.scatter(R_z0_055[start_100s:final_100s], CPU_times_z0_055[start_100s:final_100s], c = 'b',  marker="s", s=10, label = 'Zone length 0.05m')
# ax1.scatter(R_z0_055[final_100s:final_200s], CPU_times_z0_055[final_100s:final_200s], c = 'gray',  marker="D", s=10, label = 'Zone length 0.02m')
# ax1.scatter(R_z0_055[final_200s::], CPU_times_z0_055[final_200s::], c = 'r', marker="o",  s=10, label = 'Zone Length 0.01m')
# ax1.set_xlabel('S / cell length')
# ax1.set_ylabel('Wall time (s)')
# ax1.minorticks_on()
# ax1.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
# ax1.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
# ax1.set_yscale('log')
# ax1.set_xlim(0,41)
# ax1.set_ylim(1,1e5)

# plt.tight_layout()
# ax0.locator_params(axis = 'x',tight=True, nbins=6)
# ax0a.locator_params(axis = 'x',tight=True, nbins=6)
# ax1.locator_params(axis = 'x',tight=True, nbins=6)
# ax0.locator_params(axis = 'y',tight=True, nbins=6)
# ax0a.locator_params(axis = 'y',tight=True, nbins=6)
# #ax1.locator_params(axis = 'y', tight=True, nbins=6)
# fig.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_convergence_z0_055_1.pdf', format = 'pdf')

#Graph 1a ---------------------------------------------------------------------
fig, [ax0, ax0a, ax1] = plt.subplots(1,3)
fig.set_size_inches(7, 2.5)


ax0.scatter(R_z0_055[start_100s:final_100s], peak_impulse_z0_055[start_100s:final_100s]/1e3, c = 'b', marker="s", edgecolor = 'k', s=10, label = 'Zone length 0.05m')
ax0.scatter(R_z0_055[final_100s:final_200s], peak_impulse_z0_055[final_100s:final_200s]/1e3, c = 'gray',  marker="D", edgecolor = 'k', s=10, label = 'Zone length 0.02m')
ax0.scatter(R_z0_055[final_200s::], peak_impulse_z0_055[final_200s::]/1e3, c = 'r', marker="o", edgecolor = 'k', s=10, label = 'Zone length 0.01m')

temp = max(peak_impulse_z0_055[final_200s::]/1e3)
ax0.plot([0,50], [temp, temp], linewidth = 0.5, c = 'k')
ax0.plot([0,50], [temp*0.9, temp*0.9], linewidth = 0.5, linestyle = '--', c = 'k',label = '$10\%$ convergence')

ax0.set_xlabel('S / cell length')
ax0.set_ylabel('Peak specific impulse (MPa.ms)', fontsize = 'x-small')
ax0.set_ylim(4,12)
ax0.set_xlim(0,41)
ax0.minorticks_on()
ax0.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax0.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, bbox_to_anchor=(0., 1.08, 3.85, .102), loc='lower left', ncol = 4, mode = "expand", borderaxespad=0., prop={'size':6})



ax0a.scatter(R_z0_055[start_100s:final_100s], mesh_sensitivity_I[start_100s:final_100s]/1e3, c = 'b',  marker="s", edgecolor = 'k', s=10, label = 'Zone length 0.05m')
ax0a.scatter(R_z0_055[final_100s:final_200s], mesh_sensitivity_I[final_100s:final_200s]/1e3, c = 'gray',  marker="D", edgecolor = 'k', s=10, label = 'Zone length 0.02m')
ax0a.scatter(R_z0_055[final_200s::], mesh_sensitivity_I[final_200s::]/1e3, c = 'r', marker="o", edgecolor = 'k', s=10, label = 'Zone length 0.01m')
temp = max(mesh_sensitivity_I[final_200s::]/1e3)
ax0a.plot([0,50], [temp, temp], linewidth = 0.5, c = 'k')
ax0a.plot([0,50], [temp*0.9, temp*0.9], linewidth = 0.5, linestyle = '--', c = 'k',label = '$10\%$ convergence')
ax0a.set_xlabel('S / cell length')
ax0a.set_ylabel('Area integrated impulse (MN.ms)', fontsize = 'x-small')
ax0a.set_ylim(250/1e3,450/1e3)
ax0a.set_xlim(0,41)
ax0a.minorticks_on()
ax0a.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax0a.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)


ax1.scatter(R_z0_055[start_100s:final_100s], CPU_times_z0_055[start_100s:final_100s], c = 'b',  marker="s",edgecolor = 'k', s=10, label = 'Zone length 0.05m')
ax1.scatter(R_z0_055[final_100s:final_200s], CPU_times_z0_055[final_100s:final_200s], c = 'gray',  marker="D",edgecolor = 'k', s=10, label = 'Zone length 0.02m')
ax1.scatter(R_z0_055[final_200s::], CPU_times_z0_055[final_200s::], c = 'r', marker="o", edgecolor = 'k', s=10, label = 'Zone Length 0.01m')
ax1.set_xlabel('S / cell length')
ax1.set_ylabel('Wall time (s)')
ax1.minorticks_on()
ax1.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax1.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
ax1.set_yscale('log')
ax1.set_xlim(0,41)
ax1.set_ylim(1,1e5)

plt.tight_layout()
ax0.locator_params(axis = 'x',tight=True, nbins=6)
ax0a.locator_params(axis = 'x',tight=True, nbins=6)
ax1.locator_params(axis = 'x',tight=True, nbins=6)
ax0.locator_params(axis = 'y',tight=True, nbins=6)
ax0a.locator_params(axis = 'y',tight=True, nbins=6)
#ax1.locator_params(axis = 'y', tight=True, nbins=6)
fig.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_convergence_z0_055_1a.pdf', format = 'pdf')


# #Graph 1a ---------------------------------------------------------------------
# fig, [ax0, ax0a, ax1] = plt.subplots(1,3)
# fig.set_size_inches(7, 2.5)


# ax0.scatter(R_z0_055[start_100s:final_100s], peak_impulse_z0_055[start_100s:final_100s]/1e3, marker="s", facecolors = 'None', edgecolor = 'b', s=10, label = 'Zone length 0.05m')
# ax0.scatter(R_z0_055[final_100s:final_200s], peak_impulse_z0_055[final_100s:final_200s]/1e3, marker="D", edgecolor = 'gray', facecolor = 'None', s=10, label = 'Zone length 0.02m')
# ax0.scatter(R_z0_055[final_200s::], peak_impulse_z0_055[final_200s::]/1e3, marker="o", edgecolor = 'r', facecolors = 'None',s=10, label = 'Zone length 0.01m')

# temp = max(peak_impulse_z0_055[final_200s::]/1e3)
# ax0.plot([0,50], [temp, temp], linewidth = 0.5, c = 'k')
# ax0.plot([0,50], [temp*0.9, temp*0.9], linewidth = 0.5, linestyle = '--', c = 'k',label = '$10\%$ convergence')

# ax0.set_xlabel('S / cell length')
# ax0.set_ylabel('Peak specific impulse (MPa.ms)', fontsize = 'x-small')
# ax0.set_ylim(4,12)
# ax0.set_xlim(0,41)
# ax0.minorticks_on()
# ax0.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
# ax0.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
# handles, labels = ax0.get_legend_handles_labels()
# ax0.legend(handles, labels, bbox_to_anchor=(0., 1.08, 3.85, .102), loc='lower left', ncol = 4, mode = "expand", borderaxespad=0., prop={'size':6})



# ax0a.scatter(R_z0_055[start_100s:final_100s], mesh_sensitivity_I[start_100s:final_100s]/1e3, marker="s", facecolor = 'None', edgecolor = 'b', s=10, label = 'Zone length 0.05m')
# ax0a.scatter(R_z0_055[final_100s:final_200s], mesh_sensitivity_I[final_100s:final_200s]/1e3, marker="D", edgecolor = 'gray', facecolor = 'None', s=10, label = 'Zone length 0.02m')
# ax0a.scatter(R_z0_055[final_200s::], mesh_sensitivity_I[final_200s::]/1e3, marker="o", edgecolor = 'r', facecolor = 'None', s=10, label = 'Zone length 0.01m')
# temp = max(mesh_sensitivity_I[final_200s::]/1e3)
# ax0a.plot([0,50], [temp, temp], linewidth = 0.5, c = 'k')
# ax0a.plot([0,50], [temp*0.9, temp*0.9], linewidth = 0.5, linestyle = '--', c = 'k',label = '$10\%$ convergence')
# ax0a.set_xlabel('S / cell length')
# ax0a.set_ylabel('Area integrated impulse (MN.ms)', fontsize = 'x-small')
# ax0a.set_ylim(250/1e3,450/1e3)
# ax0a.set_xlim(0,41)
# ax0a.minorticks_on()
# ax0a.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
# ax0a.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)


# ax1.scatter(R_z0_055[start_100s:final_100s], CPU_times_z0_055[start_100s:final_100s], marker="s", facecolor = 'None', edgecolor = 'b', s=10, label = 'Zone length 0.05m')
# ax1.scatter(R_z0_055[final_100s:final_200s], CPU_times_z0_055[final_100s:final_200s], marker="D", edgecolor = 'gray', facecolor = 'None', s=10, label = 'Zone length 0.02m')
# ax1.scatter(R_z0_055[final_200s::], CPU_times_z0_055[final_200s::], marker="o", edgecolor = 'r', facecolor = 'None', s=10, label = 'Zone Length 0.01m')
# ax1.set_xlabel('S / cell length')
# ax1.set_ylabel('Wall time (s)')
# ax1.minorticks_on()
# ax1.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
# ax1.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
# ax1.set_yscale('log')
# ax1.set_xlim(0,41)
# ax1.set_ylim(1,1e5)

# plt.tight_layout()
# ax0.locator_params(axis = 'x',tight=True, nbins=6)
# ax0a.locator_params(axis = 'x',tight=True, nbins=6)
# ax1.locator_params(axis = 'x',tight=True, nbins=6)
# ax0.locator_params(axis = 'y',tight=True, nbins=6)
# ax0a.locator_params(axis = 'y',tight=True, nbins=6)
# #ax1.locator_params(axis = 'y', tight=True, nbins=6)
# fig.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_convergence_z0_055_1b.pdf', format = 'pdf')


















#Dataset sensitivity
#filelists
latest_1500_r3_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res3\*.txt")
latest_1500_r4_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res4\*.txt")
latest_1500_r5_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res5\*.txt")
latest_var_r4_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\VAR_ZL40mm_res4\*.txt")
testing_DMA_file = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\testing_DMA\*.txt")
#Checking Z
latest_1500_r3_z = [(pre.standoff_func(latest_1500_r3_file[i]))/(charge_mass**(1/3)) for i in range(len(latest_1500_r3_file))]
latest_1500_r4_z = [(pre.standoff_func(latest_1500_r4_file[i]))/(charge_mass**(1/3)) for i in range(len(latest_1500_r4_file))]
latest_1500_r5_z = [(pre.standoff_func(latest_1500_r5_file[i]))/(charge_mass**(1/3)) for i in range(len(latest_1500_r5_file))]
latest_1500_r5_z_centre = [(pre.standoff_func(latest_1500_r5_file[i]))/(charge_mass**(1/3)) for i in range(len(latest_1500_r5_file))]
latest_var_r4_z = [(pre.standoff_func(latest_var_r4_file[i]))/(charge_mass**(1/3)) for i in range(len(latest_var_r4_file))]
testing_DMA_z = [(pre.standoff_func(testing_DMA_file[i]))/(charge_mass**(1/3)) for i in range(len(testing_DMA_file))]
#Gtables
latest_1500_r3 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res3\*gtable",1)
latest_1500_r4 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res4\*gtable",1)
latest_1500_r5 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res5\*gtable",1)
latest_var_r4 = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\VAR_ZL40mm_res4\*gtable",1)
testing_DMA = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\testing_DMA\*gtable", 1)
latest_1500_r3 = np.asarray([latest_1500_r3[i][:,7] for i in range(len(latest_1500_r3))]).T
latest_1500_r4 = np.asarray([latest_1500_r4[i][:,7] for i in range(len(latest_1500_r4))]).T
latest_1500_r5 = np.asarray([latest_1500_r5[i][:,7] for i in range(len(latest_1500_r5))]).T
latest_var_r4 = np.asarray([latest_var_r4[i][:,7] for i in range(len(latest_var_r4))]).T
testing_DMA = np.asarray([testing_DMA[i][:,7] for i in range(len(testing_DMA))]).T

#Wall times
latest_1500_r3_t = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res3\*_log")
latest_1500_r4_t = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res4\*_log")
latest_1500_r5_t = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res5\*_log")
latest_var_r4_t = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\VAR_ZL40mm_res4\*_log")
testing_DMA_t = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\testing_DMA\*_log")
latest_1500_r3_t = [pre.CPUFinder(i) for i in latest_1500_r3_t]
latest_1500_r4_t = [pre.CPUFinder(i) for i in latest_1500_r4_t]
latest_1500_r5_t = [pre.CPUFinder(i) for i in latest_1500_r5_t]
latest_var_r4_t = [pre.CPUFinder(i) for i in latest_var_r4_t]
testing_DMA_t = [pre.CPUFinder(i) for i in testing_DMA_t]

#Gauges
latest_1500_r3_gauges = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res3\*gauges",1)
latest_1500_r4_gauges = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res4\*gauges",1)
latest_1500_r5_gauges = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\1500mm_ZL100mm_res5\*gauges",1)
latest_var_r4_gauges = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\VAR_ZL40mm_res4\*gauges",1)
testing_DMA_gauges = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16_latest\testing_DMA\*gauges", 1)

dist = 0.5/np.tan(np.deg2rad(80))
latest_1500_r3_I = [Impulse_CFD(latest_1500_r3[:,i], dist, 80, np.linspace(0,80,200)) for i in range(len(latest_1500_r3_file))]
latest_1500_r4_I = [Impulse_CFD(latest_1500_r4[:,i], dist, 80, np.linspace(0,80,200)) for i in range(len(latest_1500_r4_file))]
latest_1500_r5_I = [Impulse_CFD(latest_1500_r5[:,i], dist, 80, np.linspace(0,80,200)) for i in range(len(latest_1500_r5_file))]
latest_var_r4_I = [Impulse_CFD(latest_var_r4[:,i], dist, 80, np.linspace(0,80,200)) for i in range(len(latest_var_r4_file))]
testing_DMA_I = [Impulse_CFD(testing_DMA[:,i], dist, 80, np.linspace(0,80,200)) for i in range(len(testing_DMA_file))]

fig_dataset2, ax0 = plt.subplots(1,1)
fig_dataset2.set_size_inches(2.5,2.5)
ax0.plot(np.linspace(0,80,200),latest_var_r4[:,3]/1e3, 'c', ls = '-', lw = 1)
ax0.plot(np.linspace(0,80,200),latest_1500_r3[:,3]/1e3, 'g', ls= '--', lw = 1)
ax0.plot(np.linspace(0,80,200),latest_1500_r4[:,3]/1e3, 'm', ls = '--', lw = 0.5)
ax0.plot(np.linspace(0,80,200),latest_1500_r5[:,3]/1e3, 'r', ls = '-', lw = 0.5)
ax0.plot(np.linspace(0,80,200),testing_DMA[:,0]/1e3, 'k', ls = '-.', lw = 0.75)
ax0.plot(np.linspace(0,80,200),testing_DMA[:,1]/1e3, 'b', ls = '-.', lw = 0.75)
labels = ['A', 'B', 'C', 'D', 'E', 'F']
colors = ['c', 'g', 'm', 'r', 'k', 'b']
lws = [1, 1, 0.5, 0.5, 0.75, 0.75]
lss = ['-', '--', '--', '-', '-.', '-.']
lines = [Line2D([0], [0],  lw=lws[i], ls = lss[i], color=colors[i]) for i in range(len(labels))]
ax0.legend(lines,labels, loc='upper right', prop={'size':6})
ax0.set_ylabel('Peak specific impulse (MPa.ms)')
ax0.set_xlabel('Angle of incidence')
ax0.set_xlim(0,80)
ax0.set_ylim(0,5)
ax0.minorticks_on()
ax0.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax0.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
plt.tight_layout()
fig_dataset2.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_convergence_dataset2.pdf', format = 'pdf')

fig_dataset3, ax1 = plt.subplots(1,1)
fig_dataset3.set_size_inches(2.5,2.5)
ax1.scatter(latest_var_r4_z, np.asarray(latest_var_r4_I)/1e3, marker ="o", facecolors = 'none', edgecolors='c',s=10., label = 'A')
ax1.scatter(latest_1500_r3_z, np.asarray(latest_1500_r3_I)/1e3, marker = "D", facecolors = 'none', edgecolors='g', s=10., label = 'B')
ax1.scatter(latest_1500_r4_z, np.asarray(latest_1500_r4_I)/1e3, marker = "^",facecolors = 'none', edgecolors='m', s=10., label = 'C')
ax1.scatter(latest_1500_r5_z, np.asarray(latest_1500_r5_I)/1e3, marker = "s", facecolors = 'none', edgecolors='r', s=10., label = 'D')
ax1.scatter(testing_DMA_z[0], np.asarray(testing_DMA_I[0])/1e3, marker = '>', facecolors = 'none', edgecolors='k', s=10., label = 'E')
ax1.scatter(testing_DMA_z[1], np.asarray(testing_DMA_I[1])/1e3, marker = '<',facecolors = 'none', edgecolors='b', s=10., label = 'F')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='upper right', prop={'size':6})
ax1.set_xlabel('Scaled distance, Z ' + r'$(m/kg^{1/3}$)')
ax1.set_ylabel('Area integrated impulse (MPa.ms)', fontsize = 'x-small')
ax1.set_ylim(0.2,0.5)
#ax1.set_title(r'$1m^2$ area integrated impulse')
ax1.minorticks_on()
ax1.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax1.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
plt.tight_layout()
fig_dataset3.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_convergence_dataset3.pdf', format = 'pdf')

fig_dataset3a, ax1 = plt.subplots(1,1)
fig_dataset3a.set_size_inches(2.5,2.5)
ax1.scatter(latest_var_r4_z, latest_var_r4_t, marker ="o", facecolors = 'none', edgecolors='c',s=10., label = 'A')
ax1.scatter(latest_1500_r3_z, latest_1500_r3_t , marker = "D", facecolors = 'none', edgecolors='g', s=10., label = 'B')
ax1.scatter(latest_1500_r4_z, latest_1500_r4_t, marker = "^",facecolors = 'none', edgecolors='m', s=10., label = 'C')
ax1.scatter(latest_1500_r5_z, latest_1500_r5_t, marker = "s", facecolors = 'none', edgecolors='r', s=10., label = 'D')
ax1.scatter(testing_DMA_z[0], testing_DMA_t[0], marker = '>', facecolors = 'none', edgecolors='k', s=10., label = 'E')
ax1.scatter(testing_DMA_z[1], testing_DMA_t[1], marker = '<',facecolors = 'none', edgecolors='b', s=10., label = 'F')
ax1.set_xlabel('Scaled distance, Z ' + r'$(m/kg^{1/3}$)')
ax1.set_ylabel('Wall time (s)')
ax1.minorticks_on()
ax1.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax1.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
ax1.set_yscale('log')
plt.tight_layout()
fig_dataset3a.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_convergence_dataset3a.pdf', format = 'pdf')


# fig_dataset1, ax = plt.subplots(1,1)
# fig_dataset1.set_size_inches(2.5,2.5)
# ax.plot(np.linspace(0,80,200),latest_1500_r5/1e3, 'r', ls = '-', lw = 0.5)
# ax.set_ylabel('Peak specific impulse (MPa.ms)')
# ax.set_xlabel('Angle of incidence')
# labels = [str(round(latest_1500_r5_z_centre[i], 3)) for i in range(len(latest_1500_r5_file))]
# colors = ['r' for i in range(len(latest_1500_r5_file)) ]
# lws = [0.5 for i in range(len(latest_1500_r5_file))]
# lss = ['-' for i in range(len(latest_1500_r5_file))]
# lines = [Line2D([0], [0],  lw=lws[i], ls = lss[i], color=colors[i]) for i in range(len(labels))]
# ax.legend(lines,labels, loc='upper right', title=r'$z(m/kg^{1/3})$', title_fontsize = 'x-small', prop={'size':6})
# ax.minorticks_on()
# ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
# ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
# ax.set_xlim(0,80)
# ax.set_ylim(0,12)
# plt.tight_layout()
# fig_dataset1.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_convergence_dataset1.pdf', format = 'pdf')

# fig_dataset1, ax = plt.subplots(1,1)
# fig_dataset1.set_size_inches(2.5,2.5)
# smooth = np.asarray([savgol_filter(latest_1500_r5[:,i], 101, 3) for i in range(len(latest_1500_r5_file))]).T
# ax.plot(np.linspace(0,80,200), smooth/1e3, 'r', ls = '-', lw = 0.5)
# ax.set_ylabel('Peak specific impulse (MPa.ms)')
# ax.set_xlabel('Angle of incidence')
# labels = [str(round(latest_1500_r5_z_centre[i], 3)) for i in range(len(latest_1500_r5_file))]
# colors = ['r' for i in range(len(latest_1500_r5_file)) ]
# lws = [0.5 for i in range(len(latest_1500_r5_file))]
# lss = ['-' for i in range(len(latest_1500_r5_file))]
# lines = [Line2D([0], [0],  lw=lws[i], ls = lss[i], color=colors[i]) for i in range(len(labels))]
# ax.legend(lines,labels, loc='upper right', title=r'$z(m/kg^{1/3})$', title_fontsize = 'x-small', prop={'size':6})
# ax.minorticks_on()
# ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
# ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
# ax.set_xlim(0,80)
# ax.set_ylim(0,12)
# plt.tight_layout()
# fig_dataset1.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_convergence_dataset1a.pdf', format = 'pdf')

# fig_dataset1, ax = plt.subplots(1,1)
# fig_dataset1.set_size_inches(2.5,2.5)
# ax.plot(np.linspace(0,80,200), smooth/smooth.max(0), 'r', ls = '-', lw = 0.5)
# ax.scatter(theta_exp_80mm, np.divide(MxI_1_80mm, max(MxI_1_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = '80mm exp.')
# ax.scatter(theta_exp_80mm, np.divide(MxI_2_80mm, max(MxI_2_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
# ax.scatter(theta_exp_80mm, np.divide(MxI_3_80mm, max(MxI_3_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
# ax.scatter(theta_exp_80mm_mean, np.divide(Mx_mean_80mm, max(Mx_mean_80mm)), marker="o", s=15., label = '80mm Exp Mean')
# ax.set_ylabel('Peak specific impulse ratio')
# ax.set_xlabel('Angle of incidence')
# legend_elements = [Line2D([], [],  lw=0.5, linestyle = '-', color='r', label = 'CFD'), 
#                    Line2D([], [],  color=[0.75,0.75,0.75], marker ='x', linestyle='None', markersize=4, label='Exp. repeats'),
#                    Line2D([], [],   marker ='o', linestyle='None', markersize=4, label='Exp. mean')]
# ax.legend(handles=legend_elements, loc='upper right', prop={'size':6})
# ax.minorticks_on()
# ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
# ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25)
# ax.set_xlim(0,80)
# ax.set_ylim(0,1)
# plt.tight_layout()
#fig_dataset1.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_convergence_dataset1b.pdf', format = 'pdf')



#Graph 3 ---------------------------------------------------------------------
#NF validation graph
# #theta vs scaled I plots 
Apollo_gtable_z80mm_chosenmesh = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\80mm_validation\*gtable",1)
Apollo_gtable_z80mm_first = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Near Field Sims\Sims\Latest\80mm_with_afterburn\*gtable",1)
z80mm_first_theta = np.rad2deg(np.arctan2(Apollo_gtable_z80mm_first[0][:,2], 0.08))
Apollo_gauges_z80mm_chosenmesh = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\80mm_validation\*gauges",1)
Apollo_gauges_z80mm_first = pre.FileAddressList(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Near Field Sims\Sims\Latest\80mm_with_afterburn\*gauges",1)
#theta vs scaled I plots 
theta = np.linspace(0,80,200)
term = 0.15e-3
z80mm_chosenmesh_fin = int(np.argwhere(Apollo_gauges_z80mm_chosenmesh[0][:,0]>term)[0][0])
z80mm_chosenmesh_adjustedi = np.max(Apollo_gauges_z80mm_chosenmesh[0][0:z80mm_chosenmesh_fin,201:], axis = 0)




# fig2, ax = plt.subplots(1,1)
# fig2.set_size_inches(2.5,2.5)
# ax.scatter(theta_exp_80mm, MxI_1_80mm/1e3, marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = 'Exp.')
# ax.scatter(theta_exp_80mm, MxI_2_80mm/1e3, marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
# ax.scatter(theta_exp_80mm, MxI_3_80mm/1e3, marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
# ax.scatter(theta_exp_80mm_mean, Mx_mean_80mm/1e3, marker="o", s=15., label = 'Exp. mean')
# ax.plot(theta, savgol_filter(z80mm_chosenmesh_adjustedi/1e3,101,3), 'r', ls = '-', lw = 0.5, label = 'CFD')
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels, loc='upper right', prop={'size':6})
# ax.set_ylim(0,5)
# ax.set_xlim(0,80)
# ax.set_xlabel('theta (degrees)')
# ax.set_ylabel('Peak specific impulse (MPa.ms)')
# ax.minorticks_on()
# ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
# ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 
# plt.tight_layout()
#fig2.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_convergence_z0_055_3.pdf', format = 'pdf')



fig2, ax = plt.subplots(1,1)
fig2.set_size_inches(2.5,2.5)
ax.scatter(abs(coords)*1e3, MxI_1_80mm/1e3, marker="x", s=10., color='grey', edgecolors='none', label = 'Exp.')
ax.scatter(abs(coords)*1e3, MxI_2_80mm/1e3, marker="x", s=10., color='grey', edgecolors='none')
ax.scatter(abs(coords)*1e3, MxI_3_80mm/1e3, marker="x", s=10., color='grey', edgecolors='none')
ax.scatter(np.arange(0,101,25), Mx_mean_80mm/1e3, marker="o", s=10., edgecolor = 'k', label = 'Exp. mean')
ax.plot(0.08*np.tan(np.deg2rad(theta))*1e3, savgol_filter(z80mm_chosenmesh_adjustedi/1e3,101,3), 'r', ls = '-', lw = 1, label = 'CFD')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper right', prop={'size':6})
ax.set_ylim(0,5)
ax.set_xlim(0,100)
ax.set_xlabel('Distance from centre (mm)')
ax.set_ylabel('Peak specific impulse (MPa.ms)')
ax.minorticks_on()
ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 
plt.tight_layout()
fig2.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_80mm_val.pdf', format = 'pdf')
fig2, ax = plt.subplots(1,1)
fig2.set_size_inches(2.5,2.5)
ax.scatter(abs(coords)*1e3, MxI_1_380mm/1e3, marker="x", s=10., color='grey', edgecolors='none', label = 'Exp.')
ax.scatter(abs(coords)*1e3, MxI_2_380mm/1e3, marker="x", s=10., color='grey', edgecolors='none')
ax.scatter(abs(coords)*1e3, MxI_3_380mm/1e3, marker="x", s=10., color='grey', edgecolors='none')
ax.scatter(np.arange(0,101,25), Mx_mean_380mm/1e3, marker="o", s=10., edgecolor = 'k', label = 'Exp. mean')
ax.plot(gtable380mm[0][:,2]*1e3,gtable380mm[0][:,7]/1e3, 'r', ls = '-', lw = 1, label = 'CFD')
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper right', prop={'size':6})
ax.set_ylim(0,0.5)
ax.set_xlim(0,100)
ax.set_xlabel('Distance from centre (mm)')
ax.set_ylabel('Peak specific impulse (MPa.ms)')
ax.minorticks_on()
ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 
plt.tight_layout()
fig2.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\mesh_380mm_val.pdf', format = 'pdf')


# #-----------------------------------------------------------------------------



#specific impulse validation plots
#theta = 0 
fig7, [ax1, ax] = plt.subplots(1,2)
fig7.set_size_inches(5, 1.8)
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,4,0]/1e3, 'grey',  label = 'Exp.')
ax.plot(Apollo_gauges_z80mm_chosenmesh[0][:,0]*1000, Apollo_gauges_z80mm_chosenmesh[0][:,201]/1e3, 'r', lw = 1, label = 'CFD')

ax.set_xlim(0,0.15)
ax.set_ylim(0,5)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Impulse (MPa.ms)')
ax.minorticks_on()
ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 

ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,4,0]/1e6, 'grey', label = 'Exp: theta 0')
ax1.plot(Apollo_gauges_z80mm_chosenmesh[0][:,0]*1000, Apollo_gauges_z80mm_chosenmesh[0][:,1]/1e6, 'r', lw = 1)

ax1.set_xlim(0,0.15)
ax1.set_ylim(-50,250)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Overpressure (MPa)')
plt.tight_layout()
ax.xaxis.set_major_locator(LinearLocator(4)) 
ax.yaxis.set_major_locator(LinearLocator(6))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.xaxis.set_major_locator(LinearLocator(4)) 
ax1.yaxis.set_major_locator(LinearLocator(4))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
handles, labels = ax.get_legend_handles_labels()
ax1.minorticks_on()
ax1.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax1.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 
ax.legend(handles, labels, loc='lower right', prop={'size':6})
fig7.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\80mm_validation.pdf', format = 'pdf')

#theta = 17
fig7a, [ax1, ax] = plt.subplots(1,2)
fig7a.set_size_inches(5, 1.8)
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,3,0]/1e3, 'grey')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,5,0]/1e3, 'grey')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,12,0]/1e3, 'grey')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,13,0]/1e3, 'grey')
#ax.plot(Apollo_gauges_z80mm_first[0][:,0]*1000, Apollo_gauges_z80mm_first[0][:,125]/1e3, 'k', lw = 0.75)
ax.plot(Apollo_gauges_z80mm_chosenmesh[0][:,0]*1000, Apollo_gauges_z80mm_chosenmesh[0][:,243]/1e3, 'r', lw = 1)

ax.set_xlim(0,0.15)
ax.set_ylim(0,5)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Impulse (MPa.ms)')
ax.minorticks_on()
ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 

ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,3,0]/1e6, 'grey')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,5,0]/1e6, 'grey')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,12,0]/1e6, 'grey')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,13,0]/1e6, 'grey')
#ax1.plot(Apollo_gauges_z80mm_first[0][:,0]*1000, Apollo_gauges_z80mm_first[0][:,25]/1e6, 'k', lw = 0.75)
ax1.plot(Apollo_gauges_z80mm_chosenmesh[0][:,0]*1000, Apollo_gauges_z80mm_chosenmesh[0][:,43]/1e6, 'r', lw = 1)

ax1.set_xlim(0,0.15)
ax1.set_ylim(-50,250)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Overpressure (MPa)')
plt.tight_layout()
ax.xaxis.set_major_locator(LinearLocator(4)) 
ax.yaxis.set_major_locator(LinearLocator(6))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.xaxis.set_major_locator(LinearLocator(4)) 
ax1.yaxis.set_major_locator(LinearLocator(4))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.minorticks_on()
ax1.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax1.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 
fig7a.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\80mm_validation_a.pdf', format = 'pdf')

#theta = 32
fig7b, [ax1, ax] = plt.subplots(1,2)
fig7b.set_size_inches(5, 1.8)
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,2,0]/1e3, 'grey', label = 'Exp: theta 0')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,6,0]/1e3, 'grey')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,11,0]/1e3, 'grey')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,14,0]/1e3, 'grey')
#ax.plot(Apollo_gauges_z80mm_first[0][:,0]*1000, Apollo_gauges_z80mm_first[0][:,150]/1e3, 'k', lw = 0.75)
ax.plot(Apollo_gauges_z80mm_chosenmesh[0][:,0]*1000, Apollo_gauges_z80mm_chosenmesh[0][:,280]/1e3, 'r', lw = 1)

ax.set_xlim(0,0.15)
ax.set_ylim(0,3)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Impulse (MPa.ms)')
ax.minorticks_on()
ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 

ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,2,0]/1e6, 'grey', label = 'Exp: theta 0')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,6,0]/1e6, 'grey')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,11,0]/1e6, 'grey')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,14,0]/1e6, 'grey')
#ax1.plot(Apollo_gauges_z80mm_first[0][:,0]*1000, Apollo_gauges_z80mm_first[0][:,50]/1e6, 'k', lw = 0.75)
ax1.plot(Apollo_gauges_z80mm_chosenmesh[0][:,0]*1000, Apollo_gauges_z80mm_chosenmesh[0][:,80]/1e6, 'r', lw = 1)

ax1.set_xlim(0,0.15)
ax1.set_ylim(-50,150)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Overpressure (MPa)')
plt.tight_layout()
ax.xaxis.set_major_locator(LinearLocator(4)) 
ax.yaxis.set_major_locator(LinearLocator(4))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.xaxis.set_major_locator(LinearLocator(4)) 
ax1.yaxis.set_major_locator(LinearLocator(5))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.minorticks_on()
ax1.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax1.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 
fig7b.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\80mm_validation_b.pdf', format = 'pdf')

#theta = 43
fig7c, [ax1, ax] = plt.subplots(1,2)
fig7c.set_size_inches(5, 1.8)
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,1,0]/1e3, 'grey', label = 'Exp: theta 0')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,7,0]/1e3, 'grey')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,10,0]/1e3, 'grey')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,15,0]/1e3, 'grey')
#ax.plot(Apollo_gauges_z80mm_first[0][:,0]*1000, Apollo_gauges_z80mm_first[0][:,174]/1e3, 'k', lw = 0.75)
ax.plot(Apollo_gauges_z80mm_chosenmesh[0][:,0]*1000, Apollo_gauges_z80mm_chosenmesh[0][:,307]/1e3, 'r', lw = 1)

ax.set_xlim(0,0.15)
ax.set_ylim(0,2)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Impulse (MPa.ms)')
ax.minorticks_on()
ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 

ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,1,0]/1e6, 'grey', label = 'Exp: theta 0')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,7,0]/1e6, 'grey')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,10,0]/1e6, 'grey')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,15,0]/1e6, 'grey')
#ax1.plot(Apollo_gauges_z80mm_first[0][:,0]*1000, Apollo_gauges_z80mm_first[0][:,74]/1e6, 'k', lw = 0.75)
ax1.plot(Apollo_gauges_z80mm_chosenmesh[0][:,0]*1000, Apollo_gauges_z80mm_chosenmesh[0][:,107]/1e6, 'r', lw = 1)

ax1.set_xlim(0,0.15)
ax1.set_ylim(-25,100)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Overpressure (MPa)')
ax.minorticks_on()
ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 
plt.tight_layout()
ax.xaxis.set_major_locator(LinearLocator(4)) 
ax.yaxis.set_major_locator(LinearLocator(3))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.xaxis.set_major_locator(LinearLocator(4)) 
ax1.yaxis.set_major_locator(LinearLocator(6))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.minorticks_on()
ax1.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax1.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 
fig7c.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\80mm_validation_c.pdf', format = 'pdf')

#theta = 51
fig7d, [ax1, ax] = plt.subplots(1,2)
fig7d.set_size_inches(5, 1.8)
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,0,0]/1e3, 'grey', label = 'Exp: theta 0')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,8,0]/1e3, 'grey')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,9,0]/1e3, 'grey')
ax.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['impulse'][:,16,0]/1e3, 'grey')
#ax.plot(Apollo_gauges_z80mm_first[0][:,0]*1000, Apollo_gauges_z80mm_first[0][:,198]/1e3, 'k', lw = 0.75)
ax.plot(Apollo_gauges_z80mm_chosenmesh[0][:,0]*1000, Apollo_gauges_z80mm_chosenmesh[0][:,327]/1e3, 'r', lw = 1)

ax.set_xlim(0,0.15)
ax.set_ylim(0,1)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Impulse (MPa.ms)')
ax.minorticks_on()
ax.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 


ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,0,0]/1e6, 'grey', label = 'Exp: theta 0')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,8,0]/1e6, 'grey')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,9,0]/1e6, 'grey')
ax1.plot(NF_80mm_exp['time'][:,0]*1000, NF_80mm_exp['pressure'][:,16,0]/1e6, 'grey')
#ax1.plot(Apollo_gauges_z80mm_first[0][:,0]*1000, Apollo_gauges_z80mm_first[0][:,98]/1e6, 'k', lw = 0.75)
ax1.plot(Apollo_gauges_z80mm_chosenmesh[0][:,0]*1000, Apollo_gauges_z80mm_chosenmesh[0][:,127]/1e6, 'r', lw = 1)

ax1.set_xlim(0,0.15)
ax1.set_ylim(-25,100)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Overpressure (MPa)')
plt.tight_layout()
ax.xaxis.set_major_locator(LinearLocator(4)) 
ax.yaxis.set_major_locator(LinearLocator(3))
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.xaxis.set_major_locator(LinearLocator(4)) 
ax1.yaxis.set_major_locator(LinearLocator(6))
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.minorticks_on()
ax1.grid(which='major', ls = '-', color = [0.15, 0.15, 0.15], alpha=0.15)
ax1.grid(which='minor', ls=':',  dashes=(1,5,1,5), color = [0.1, 0.1, 0.1], alpha=0.25) 
fig7d.savefig(os.environ['USERPROFILE'] + r'\Dropbox\Papers\Paper1\Graphs\80mm_validation_d.pdf', format = 'pdf')

#------------------------------------------------------------------------------
