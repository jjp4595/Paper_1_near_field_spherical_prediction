"""
Peak impulse predictor for theta. 
"""

import preamble_functions as pre
import matplotlib.pyplot as plt #3.0.2
from matplotlib import cm
import numpy as np #1.15.4
from mpl_toolkits.mplot3d import Axes3D
import scipy.io as sio
import lmfit as lm
from sklearn.metrics import mean_squared_error
from MCEER_curves import MCEER

#plt.rcParams["font.family"] = "cmr10" #Set Graph fonts to cmr10
params = {'font.family':'serif',
        'axes.labelsize':'small',
        'xtick.labelsize':'x-small',
        'ytick.labelsize':'x-small', 
        'legend.fontsize':'small',
        'legend.title_fontsize':'small',
        'legend.fancybox': True,
        'legend.framealpha': 0.5,
        'legend.shadow': False,
        'legend.frameon': True,
        'grid.linestyle':'--',
        'grid.linewidth':'0.5',
        'lines.linewidth':'0.5'}
plt.rcParams.update(params)

#Import Apollo data
Apollo_FileList = pre.FileAddressList(r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*.txt")
Apollo_gtable = pre.FileAddressList(r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*gtable",1)
Apollo_gauges = pre.FileAddressList(r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*gauges",1)
#Apollo_FileList = pre.FileAddressList(r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_dataset\*.txt")
#Apollo_gtable = pre.FileAddressList(r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_dataset\*gtable",1)
#Apollo_gauges = pre.FileAddressList(r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_dataset\*gauges",1)


#Some charge properties
charge_rad = 0.0246
charge_mass = 0.1


#Finding theta and creating data structures
clear_standoff=np.zeros((len(Apollo_FileList),1))
peak_impulse=[]
Icr_Ir = []

theta = []
for i in range(len(Apollo_FileList)):
    clear_standoff[i] = pre.standoff_func(Apollo_FileList[i]) - charge_rad
    peak_impulse.append(Apollo_gtable[i][:,7])
    #incident_impulse.append(Apollo_gtable[i][incident_index::,7])
    #ref_factor.append(peak_impulse[i]/incident_impulse[i])
    Icr_Ir.append(peak_impulse[i] / max(peak_impulse[i]))
    theta.append(np.rad2deg(np.arctan(np.divide(Apollo_gtable[i][:,2], clear_standoff[i] + charge_rad))))

theta = np.stack(theta, axis=1)
peak_impulse = np.stack(peak_impulse, axis = 1)
Icr_Ir = np.stack(Icr_Ir, axis =1)
clear_standoff = np.transpose(np.repeat(clear_standoff, len(theta), axis=1))
clear_standoff = np.divide(clear_standoff, charge_rad)





#Load Data for 80mm and 380mm Apollo Experimental
fileID_NF_80mm_gtable = r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Near Field Sims\Sims\Latest\80mm_with_afterburn\*gtable"
gtable_80mm = pre.FileAddressList(fileID_NF_80mm_gtable,1)
Irmax_Ii_80mm = np.divide(gtable_80mm[0][:,7], gtable_80mm[0][:,7].max())
NF_80mm_exp = sio.loadmat(r"C:\Users\cip18jjp\Google Drive\Apollo Sims\Near Field Sims\100gPE4Sphere_80mm") 
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




#Graphs for paper -------------------------------------------------------------
#Graph 1 ----------------------------------------------------------------------
fig0, ax = plt.subplots(1,1)
fig0.set_size_inches(3, 2.5)
CS = ax.contourf(theta, clear_standoff, peak_impulse/1e3, levels = [0 , 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25], cmap = plt.cm.cubehelix)
cbar = fig0.colorbar(CS)
cbar.ax.invert_yaxis()
cbar.ax.set_ylabel('peak specific impulse (MPa.ms)')
ax.set_ylabel('standoff (clear charge radii)')
ax.set_xlabel('incident wave angle (degrees)')
plt.tight_layout()
fig0.savefig(r'C:\Users\cip18jjp\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\theta_peak_impulse_0.pdf', format = 'pdf')


#Graph 2 ----------------------------------------------------------------------
fig1, [ax0, ax1] = plt.subplots(1,2)
fig1.set_size_inches(3, 2.5)
ax0.scatter(clear_standoff[0,:], peak_impulse[0,:]/1e3, c = 'k', marker=".", s=10)
ax0.set_xlim(0, 2.5)
ax0.set_ylim(0, 30)
ax0.set_ylabel('peak specific impulse (MPa.ms)')
ax0.set_xlabel('standoff (clear charge radii)')
ax1.scatter(np.log10(clear_standoff[0,:]), np.log10(peak_impulse[0,:]/1e3), c = 'k', marker=".", s=10)
ax1.set_ylabel('peak specific impulse (MPa.ms)')
ax1.set_xlabel('standoff (clear charge radii)')
ax1.set_xlim(-0.6, 2)
ax1.set_ylim(-0.6, 2)
#------------------------------------------------------------------------------

#Henrych Model-- --------------------------------------------------------------
Qw = 4.79e6/4184 #Specific energy of PE4, 4.79 assumed from EOS.dat
Ux = (2 * Qw)**0.5
A0 = (Ux * 1.8) / (4*np.pi)
def Henrych(A0, W, a, theta):
    i = ((A0 * W) / (a**2)) * np.cos(np.deg2rad(theta))**4
    return i
def Henrych_I(A0, W, theta):
    return np.pi * A0 * W * np.sin(np.deg2rad(theta))**2


#MCEER information-------------------------------------------------------------
"""
Generating impulse from MCEER curves and modelling distribution via RPB model.
"""
TNTeq = 1
W = 0.1 * TNTeq
R = 0.08

def RPB_MCEER_i(W, R):
    #Wroot = W**(1/3)
    A = R * np.tan(np.deg2rad(80))
    B = A
    res = 100
    x = np.linspace(-A/2, A/2, res)
    y = np.linspace(-B/2, B/2, res)
    [X,Y] = np.meshgrid(x,y)
    Rs = (np.power(X,2) + np.power(Y,2) + R**2)**0.5
    #AOIS = np.arctan( np.divide((np.power(X,2) + np.power(Y,2))**0.5, R) )
    #Z = np.divide(Rs, Wroot)
    theta_MCEER = np.rad2deg(np.arccos(np.divide(R, Rs)))
    
    is_max = np.zeros((res,res))
    ir_max = np.zeros((res,res))
    for i in range(len(is_max)):
        for j in range(len(is_max)):
                is_max[i,j], ir_max[i,j] = MCEER(Rs[i,j], W)
    imp =( np.multiply(ir_max, np.divide(np.power(R,2), np.power(Rs,2)) ) + 
           np.multiply(is_max, (1 + np.divide(np.power(R,2), np.power(Rs,2)) - (2 * np.divide(R,Rs)) ))
          )/1e3
    Area = np.ones_like(imp)
    Area = Area * (x[1]-x[0]) * (y[1]-y[0])
    Area[:,0]  *= 0.5
    Area[0,:]  *= 0.5
    Area[:,-1] *= 0.5
    Area[-1,:] *= 0.5
    Imp = np.multiply(imp, Area)
    I = np.sum(Imp)
    return [X,Y, theta_MCEER, imp, I]
  
RPB_MCEER_exp = RPB_MCEER_i(0.1, 0.08)
      
fig2, ax = plt.subplots(1,1)
fig2.set_size_inches(3, 2.5)
CS = ax.contourf(RPB_MCEER_exp[0], RPB_MCEER_exp[1], RPB_MCEER_exp[3], cmap = plt.cm.cubehelix)
cbar = fig2.colorbar(CS)
cbar.ax.set_ylabel('peak specific impulse (MPa.ms)')
ax.set_ylabel('x-position')
ax.set_xlabel('y-position')
plt.tight_layout()

fig3, [ax0,ax1] = plt.subplots(1,2)
fig3.set_size_inches(5,2.5)
ax0.scatter(theta_exp_80mm, MxI_1_80mm/1e3, marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = 'Exp')
ax0.scatter(theta_exp_80mm, MxI_2_80mm/1e3, marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax0.scatter(theta_exp_80mm, MxI_3_80mm/1e3, marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax0.scatter(theta_exp_80mm_mean, Mx_mean_80mm/1e3, marker="o", s=15., label = 'Exp - mean')
ax0.plot(RPB_MCEER_exp[2][int(len(RPB_MCEER_exp[2])/2), 0:int(len(RPB_MCEER_exp[2])/2)], RPB_MCEER_exp[3][int(len(RPB_MCEER_exp[2])/2),0:int(len(RPB_MCEER_exp[2])/2)], 'k', label = 'MCEER/RPB')
ax0.plot(theta_80mm_mesh, gtable_80mm[0][:,7]/1e3, 'k-.', dashes=[12,6,12,6,3,6], label = 'CFD - 1.25mm')
ax0.plot(np.linspace(0,80, num=80), Henrych(255,0.1,0.08, np.linspace(0,80, num=80))/1e3, 'k:', label = 'Henrych')
ax0.set_xlabel('theta (degrees)')
ax0.set_ylabel('peak specific impulse (MPa.ms)')
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.72, 0.78), prop={'size':6})

ax1.set_xlabel('theta (degrees)')
ax1.set_ylabel('peak specific impulse ratio')
ax1.scatter(theta_exp_80mm, np.divide(MxI_1_80mm, max(MxI_1_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = '80mm Exp')
ax1.scatter(theta_exp_80mm, np.divide(MxI_2_80mm, max(MxI_2_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax1.scatter(theta_exp_80mm, np.divide(MxI_3_80mm, max(MxI_3_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax1.scatter(theta_exp_80mm_mean, np.divide(Mx_mean_80mm, max(Mx_mean_80mm)), marker="o", s=15., label = '80mm Exp Mean')
ax1.plot(RPB_MCEER_exp[2][int(len(RPB_MCEER_exp[2])/2), 0:int(len(RPB_MCEER_exp[2])/2)], RPB_MCEER_exp[3][int(len(RPB_MCEER_exp[2])/2), 0:int(len(RPB_MCEER_exp[2])/2)]/ max(RPB_MCEER_exp[3][int(len(RPB_MCEER_exp[2])/2), 0:int(len(RPB_MCEER_exp[2])/2)]), 'k')
ax1.plot(theta_80mm_mesh, gtable_80mm[0][:,7]/max(gtable_80mm[0][:,7]), 'k-.', dashes=[12,6,12,6,3,6])
ax1.plot(np.linspace(0,80, num=80), Henrych(255,0.1,0.08, np.linspace(0,80, num=80)) / max(Henrych(255,0.1,0.08, np.linspace(0,80, num=80))), 'k:')
ax1.plot(theta[:,1], Icr_Ir[:,1])
plt.tight_layout()  
#------------------------------------------------------------------------------



















def jang(theta, Ir, beta):
    """
    Beta is a function of Z, Ir is maximum reflected impulse at theta = 0. Theta inputted in degrees.
    """
    return Ir * np.exp( (-beta) * np.tan(np.deg2rad(theta)) * np.tan(np.deg2rad(theta)))

def dharmasena(theta, Ir):
    """
    Beta is a function of Z, Ir is maximum reflected impulse at theta = 0. Theta inputted in degrees.
    """
    return Ir * np.exp( (-1) * np.tan(np.deg2rad(theta)) * np.tan(np.deg2rad(theta)))

def RPB(theta, Ir):
    Ii = 0.1 * Ir
    return Ir * np.cos(np.deg2rad(theta)) * np.cos(np.deg2rad(theta)) + (Ii/Ir)*(1 + (np.cos(np.deg2rad(theta)) * np.cos(np.deg2rad(theta)) ) - 2*np.cos(np.deg2rad(theta)))
    


fig, [ax0,ax1] = plt.subplots(1,2)
ax0.plot(theta, peak_impulse)
ax0.set_xlabel('theta (degrees)')
ax0.set_ylabel('peak specific impulse')
#Plotting Apollo data
ax1.plot(theta, Icr_Ir)
#Plotting Experimental data
ax1.scatter(theta_exp_80mm, np.divide(MxI_1_80mm, max(MxI_1_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = '80mm Exp')
ax1.scatter(theta_exp_80mm, np.divide(MxI_2_80mm, max(MxI_2_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax1.scatter(theta_exp_80mm, np.divide(MxI_3_80mm, max(MxI_3_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax1.scatter(theta_exp_80mm_mean, np.divide(Mx_mean_80mm, max(Mx_mean_80mm)), marker="o", s=15., label = '80mm Exp Mean')
ax1.set_xlabel('theta (degrees)')
ax1.set_ylabel('Ir / Ir Max')
#Also Plotting model distributions
choice = 0
beta = 0.7 
ax1.plot(theta[:, choice],  jang(theta[:, choice], Icr_Ir[:, choice], beta) , label = 'Jang')
ax1.plot(theta[:, choice], RPB(theta[:, choice], Icr_Ir[:, choice]), label = 'RPB')
ax1.plot(theta[:, choice], dharmasena(theta[:, choice], Icr_Ir[:, choice]), label = 'Dharmasena')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='center', bbox_to_anchor=(0.65, 0.90), prop={'size':6})
plt.tight_layout()




choice = 0


def RPB2(theta, Ir):
    Ii = 0 * Ir
    return Ir * np.cos(np.deg2rad(theta)) * np.cos(np.deg2rad(theta)) + (Ii/Ir)*(1 + (np.cos(np.deg2rad(theta)) * np.cos(np.deg2rad(theta)) ) - 2*np.cos(np.deg2rad(theta)))
def RPB2nd_part(theta): 
    return (1 + (np.cos(np.deg2rad(theta)) * np.cos(np.deg2rad(theta)) ) - 2*np.cos(np.deg2rad(theta)))

    
loss = np.subtract(Icr_Ir[:, choice], RPB2(theta[:,choice], max(Icr_Ir[:, choice])))
Cr = np.divide(RPB2nd_part(theta[:,choice]), loss)


fig, [ax0, ax1] = plt.subplots(1,2)
ax1.set_xlabel('theta (degrees)')
ax1.set_ylabel('peak specific impulse')
#Plotting Experimental data
ax1.scatter(theta_exp_80mm, np.divide(MxI_1_80mm, max(MxI_1_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = '80mm Exp')
ax1.scatter(theta_exp_80mm, np.divide(MxI_2_80mm, max(MxI_2_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax1.scatter(theta_exp_80mm, np.divide(MxI_3_80mm, max(MxI_3_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax1.scatter(theta_exp_80mm_mean, np.divide(Mx_mean_80mm, max(Mx_mean_80mm)), marker="o", s=15., label = '80mm Exp Mean')
ax1.set_xlabel('theta (degrees)')
ax1.set_ylabel('Ir / Ir Max')
#Also Plotting model distributions

beta = 0.7 
ax0.plot(theta[:,choice], loss)
ax1.plot(theta[:, choice], Icr_Ir[:, choice], label = 'Apollo')
ax1.plot(theta[:, choice],  jang(theta[:, choice], Icr_Ir[:, choice], beta) , label = 'Jang')
ax1.plot(theta[:, choice], RPB(theta[:, choice], Icr_Ir[:, choice]), label = 'RPB')
ax1.plot(theta[:, choice], RPB2(theta[:, choice], Icr_Ir[:, choice]), label = 'RPB2')
ax1.plot(theta[:, choice], dharmasena(theta[:, choice], Icr_Ir[:, choice]), label = 'Dharmasena')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='center', bbox_to_anchor=(0.65, 0.90), prop={'size':6})
plt.tight_layout()



#--------------------------------------------
#Creating and fitting a gaussian model
def gauss_curve(x, *params):
    y1 = np.zeros_like(x)

    y1 = params[1] * np.exp( -((x - params[0])**2 / (2*(params[2]**2)) ) )

    y = y1
    return y

#Objective function to be minimised. Loss function is MSE.
def lm_residual(params, x, ytrue):
    cen1 = params['cen1'].value
    amp1 = params['amp1'].value
    wid1 = params['wid1'].value

    args = [cen1, amp1, wid1]
    
    global current_cost
    current_cost = mean_squared_error(ytrue, gauss_curve(x, *args))  
    return current_cost

#Create a set of parameters
params = lm.Parameters()
params['cen1'] = lm.Parameter(name='cen1', value = 0.5, vary = 'false', min = 0.4, max = 0.6)
params['amp1'] = lm.Parameter(name='amp1', value = 0.5, min=0, max=1)
params['wid1'] = lm.Parameter(name='wid1', value = 0.5, min=0, max=1)   

#
test = 0
x = np.concatenate([-np.flipud(theta[:,test]),theta[:,test]])
x = (x - min(x)) / (max(x) - min(x))
data = np.concatenate([np.flipud(Icr_Ir[:,test]), Icr_Ir[:,test]])
#data = np.concatenate([np.flipud(Icr_Ir.max(1)), Icr_Ir.max(1)])
data = (data - min(data)) / (max(data) - min(data)) 

#Nonlinear regression
result = lm.minimize(lm_residual, params, method = 'least_squares', args = (x, data))
gaussmod = gauss_curve(x, result.params['cen1'].value, result.params['amp1'].value, result.params['wid1'].value)





#Plotting some more model functions
ax0.scatter(theta, Icr_Ir, marker='o', s=0.1, color = 'red' , label = 'CFD')
#Plotting Experimental data
ax0.scatter(theta_exp_80mm, np.divide(MxI_1_80mm, max(MxI_1_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = '80mm Exp')
ax0.scatter(theta_exp_80mm, np.divide(MxI_2_80mm, max(MxI_2_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax0.scatter(theta_exp_80mm, np.divide(MxI_3_80mm, max(MxI_3_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
ax0.scatter(theta_exp_80mm_mean, np.divide(Mx_mean_80mm, max(Mx_mean_80mm)), marker="o", s=15., label = '80mm Exp Mean')
ax0.set_xlabel('theta (degrees)')
ax0.set_ylabel('Ir / Ir Max')
#Gaussian
ax0.plot(theta[:,test], gaussmod[300::], color = 'black', label = 'Gaussian')
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles, labels, loc='center', bbox_to_anchor=(0.65, 0.90), prop={'size':6})
plt.tight_layout()