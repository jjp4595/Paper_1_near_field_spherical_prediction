"""
Peak impulse predictor for theta. 
"""

import preamble_functions as pre
import matplotlib.pyplot as plt #3.0.2
from matplotlib import cm
import numpy as np #1.15.4
import blast as blst
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import scipy.io as sio
import lmfit as lm
from sklearn.metrics import r2_score, mean_squared_error
import lmfit as lm

plt.rcParams["font.family"] = "cmr10" #Set Graph fonts to cmr10

#Import Apollo data
#Apollo_FileList = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\2cr\*.txt")
#Apollo_gtable = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\2cr\*gtable",1)
#Apollo_gauges = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\2cr\*gauges",1)
Apollo_FileList = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*.txt")
Apollo_gtable = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*gtable",1)
Apollo_gauges = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*gauges",1)


#Quick test to see if enough time in simulation
#fig00, [ax_test1, ax_test2] = plt.subplots(1,2)
#ax_test1.plot(Apollo_gauges[4][:,0],Apollo_gauges[4][:,200])
#ax_test2.plot(Apollo_gauges[4][:,0],Apollo_gauges[4][:,400])

#Some charge properties
charge_rad = 0.0246
charge_mass = 0.1


#Finding theta and creating data structures


clear_standoff=np.zeros((len(Apollo_FileList),1))
peak_impulse=[]
#incident_impulse = []
#ref_factor = [] 
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
#incident_impulse = np.stack(incident_impulse, axis = 1)
#ref_factor = np.stack(ref_factor, axis = 1)
Icr_Ir = np.stack(Icr_Ir, axis =1)
clear_standoff = np.transpose(np.repeat(clear_standoff, len(theta), axis=1))





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



#Model Functions
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
    



#Data Visualisations
#surface
fig = plt.figure()
ax = Axes3D(fig)
#ax.scatter(standoff, theta, all_impulses);
ax.plot_surface(clear_standoff, theta, peak_impulse, cmap=cm.coolwarm);
ax.set_title('peak specific impulse surface for theta')
ax.set_xlabel('standoff (m)')
ax.set_ylabel('incident wave angle')
ax.set_zlabel('peak specific impulse')
plt.show()

#
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










#---------------------Current playing around!!!!!----------------------------------------------

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