"""
Peak impulse predictor for theta. 
"""

import preamble_functions as pre
import matplotlib.pyplot as plt #3.0.2
from matplotlib import cm
import numpy as np #1.15.4
import blast as blst
from mpl_toolkits.mplot3d import Axes3D
import scipy.io as sio
import lmfit as lm
from sklearn.metrics import r2_score, mean_squared_error

plt.rcParams["font.family"] = "cmr10" #Set Graph fonts to cmr10

#Import Apollo data
Apollo_FileList = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\PE4_100g_theta80_z055_16\*.txt")
Apollo_gtable = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\PE4_100g_theta80_z055_16\*gtable",1)

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
    Icr_Ir.append(peak_impulse[i] / max(peak_impulse[i]))
    theta.append(np.rad2deg(np.arctan(np.divide(Apollo_gtable[i][:,2], clear_standoff[i] + charge_rad))))
theta = np.stack(theta, axis=1)
peak_impulse = np.stack(peak_impulse, axis = 1)
Icr_Ir = np.stack(Icr_Ir, axis =1)
clear_standoff = np.transpose(np.repeat(clear_standoff, len(theta), axis=1))


#Load Data for 80mm and 380mm Apollo Experimental
NF_80mm_exp = sio.loadmat(r"C:\Users\jorda\Google Drive\Apollo Sims\Near Field Sims\100gPE4Sphere_80mm") 
NF_380mm_exp = sio.loadmat(r"C:\Users\jorda\Google Drive\Apollo Sims\Near Field Sims\100gPE4Sphere_380mm")
coords = np.append(np.arange(-100,101,25), np.arange(-100,-24,25)) 
coords = np.append(coords, np.arange(25,101,25))
coords = np.divide(coords,1000)
coords_mean = np.arange(0,0.101,0.025)
theta_exp_80mm = np.rad2deg(np.arctan(np.abs(coords)/0.08))
theta_exp_80mm_mean = np.rad2deg(np.arctan(np.abs(coords_mean)/0.08))
theta_exp_380mm = np.rad2deg(np.arctan(np.abs(coords)/0.38))
theta_exp_380mm_mean = np.rad2deg(np.arctan(np.abs(coords_mean)/0.38))
MxI_1_80mm = np.transpose(NF_80mm_exp['MxI'][:,:,0])
MxI_2_80mm = np.transpose(NF_80mm_exp['MxI'][:,:,1])
MxI_3_80mm = np.transpose(NF_80mm_exp['MxI'][:,:,2])
Mx_mean_80mm = np.transpose(NF_80mm_exp['MEANI'])
MxI_1_380mm = np.transpose(NF_380mm_exp['MxI'][:,:,0])
MxI_2_380mm = np.transpose(NF_380mm_exp['MxI'][:,:,1])
MxI_3_380mm = np.transpose(NF_380mm_exp['MxI'][:,:,2])
Mx_mean_380mm = np.transpose(NF_380mm_exp['MEANI'])


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
    
 




#Creating and fitting a gaussian model ----------------------------------------
def gauss_two_curve(x, *params):
    y1 = np.zeros_like(x)
    y2 = np.zeros_like(x)
    y1 = params[1] * np.exp( -((x - params[0])**2 / (2*(params[2]**2)) ) )
    y2 = params[4] * np.exp( -((x - params[3])**2 / (2*(params[5]**2)) ) )
    y = y1 + y2
    return y
def pen1(param1, param2):
    delta = param2 - param1
    out = max(0,delta)
    return out**0.5  
#Objective function to be minimised. Loss function is MSE.
def lm_residual(params, x, ytrue):
    cen1 = params['cen1'].value
    amp1 = params['amp1'].value
    wid1 = params['wid1'].value
    cen2 = params['cen2'].value
    amp2 = params['amp2'].value
    wid2 = params['wid2'].value

    args = [cen1, amp1, wid1, cen2, amp2, wid2]
    
    global current_cost
    current_cost = ( (mean_squared_error(ytrue, gauss_two_curve(x, *args)))**2 + (reg1 * pen1(amp1, amp2)**2) ) **0.5
    #current_cost = ( (mean_absolute_error(ytrue, gauss_two_curve(x, *args)))**2 + (reg1 * pen1(amp1, amp2)**2) ) **0.5
    return current_cost
#Create a set of parameters
params = lm.Parameters()
params['cen1'] = lm.Parameter(name='cen1', value = 0.5, vary = 'false', min = 0.49, max = 0.51)
params['amp1'] = lm.Parameter(name='amp1', value = 0.5, min=0, max=1)
params['wid1'] = lm.Parameter(name='wid1', value = 0.5, min=0, max=1)   
params['cen2'] = lm.Parameter(name='cen2', value = 0.5, vary = 'false', min = 0.49, max = 0.51)
params['amp2'] = lm.Parameter(name='amp2', value = 0.3, min=0, max=1)
params['wid2'] = lm.Parameter(name='wid2', value = 0.5, min=0, max=1)

reg1 = 1
gauss_theta = np.concatenate([-np.flipud(theta[:,0]),theta[:,0]])
gauss_theta_norm = (gauss_theta - min(gauss_theta)) / (max(gauss_theta) - min(gauss_theta))
gauss_Icr_Ir = np.concatenate([np.flipud(Icr_Ir[:,0]), Icr_Ir[:,0]])
result = lm.minimize(lm_residual, params, method = 'least_squares', args = (gauss_theta_norm, gauss_Icr_Ir))
opt = result.params

gauss_mod = gauss_two_curve(gauss_theta, opt['cen1'].value, opt['amp1'].value, opt['wid1'].value, opt['cen2'].value, opt['amp2'].value, opt['wid2'].value)












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
#ax1.scatter(theta_exp_380mm, np.divide(MxI_1_380mm, max(MxI_1_380mm)), marker="+", s=15., color=[0.5,0.5,0.5], edgecolors='none', label = '380mm Exp')
#ax1.scatter(theta_exp_380mm, np.divide(MxI_2_380mm, max(MxI_2_380mm)), marker="+", s=15., color=[0.5,0.5,0.5], edgecolors='none')
#ax1.scatter(theta_exp_380mm, np.divide(MxI_3_380mm, max(MxI_3_380mm)), marker="+", s=15., color=[0.5,0.5,0.5], edgecolors='none')
#ax1.scatter(theta_exp_380mm_mean, np.divide(Mx_mean_380mm, max(Mx_mean_380mm)), marker="p", s=15., label = '380mm Exp Mean')
ax1.set_xlabel('theta (degrees)')
ax1.set_ylabel('Ir / Ir Max')
#Also Plotting model distributions
choice = 0
beta = 0.8 
ax1.plot(theta[:, choice],  jang(theta[:, choice], 1, beta) , label = 'Jang')
ax1.plot(theta[:, choice], RPB(theta[:, choice], 1), label = 'RPB')
ax1.plot(theta[:, choice], dharmasena(theta[:, choice], 1), label = 'Dharmasena')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='center', bbox_to_anchor=(0.65, 0.90), prop={'size':6})
plt.tight_layout()


