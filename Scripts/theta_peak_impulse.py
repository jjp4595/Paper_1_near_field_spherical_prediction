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
from sklearn.metrics import mean_squared_error, r2_score
from MCEER_curves import MCEER
from scipy import ndimage
import os
from scipy import stats


#plt.rcParams["font.family"] = "cmr10" #Set Graph fonts to cmr10
params = {'font.family':'serif',
        'axes.labelsize':'small',
        'xtick.labelsize':'x-small',
        'ytick.labelsize':'x-small',
#        
        'lines.markersize': 10,
        'scatter.marker': 'o',
#        
        'legend.fontsize':'small',
        'legend.title_fontsize':'small',
        'legend.fancybox': True,
        'legend.framealpha': 0.5,
        'legend.shadow': False,
        'legend.frameon': True,
#        
#
        'grid.linestyle':'--',
        'grid.linewidth':'0.5',
        'lines.linewidth':'0.5'}

plt.rcParams.update(params)
#plt.style.use('paper_stylesheet.mpltstyle')

#Import Apollo data
# Apollo_FileList = pre.FileAddressList(os.path.join(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16\*.txt"))
# Apollo_gtable = pre.FileAddressList(os.path.join(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\main_z055_16\*gtable"),1)

Apollo_FileList = pre.FileAddressList(os.path.join(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*.txt"))
Apollo_gtable = pre.FileAddressList(os.path.join(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\Sphere\original_PE4_100g_theta80_z055_16\*gtable"),1)
Apollo_FileList = Apollo_FileList[4::]
Apollo_gtable = Apollo_gtable[4::]


#Some charge properties
charge_rad = 0.0246
charge_mass = 0.1

#Finding theta and creating data structures ----------------------------------
clear_standoff=np.zeros((len(Apollo_gtable),1))
peak_impulses=[]
Icr_Ir = []
z_dataset = []

theta = []
for i in range(len(Apollo_gtable)):
    clear_standoff[i] = pre.standoff_func(Apollo_FileList[i]) - charge_rad
    peak_impulses.append(Apollo_gtable[i][:,7])
    #incident_impulse.append(Apollo_gtable[i][incident_index::,7])
    #ref_factor.append(peak_impulse[i]/incident_impulse[i])
    Icr_Ir.append(peak_impulses[i] / max(peak_impulses[i]))
    theta.append(np.rad2deg(np.arctan(np.divide(Apollo_gtable[i][:,2], clear_standoff[i] + charge_rad))))

theta = np.stack(theta, axis=1)
peak_impulse = np.stack(peak_impulses, axis = 1)
Icr_Ir = np.stack(Icr_Ir, axis =1)
clear_standoff = np.transpose(np.repeat(clear_standoff, len(theta), axis=1))
z_dataset = np.divide(clear_standoff, charge_mass**(1/3)) 
clear_standoff = np.divide(clear_standoff, charge_rad) #multiples of clear charge radii
#-----------------------------------------------------------------------------


def import_exp():
    gtable_80mm = pre.FileAddressList(os.path.join(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Near Field Sims\Sims\Latest\80mm_with_afterburn\*gtable"),1)
    theta_80mm_mesh = np.rad2deg(np.arctan(gtable_80mm[0][:,2]/0.08))
    return gtable_80mm, theta_80mm_mesh
gtable_80mm, theta_80mm_mesh = import_exp()

def overlay_exps(ax, is_crit = None):
    #Load Data for 80mm and 380mm Apollo Experimental ----------------------------
    NF_80mm_exp = sio.loadmat(os.path.join(os.environ['USERPROFILE'] + r"\Google Drive\Apollo Sims\Near Field Sims\100gPE4Sphere_80mm") )

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
    if is_crit is None:
        ax.scatter(theta_exp_80mm, MxI_1_80mm/1e3, marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = 'Exp')
        ax.scatter(theta_exp_80mm, MxI_2_80mm/1e3, marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
        ax.scatter(theta_exp_80mm, MxI_3_80mm/1e3, marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
        ax.scatter(theta_exp_80mm_mean, Mx_mean_80mm/1e3, marker="o", s=15., label = 'Exp - mean')
    else:
        ax.scatter(theta_exp_80mm, np.divide(MxI_1_80mm, max(MxI_1_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none', label = 'Exp')
        ax.scatter(theta_exp_80mm, np.divide(MxI_2_80mm, max(MxI_2_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
        ax.scatter(theta_exp_80mm, np.divide(MxI_3_80mm, max(MxI_3_80mm)), marker="x", s=15., color=[0.75,0.75,0.75], edgecolors='none')
        ax.scatter(theta_exp_80mm_mean, np.divide(Mx_mean_80mm, max(Mx_mean_80mm)), marker="o", s=15., label = 'Exp - mean')
#------------------------------------------------------------------------------    




#Graph1 - dataset- overview ---------------------------------------------------
def graph_dataset_overview():   
    fig0, ax = plt.subplots(1,1)
    fig0.set_size_inches(3, 2.5)
    CS = ax.contourf(theta, clear_standoff, peak_impulse/1e3, levels = [0 , 2.5, 5, 7.5, 10, 12.5, 15], cmap = plt.cm.cubehelix)
    cbar = fig0.colorbar(CS)
    cbar.ax.invert_yaxis()
    cbar.ax.set_ylabel('peak specific impulse (MPa.ms)')
    ax.set_ylabel('standoff (clear charge radii)')
    ax.set_xlabel('incident wave angle (degrees)')
    plt.tight_layout()
    fig0.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\theta_peak_impulse_0.pdf"), format = 'pdf')
    return fig0
graph_dataset_overview()


#Graph 2 Power law-------------------------------------------------------------
def graph_powerlaw():
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(clear_standoff[0,:]), np.log10(peak_impulse[0,:]/1e3))
    linear_model = intercept + slope*np.log10(clear_standoff[0,:])
    const = np.divide(peak_impulse[0,:]/1e3, clear_standoff[0,:]**slope)
    const = const.sum()/len(const)
    residuals_power = np.log10(peak_impulse[0,:]/1e3) - linear_model
    RSS_power = np.power(residuals_power,2).sum()
    RSE_power = ( (1/(len(residuals_power)-2)) * RSS_power   )
    text_ax1 = "$R^2 =$" + str(round(r2_score(np.log10(peak_impulse[0,:]/1e3), linear_model), 3))
    
    if p_value < 0.0001:
        text_ax3_p = "$p < 0.0001$"
    else:
        text_ax3_p = "$p = {:.3f}$".format(p_value)
    
    if RSE_power < 0.0001:
        text_ax3_se = "$RSE < 0.0001$"
    else:
        text_ax3_se = "$RSE = {:.3f}$".format(RSE_power)
    
    
    text_ax2 = "$f(x) = {%.3f}.x^{%.3f}$" % (intercept, slope)
    
    fig1, [ax1, ax3, ax2] = plt.subplots(1,3)
    fig1.set_size_inches(7.5, 2.5)
    
    ax1.scatter(np.log10(clear_standoff[0,:]), np.log10(peak_impulse[0,:]/1e3), label = 'CFD data')
    ax1.plot(np.log10(clear_standoff[0,:]), linear_model, label='fitted model')
    ax1.plot(np.log10(clear_standoff[0,:]), linear_model + (2 * RSE_power), '--k', alpha = 0.2, label='95% PI')
    ax1.plot(np.log10(clear_standoff[0,:]), linear_model - (2 * RSE_power), '--k', alpha = 0.2)
    ax1.text(0.1, 0.4, text_ax1, transform=ax1.transAxes)
    ax1.text(0.1, 0.3, text_ax3_p, transform=ax1.transAxes)
    ax1.text(0.1, 0.2, text_ax3_se, transform=ax1.transAxes)
    ax1.set_ylabel('log(peak specific impulse (MPa.ms))')
    ax1.set_xlabel('log(standoff (clear charge radii))')
    
    #marker="o", s=10, edgecolors = (0,0,0,0.4), facecolor=(0,0,0,0.4),
    ax2.scatter(clear_standoff[0,:], peak_impulse[0,:]/1e3,  label = 'CFD data')
    ax2.plot(clear_standoff[0,:], const * clear_standoff[0,:]**slope, 'r', label='fitted model')
    ax2.text(0.35, 0.7, text_ax2, fontsize = 'small', transform=ax2.transAxes)
    ax2.set_ylabel('peak specific impulse (MPa.ms)')
    ax2.set_xlabel('standoff (clear charge radii)')
    
    ax3.scatter(np.log10(clear_standoff[0,:]), residuals_power, label = 'Residuals')
    ax3.set_ylim(-0.1,0.1)
    ax3.set_ylabel('Residual')
    ax3.set_xlabel('log(standoff (clear charge radii))')
    
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)
    ax1.locator_params(axis = 'both',tight=True, nbins=6)
    ax2.locator_params(axis = 'both',tight=True, nbins=6)
    ax3.locator_params(axis = 'both',tight=True, nbins=6)
    plt.tight_layout()
    fig1.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\power_fit.pdf"), format = 'pdf')
    return fig1, slope, intercept
fig_power, slope, intercept = graph_powerlaw()
#-----------------------------------------------------------------------------


#Models in the literature: ---------------------------------------------------
#1) Henrych Model
Qw = 4.79e6/4184 #Specific energy of PE4, 4.79 assumed from EOS.dat
Ux = (2 * Qw)**0.5
A0 = (Ux * 1.8) / (4*np.pi)
def Henrych(A0, W, a, theta):
    """
    Calculating specific impulse
    """
    i = ((A0 * W) / (a**2)) * np.cos(np.deg2rad(theta))**4
    return i
def Henrych_I(A0, W, theta_lim):
    """
    Calculating total impulse across a circle
    """
    return np.pi * A0 * W * np.sin(np.deg2rad(theta_lim))**2



#2) RPB - MCEER Model hybrid
"""
Generating impulse from MCEER curves and modelling distribution via RPB model.
"""
def RPB_MCEER_i(W, R, theta_lim):
    """
    imp and I output are in MPa.ms
    """
    #Wroot = W**(1/3)
    A = R * np.tan(np.deg2rad(theta_lim)) * 2
    B = A
    res = 200
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
    #Rectangular area
    Area = np.ones_like(imp)
    Area = Area * (x[1]-x[0]) * (y[1]-y[0])
    #Boundary conds for rectangular target
    ##Area[:,0]  *= 0.5
    ##Area[0,:]  *= 0.5
    ##Area[:,-1] *= 0.5
    ##Area[-1,:] *= 0.5
    
    #Circular mask over area for circular target 
    r = int(res/2)
    cen_x, cen_y = r, r
    y1, x1 = np.ogrid[-cen_x:res-cen_x, -cen_y:res-cen_y]
    mask = x1 * x1 + y1 * y1 <= r*r
    
    array = np.zeros((res, res))
    array[mask] = 1
    struct = ndimage.generate_binary_structure(2, 2)
    erode = ndimage.binary_erosion(mask, struct)
    edges = mask ^ erode
    array[edges] = 0.5
    
    Area = np.multiply(array, Area)
    
    Imp = np.multiply(imp, Area)
    I = np.sum(Imp)
    return [X,Y, theta_MCEER, imp, I]


#3) CFD Impulse generation for circular target --------------------------------
def Impulse_CFD(specific_impulse, R, theta_lim, theta):
    """
    i is impulse distribution to be interpolated
    R is clear standoff
    theta_lim is the upper limit of theta for size of the circle (from bomb POV)
    theta is the theta coords for interpolation
    output is in Pa.s
    """
	#Plate dimensions
    A = R * np.tan(np.deg2rad(theta_lim)) * 2
    B = A
    res = 200
    x = np.linspace(-A/2, A/2, res)
    y = np.linspace(-B/2, B/2, res)
    [X,Y] = np.meshgrid(x,y)
    theta_array = np.rad2deg(np.arctan(np.divide(np.power(np.add(np.power(X,2), np.power(Y,2)), 0.5), R))) 
	
    ir = np.zeros((res,res))

    for i in range(len(ir)):
        for j in range(len(ir)):
            ir[i,j] = np.interp(theta_array[i,j], theta, specific_impulse) 
				
	#Rectangular area
    Area = np.ones_like(ir)
    Area = Area * (x[1]-x[0]) * (y[1]-y[0])

    
    #Circular mask over area for circular target 
    r = int(res/2)
    cen_x, cen_y = r, r
    y1, x1 = np.ogrid[-cen_x:res-cen_x, -cen_y:res-cen_y]
    mask = x1 * x1 + y1 * y1 <= r*r
    
    array = np.zeros((res, res))
    array[mask] = 1
    struct = ndimage.generate_binary_structure(2, 2)
    erode = ndimage.binary_erosion(mask, struct)
    edges = mask ^ erode
    array[edges] = 0.5
    
    Area = np.multiply(array, Area)
    
    Imp = np.multiply(ir, Area)
    I = np.sum(Imp)
    return I

#4) Jang Model for Icr --------------------------------------------------------
def jang(theta, Ir, beta):
    """
    Beta is a function of Z, Ir is maximum reflected impulse at theta = 0. Theta inputted in degrees.
    """
    return Ir * np.exp( (-beta) * np.tan(np.deg2rad(theta)) * np.tan(np.deg2rad(theta)))

#5) Dharmasena model for Icr --------------------------------------------------
def dharmasena(theta, Ir):
    """
    Beta is 1
    """
    return Ir * np.exp( (-1) * np.tan(np.deg2rad(theta)) * np.tan(np.deg2rad(theta)))

#6) RPB Model -----------------------------------------------------------------
def RPB(theta, Ir, frac):
    """
    Parameters
    ----------
    theta : TYPE
        DESCRIPTION.
    Ir : TYPE
        DESCRIPTION.
    frac : Ii = frac * Ir.


    """
    Ii = frac * Ir
    return Ir * np.cos(np.deg2rad(theta)) * np.cos(np.deg2rad(theta)) + (Ii/Ir)*(1 + (np.cos(np.deg2rad(theta)) * np.cos(np.deg2rad(theta)) ) - 2*np.cos(np.deg2rad(theta)))
 

#Compare 
def graph_impulse_comparisons():      
    
    RPB_MCEER_exp = RPB_MCEER_i(0.1, 0.08, 80)
    i_jang = jang(np.linspace(0,80, num=80), 1, 0.8)
    i_dharmasena = dharmasena(np.linspace(0,80,num=80), 1)
    
    #Graph of RPB & MCEER contour
    fig2, ax = plt.subplots(1,1)
    fig2.set_size_inches(3, 2.5)
    CS = ax.contourf(RPB_MCEER_exp[0], RPB_MCEER_exp[1], RPB_MCEER_exp[3], cmap = plt.cm.cubehelix)
    cbar = fig2.colorbar(CS)
    cbar.ax.set_ylabel('peak specific impulse (MPa.ms)')
    ax.set_ylabel('x-position')
    ax.set_xlabel('y-position')
    plt.tight_layout()
    fig2.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\theta_peak_impulse_RPB_MCEER_plate.pdf"), format = 'pdf')



    #Graph comparing different models    
    fig3, [ax0,ax1] = plt.subplots(1,2)
    fig3.set_size_inches(5,2.5)
    overlay_exps(ax0)
    ax0.plot(RPB_MCEER_exp[2][int(len(RPB_MCEER_exp[2])/2), 0:int(len(RPB_MCEER_exp[2])/2)], RPB_MCEER_exp[3][int(len(RPB_MCEER_exp[2])/2),0:int(len(RPB_MCEER_exp[2])/2)], 'k', label = 'MCEER-RPB')
    ax0.plot(theta_80mm_mesh, gtable_80mm[0][:,7]/1e3, 'k-.', dashes=[12,6,12,6,3,6], label = 'CFD - 1.25mm')
    ax0.plot(np.linspace(0,80, num=80), Henrych(A0,0.1,0.08, np.linspace(0,80, num=80))/1e3, 'k:', label = 'Henrych')
    ax0.set_xlabel('theta (degrees)')
    ax0.set_ylabel('peak specific impulse (MPa.ms)')
    handles, labels = ax0.get_legend_handles_labels()
    ax0.legend(handles, labels)
    
    ax1.set_xlabel('theta (degrees)')
    ax1.set_ylabel('peak specific impulse ratio')
    overlay_exps(ax1, 1)
    ax1.plot(RPB_MCEER_exp[2][int(len(RPB_MCEER_exp[2])/2), 0:int(len(RPB_MCEER_exp[2])/2)], RPB_MCEER_exp[3][int(len(RPB_MCEER_exp[2])/2), 0:int(len(RPB_MCEER_exp[2])/2)]/ max(RPB_MCEER_exp[3][int(len(RPB_MCEER_exp[2])/2), 0:int(len(RPB_MCEER_exp[2])/2)]), 'k')
    ax1.plot(theta_80mm_mesh, gtable_80mm[0][:,7]/max(gtable_80mm[0][:,7]), 'k-.', dashes=[12,6,12,6,3,6])
    ax1.plot(np.linspace(0,80, num=80), Henrych(255,0.1,0.08, np.linspace(0,80, num=80)) / max(Henrych(255,0.1,0.08, np.linspace(0,80, num=80))), 'k:')
    ax1.plot(np.linspace(0,80, num=80), i_jang, 'k', marker='o', markevery=16,  label = 'Jang')
    ax1.plot(np.linspace(0,80, num=80), i_dharmasena, 'k-.', marker='D', markevery=16,  label = 'Dharmasena')
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)
    plt.tight_layout()  
    fig3.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\theta_peak_impulse_i_theta_comparisons.pdf"), format = 'pdf')
graph_impulse_comparisons()
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
#Generating %Difference for where z_dataset is greater than lower limit of MCEER
# z_I_percent = z_dataset[0,4::]
# theta_I_percent = np.linspace(0,80,num=40)

# Imp_Henrych = np.zeros((len(theta_I_percent), len(z_I_percent)))
# Imp_RPB_MCEER = np.zeros_like(Imp_Henrych)
# Imp_CFD = np.zeros_like(Imp_Henrych)

# for j in range(np.size(Imp_Henrych,1)):
#     for i in range(np.size(Imp_Henrych,0)):       
#         Imp_RPB_MCEER[i,j] = RPB_MCEER_i(charge_mass, z_I_percent[j]*(charge_mass**(1/3)), theta_I_percent[i])[4]
#         Imp_Henrych[i,j] = Henrych_I(A0, charge_mass, theta_I_percent[i])
#         Imp_CFD[i,j] = Impulse_CFD(peak_impulses[j+4], z_I_percent[j]*(charge_mass**(1/3)), theta_I_percent[i], theta[:,j+4])
#     print(j/np.size(Imp_Henrych,1))


# fig4a, ax0 = plt.subplots(1,1)
# fig4a.set_size_inches(3, 2.5)
# #CS0 = ax0.contourf(theta_I_percent, z_I_percent, np.transpose(Imp_RPB_MCEER), cmap = plt.cm.cubehelix)
# #cbar = fig4a.colorbar(CS0, ax=ax0)
# cbar.ax.set_ylabel('peak specific impulse (MPa.ms)')
# ax0.set_ylabel('Clear scaled distance, z$(m/kg^{1/3})$')
# ax0.set_xlabel('$Theta_{lim}$')
# plt.tight_layout()
# #fig4a.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\theta_peak_impulse_i_percent_a.pdf"), format = 'pdf')

# fig4b, ax0 = plt.subplots(1,1)
# fig4b.set_size_inches(3, 2.5)
# #CS0 = ax0.contourf(theta_I_percent, z_I_percent, np.transpose(Imp_Henrych), cmap = plt.cm.cubehelix)
# #cbar = fig4b.colorbar(CS0, ax=ax0)
# cbar.ax.set_ylabel('peak specific impulse (MPa.ms)')
# ax0.set_ylabel('Clear scaled distance, z$(m/kg^{1/3})$')
# ax0.set_xlabel('$Theta_{lim}$')
# plt.tight_layout()
# #fig4b.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\theta_peak_impulse_i_percent_b.pdf"), format = 'pdf')

# fig4c, ax0 = plt.subplots(1,1)
# fig4c.set_size_inches(3, 2.5)
# #CS0 = ax0.contourf(theta_I_percent, z_I_percent, np.transpose(Imp_CFD)/1e3, cmap = plt.cm.cubehelix)
# #cbar = fig4c.colorbar(CS0, ax=ax0)
# cbar.ax.set_ylabel('peak specific impulse (MPa.ms)')
# ax0.set_ylabel('Clear scaled distance, z$(m/kg^{1/3})$')
# ax0.set_xlabel('$Theta_{lim}$')
# plt.tight_layout()
# #fig4c.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\theta_peak_impulse_i_percent_c.pdf"), format = 'pdf')
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#Model fitting Stuff that needs refinement 
#------------------------------------------------------------------------------

#1) One Gaussian
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
params = lm.Parameters()
params['cen1'] = lm.Parameter(name='cen1', value = 0.5, vary = 'false', min = 0.4, max = 0.6)
params['amp1'] = lm.Parameter(name='amp1', value = 0.5, min=0, max=1)
params['wid1'] = lm.Parameter(name='wid1', value = 0.5, min=0, max=1)   



#2) Two Gaussian 
params2 = lm.Parameters()
params2['cen1'] = lm.Parameter(name='cen1', value = 0.5, vary = 'false', min = 0.4, max = 0.6)
params2['amp1'] = lm.Parameter(name='amp1', value = 0.5, min=0, max=1)
params2['wid1'] = lm.Parameter(name='wid1', value = 0.5, min=0, max=1)  
params2['cen2'] = lm.Parameter(name='cen2', value = 0.5, vary = 'false', min = 0.4, max = 0.6)
params2['amp2'] = lm.Parameter(name='amp2', value = 0.5, min=0, max=1)
params2['wid2'] = lm.Parameter(name='wid2', value = 0.5, min=0, max=1)  
def gauss_two_curve(x, *params):
    y1 = np.zeros_like(x)
    y2 = np.zeros_like(x)
    y1 = params[1] * np.exp( -((x - params[0])**2 / (2*(params[2]**2)) ) )
    y2 = params[4] * np.exp( -((x - params[3])**2 / (2*(params[5]**2)) ) )
    y = y1 + y2
    return y

#Amplitude penalty function
def pen1(param1, param2):
    delta = param2 - param1
    out = max(0,delta)
    return out**0.5  

#Objective function to be minimised. Loss function is MSE.
def lm_residual2(params, x, ytrue):
    cen1 = params['cen1'].value
    amp1 = params['amp1'].value
    wid1 = params['wid1'].value
    cen2 = params['cen2'].value
    amp2 = params['amp2'].value
    wid2 = params['wid2'].value


    args = [cen1, amp1, wid1, cen2, amp2, wid2]
    reg1 = 5
    current_cost = ( (mean_squared_error(ytrue, gauss_two_curve(x, *args)))**2 + (reg1 * pen1(amp1, amp2)**2) ) **0.5
    return current_cost

#------------------------------------------------------------------------------




#gauss model current fitted to the closest test
test = 0


#data_o = np.concatenate([np.flipud(Icr_Ir[:,test]), Icr_Ir[:,test]]) #Fit data to closest sample
#data_o = np.concatenate([np.flipud(Icr_Ir.max(1)), Icr_Ir.max(1)]) #Fit data to max Icr_Ir
data_o = np.concatenate([np.flipud(Icr_Ir.mean(1)), Icr_Ir.mean(1)]) #Fit data to mean
data = (data_o - min(data_o)) / (max(data_o) - min(data_o))


x = np.linspace(0,80,len(data))
x = (x - min(x)) / (max(x) - min(x))


#Nonlinear regression -gaussians
result = lm.minimize(lm_residual, params, method = 'least_squares', args = (x, data))
gaussmod = gauss_curve(x, result.params['cen1'].value, result.params['amp1'].value, result.params['wid1'].value)

result2 = lm.minimize(lm_residual2, params2, method = 'least_squares', args = (x, data))
gaussmod2 = gauss_two_curve(x, result2.params['cen1'].value, result2.params['amp1'].value, result2.params['wid1'].value, result2.params['cen2'].value, result2.params['amp2'].value, result2.params['wid2'].value)

#Linearregressions
poly6 = np.polynomial.polynomial.Polynomial.fit(x, data, 6)


# #FIGURE OUT WHAT IS GOING ON HERE
# y = peak_impulse.mean(1)
# recipr_1 = np.polynomial.polynomial.Polynomial.fit(x, 1/y, 1)
# recipr_1 = recipr_1.__call__(1/x)

# fig, [ax, ax2, ax3] = plt.subplots(1,3)
# ax.plot(1/x, recipr_1)
# ax2.plot(x, 1/recipr_1)
# ax3.plot(theta.mean(1), 1/recipr_1)
# #recipr_1 = recipr_1[int(len(gaussmod)/2)::]
# recipr_1 = 1/recipr_1
# recipr_1 = (recipr_1 - min(recipr_1)) / (max(recipr_1) - min(recipr_1))




#impulse distribution in theta -----------------------------------------------
def my_model_graphs():
    fig, [ax0, ax1] = plt.subplots(1,2)
    fig.set_size_inches(10,6)
    
    #Plotting some more model functions
    
    ax0.plot(theta.mean(1), Icr_Ir.mean(1), markevery=10, label = 'CFD - mean')
    ax0.fill_between(theta.mean(1), Icr_Ir.min(1), Icr_Ir.max(1), alpha = 0.2, label = 'CFD - range')
    
    #Plotting Experimental data
    overlay_exps(ax0, 1)
    
    #Gaussian models
    ax0.plot(theta.mean(1), gaussmod[int(len(gaussmod)/2)::], label = 'One gaussian')
    ax0.plot(theta.mean(1), gaussmod2[int(len(gaussmod)/2)::], '--',label = 'Two gaussian')
    text_gaussian = r"$f(\theta) = exp \left( \frac{-(\theta-0.5)^2}{2 \times {%.3f}^2} \right)$" % (result.params['wid1'].value)
    
    #poly
    ax0.plot(theta.mean(1), poly6.__call__(x)[int(len(gaussmod)/2)::])
    ax0.text(0.6, 0.7, text_gaussian, fontsize = 'small', transform=ax0.transAxes)
    
    ##recip
    #ax0.plot(theta.mean(1), recipr_1.__call__(1/x)[int(len(gaussmod)/2)::])
    #https://statisticsbyjim.com/regression/curve-fitting-linear-nonlinear-regression/
    #residual plots
    ax1.scatter(theta.mean(1), Icr_Ir.mean(1) - poly6.__call__(x)[int(len(gaussmod)/2)::], label = 'Poly')
    ax1.scatter(theta.mean(1), Icr_Ir.mean(1) - gaussmod[int(len(gaussmod)/2)::],  label = 'One gauss')
    ax1.scatter(theta.mean(1), Icr_Ir.mean(1) - gaussmod2[int(len(gaussmod)/2)::], label = 'Two gauss')
    ax1.set_ylim(-0.1, 0.1)
    
    #axis settings
    ax0.set_xlabel('theta (degrees)')
    ax0.set_ylabel('Ir / Ir Max')
    handles, labels = ax0.get_legend_handles_labels()
    ax0.legend(handles, labels)
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles, labels)
    ax1.set_xlabel('theta (degrees)')
    ax1.set_ylabel('Residual')
    ax0.locator_params(axis = 'both',tight=True, nbins=6)
    ax1.locator_params(axis = 'both',tight=True, nbins=6)
    plt.tight_layout()
    fig.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Temp\theta_predictor.pdf"), format = 'pdf')
my_model_graphs()




#--TESTING---
#Build surface from Guassian Eq and power law.
crs = np.linspace(clear_standoff.min(), clear_standoff.max(), num = len(Apollo_gtable))

peak_I_test = np.multiply(np.power(crs, slope), 10**intercept)

y_gauss = np.repeat(gaussmod[int(len(gaussmod)/2)::][:,np.newaxis], len(crs), 1) #use numpy broadcasting rules

theta_test = np.linspace(0,80, num=len(y_gauss))
theta_test = np.repeat(theta_test[:,np.newaxis], len(crs), 1)

test_impulses = np.multiply( y_gauss , peak_I_test )
crs = np.repeat(crs[:,np.newaxis], len(y_gauss), 1)
crs = crs.transpose()

def plot_model_surfaces():
    fig_test, ax = plt.subplots(1,1)
    fig_test.set_size_inches(3, 2.5)
    CS = ax.contourf(theta_test, crs, test_impulses, levels = [0 , 2.5, 5, 7.5, 10, 12.5, 15], cmap = plt.cm.cubehelix)
    cbar = fig_test.colorbar(CS)
    cbar.ax.invert_yaxis()
    cbar.ax.set_ylabel('peak specific impulse (MPa.ms)')
    ax.set_ylabel('standoff (clear charge radii)')
    ax.set_xlabel('incident wave angle (degrees)')
    plt.tight_layout()
    
    
    #Percentage change from CFD
    model_diff = abs(np.divide(np.subtract(test_impulses, peak_impulse/1e3), test_impulses)*100)
    fig_change, ax = plt.subplots(1,1)
    fig_change.set_size_inches(6, 5)
    CS = ax.contourf(theta_test, crs, model_diff, levels = [0,5,10,15,20,30,40,50,100], cmap = plt.cm.cubehelix)
    cbar = fig_change.colorbar(CS)
    cbar.ax.invert_yaxis()
    cbar.ax.set_ylabel('Relative change from CFD (%)')
    ax.set_ylabel('standoff (clear charge radii)')
    ax.set_xlabel('incident wave angle (degrees)')
    plt.tight_layout()
    fig_change.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Temp\surface_percentage.pdf"), format = 'pdf')
plot_model_surfaces()