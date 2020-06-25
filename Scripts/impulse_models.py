 # -*- coding: utf-8 -*-
"""
Impulse models
"""
import numpy as np
from MCEER_curves import MCEER
from scipy import ndimage

#Models in the literature: ---------------------------------------------------
#1) Henrych/Rem Model
def Henrych_i(a, r0, D0, Q0, rho, theta):
    """

    Parameters
    ----------
    a : TYPE
        standoff distance (m)
    r0 : TYPE
        charge radius (m).
    D0 : TYPE
        detonation velocity (m/s).
    Q0 : TYPE
        Specific Energy (MJ/kg).
    rho : TYPE
        Density (kg/m3).
    theta : TYPE
        .
    Returns
    -------
    i : array
        specific impulse distribution.

    """
    P0 = rho * (D0**2) /2 / (3+1) #detonation pressure, Pa
    u0 = (2*Q0*1000000)**0.5 #detonation products velocity, m/s
    w0 = P0/rho/u0 #velocity of hypothetical surface, m/s
    A0 = (u0 + w0)/(4*np.pi) #explosive constant m/s
    
    tau = (1/u0 + 1/w0) * r0 #DP loading duration (s)
    Pm = P0 * (r0/a)**2 #Max pressure at barrier
    
    W = (4/3) * np.pi * (r0**3) * rho
    
    i = ((A0 * W) / (a**2)) * np.cos(np.deg2rad(theta))**4
    
    return i

# def Henrych_I(A0, W, theta_lim):
#     """
#     Calculating total impulse across a circle
#     """
#     return np.pi * A0 * W * np.sin(np.deg2rad(theta_lim))**2

# def Henrych_I(a, r0, D0, Q0, rho, theta_lim):
#     """

#     Parameters
#     ----------
#     a : TYPE
#         standoff distance (m)
#     r0 : TYPE
#         charge radius (m).
#     D0 : TYPE
#         detonation velocity (m/s).
#     Q0 : TYPE
#         Specific Energy (MJ/kg).
#     rho : TYPE
#         Density (kg/m3).
#     theta : TYPE
#         .
#     Returns
#     -------
#     i : array
#         specific impulse distribution.

#     """
#     P0 = rho * (D0**2) /2 / (3+1) #detonation pressure, Pa
#     u0 = (2*Q0*1000000)**0.5 #detonation products velocity, m/s
#     w0 = P0/rho/u0 #velocity of hypothetical surface, m/s
#     A0 = (u0 + w0)/(4*np.pi) #explosive constant m/s
    
#     tau = (1/u0 + 1/w0) * r0 #DP loading duration (s)
#     Pm = P0 * (r0/a)**2 #Max pressure at barrier
    
#     W = (4/3) * np.pi * (r0**3) * rho
    
#     I = np.pi * A0 * W * np.sin(np.deg2rad(theta_lim))**2
    
#     return I

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
    res = 401
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
    return [X,Y, theta_MCEER, imp, I, is_max, ir_max]


#3) CFD Impulse generation for circular target --------------------------------
def Impulse_CFD(specific_impulse, R, theta_lim, theta):
    """
    i is impulse distribution to be interpolated
    R is  standoff (from centre of charge) - though this only dictates size of the target.
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

def TotalImpulseCalc(specific_impulse, R, theta_lim, theta):
    """
    i is impulse distribution to be interpolated
    R is  standoff (from centre of charge) - though this only dictates size of the target.
    theta_lim is the upper limit of theta for size of the circle (from bomb POV)
    theta is the theta coords for interpolation
    output is in Pa.s
    """
	#Plate dimensions
    A = R * np.tan(np.deg2rad(theta_lim))
    B = A
    res = 400
    x = np.linspace(0, A, res)
    y = np.linspace(0, B, res)
    [X,Y] = np.meshgrid(x,y)
    theta_array = np.rad2deg(np.arctan(np.divide(np.power(np.add(np.power(X,2), np.power(Y,2)), 0.5), R))) 
	
    ir = np.interp(theta_array, theta, specific_impulse) 

	#Rectangular area
    Area = np.ones_like(ir)
    Area = Area * (x[1]-x[0]) * (y[1]-y[0])

    
    #Circular mask over area for circular target 
    r = res
    cen_x, cen_y = 0,0
    y1, x1 = np.ogrid[:res, :res]
    mask = x1 * x1 + y1 * y1 <= r*r
    
    array = np.zeros((res, res))
    array[mask] = 1
    struct = ndimage.generate_binary_structure(2, 2)
    erode = ndimage.binary_erosion(mask, struct)
    edges = mask ^ erode
    array[edges] = 0.5
    array[:,0], array[0,:] = 1,1
    
    Area = np.multiply(array, Area)
    
    Imp = np.multiply(ir, Area)
    I = np.sum(Imp) * 4
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
def RPB_MCEER(theta, so, W):
    R = so / np.cos(np.deg2rad(theta))
    RPB = np.asarray(MCEER(R, W)).T
    RPB = RPB[:,1]*np.cos(np.deg2rad(theta))**2 + RPB[:,0]*(1+np.cos(np.deg2rad(theta))**2 - 2*np.cos(np.deg2rad(theta)))
    return RPB

# def RPB_MCEER_norm(theta, so, W):
#     R = so / np.cos(np.deg2rad(theta))
#     RPB = np.asarray(MCEER(R, W)).T
#     RPB = np.cos(np.deg2rad(theta))**2 + (RPB[:,0]/RPB[:,1])*(1+np.cos(np.deg2rad(theta))**2 - 2*np.cos(np.deg2rad(theta)))
#     return RPB

# theta = np.linspace(0,80,200)
# so = 0.1
# W = 0.12
# imp1 = RPB_MCEER(theta, so, W)
# imp1 = imp1/imp1.max()
# imp2 = RPB_MCEER_norm(theta, so, W)

# plt.plot(theta,imp1)
# plt.plot(theta,imp2)
# plt.legend(['imp1 =  Normalised by ir_max', 'imp2 Normalised by ir'])