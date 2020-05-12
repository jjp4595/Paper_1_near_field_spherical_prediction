# -*- coding: utf-8 -*-
"""
Impulse models
"""
import numpy as np
from MCEER_curves import MCEER
from scipy import ndimage

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
    R is  standoff (from centre of charge) - though this only dictates size of the target.
    theta_lim is the upper limit of theta for size of the circle (from bomb POV)
    theta is the theta coords for interpolation
    output is in Pa.s
    """
	#Plate dimensions
    A = R * np.tan(np.deg2rad(theta_lim)) * 2
    B = A
    res = 500
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



