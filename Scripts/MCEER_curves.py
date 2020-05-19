# -*- coding: utf-8 -*-
"""
MCEER 
Coding up the MCEER polynomial charts
"""
import numpy as np
import matplotlib.pyplot as plt
import os

if __name__ == "__main__":
    
    params = {'font.family':'serif',
        'axes.labelsize':'small',
        'xtick.labelsize':'x-small',
        'ytick.labelsize':'x-small', 
        'lines.markersize': 2.5,
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
    
    #Graph ------------------------------------------------------------------------
    #polynomial constants for 0.0553 < Z < 40
    """
    I_s is kPa.ms/kg^(1/3)
    I_r is kPa.ms/kg^(1/3)
    """
    #Z = np.logspace(np.log10(0.0553), np.log10(40), 100)
    
    Z = np.logspace(np.log10(0.0553), np.log10(40), 100)
    
    U_i = 0.05671 + (0.8363*np.log10(Z))
    Y_i = (            2.21 +     -0.5499*U_i  +  -2.981*U_i**2 +     2.746*U_i**3 + 9.522*U_i**4 
         +  -13.01*U_i**5 +   -8.244*U_i**6  +   12.51*U_i**7 +   -0.5068*U_i**8 + 2.939*U_i**9
         +  1.023*U_i**10 +  -4.748*U_i**11  +  1.028*U_i**12 +   -3.996*U_i**13 + 1.703*U_i**14 
         + 0.5049*U_i**15 +  0.8065*U_i**16  +  2.006*U_i**17 +    -1.44*U_i**18 + 0.6624*U_i**19
         + -2.118*U_i**20 + 0.08177*U_i**21  +  1.678*U_i**22 +  -0.6299*U_i**23 + -0.1429*U_i**24
         + 0.07183*U_i**25)
    imp_i_scaled = 10**Y_i
    
    U_r = 0.2044 + (0.6997*np.log10(Z))
    Y_r = (  3.105      +     -2.032*U_r  +  1.36*U_r**2    +     -1.081*U_r**3 + -2.305*U_r**4 
                        +  3.143*U_r**5   +   1.151*U_r**6  +   -2.901*U_r**7 +   1.003*U_r**8)
    imp_r_scaled = 10**Y_r
    
    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(3, 2.5)
    ax.plot(Z, imp_i_scaled, 'k', label = 'incident impulse / $W^{1/3}$')
    ax.plot(Z, imp_r_scaled, 'k--',label = 'reflected impulse/ $W^{1/3}$')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('$Z = R/W^{1/3} (m/kg^{1/3})$')
    ax.set_ylabel('Specific impulse / $W^{1/3}$')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='center', bbox_to_anchor=(0.65, 0.8), prop={'size':6})
    plt.tight_layout()
    fig.savefig(os.path.join(os.environ['USERPROFILE'] + r"\Dropbox\Papers\Paper_1_near_field_spherical_prediction\Graphs\theta_peak_impulse_MCEER_curves.pdf"), format = 'pdf')

def MCEER(R, W):
    """
    Requires clear standoff, R and charge mass, W.
    Output is in kpa.ms
    """
    Wroot = W ** (1/3)
    Z = R / Wroot
    U_i = 0.05671 + (0.8363*np.log10(Z))
    
    Y_i = (            2.21 +     -0.5499*U_i  +  -2.981*U_i**2 +     2.746*U_i**3 + 9.522*U_i**4 
         +  -13.01*U_i**5 +   -8.244*U_i**6  +   12.51*U_i**7 +   -0.5068*U_i**8 + 2.939*U_i**9
         +  1.023*U_i**10 +  -4.748*U_i**11  +  1.028*U_i**12 +   -3.996*U_i**13 + 1.703*U_i**14 
         + 0.5049*U_i**15 +  0.8065*U_i**16  +  2.006*U_i**17 +    -1.44*U_i**18 + 0.6624*U_i**19
         + -2.118*U_i**20 + 0.08177*U_i**21  +  1.678*U_i**22 +  -0.6299*U_i**23 + -0.1429*U_i**24
         + 0.07183*U_i**25)
    imp_i_scaled = 10**Y_i
    imp_i = imp_i_scaled * Wroot
    
    U_r = 0.2044 + (0.6997*np.log10(Z))
    
    Y_r = (  3.105      +     -2.032*U_r  +  1.36*U_r**2    +     -1.081*U_r**3 + -2.305*U_r**4 
                        +  3.143*U_r**5   +   1.151*U_r**6  +   -2.901*U_r**7 +   1.003*U_r**8)
    imp_r_scaled = 10**Y_r
    imp_r = imp_r_scaled * Wroot
    return imp_i, imp_r




