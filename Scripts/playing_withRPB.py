# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 16:02:13 2020

@author: cip18jjp
"""


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
    return Ir * np.cos(np.deg2rad(theta))**2  + (Ii/Ir)*(1 + (np.cos(np.deg2rad(theta))**2 ) - 2*np.cos(np.deg2rad(theta)))

plt.figure()
theta = np.linspace(0,80,80)
frac = np.linspace(0,1,10)
for i in frac:
    plt.plot(theta,RPB(theta, 1, i))    
plt.plot(np.linspace(0,80,200), small['icr'].mean(1), 'r-', lw=2)
    

smoothRPB = np.asarray([large['RPB_MCEER_surf_1_2'][:,i]/ (max(large['RPB_MCEER_surf_1_2'][:,i])) for i in range(len(large['z']))]).T
plt.figure()
plt.plot(np.linspace(0,80,200), smoothRPB)
plt.plot(np.linspace(0,80,200), small['icr'].mean(1), 'k--', lw=2)



plt.figure() #Henrych
plt.plot(theta, np.cos(np.deg2rad(theta))**4)

#RPB
run = RPB_MCEER_i(0.1, 0.154, 80)
A = run[5][200,200::]/run[6][200,200::]
f = interp1d(run[2][200,200::], A)
A = f(np.linspace(0,80,80))

plt.figure()
plt.plot(np.linspace(0,80,200), small['icr'].mean(1), 'r-', lw=2)

#plt.plot(np.linspace(0,80,80), A)
plt.plot(np.linspace(0,80,80), RPB(np.linspace(0,80,80), 1, A.mean()))
#plt.plot(np.linspace(0,80,80), RPB(np.linspace(0,80,80),1, A))
plt.plot(run[2][200,200::],run[3][200,200::]/max(run[3][200,200::]))

#plt.plot(np.linspace(0,80,200),smoothRPB.mean(1), 'k--', lw = 2)    
runs=[]
for i in large['so']:
    runs.append(RPB_MCEER_i(0.1, i, 80))
As = [runs[i][5][200,200::]/runs[i][6][200,200::] for i in range(len(runs))]
thetaAs = np.asarray([runs[i][2][200,200::] for i in range(len(runs))]).T
plt.figure()
plt.plot(thetaAs, As)
plt.xlabel('theta')
plt.ylabel('A')
