"""
Quick checking script to plot gauge history vs time
"""
import preamble_functions as pre
import matplotlib.pyplot as plt #3.0.2


#Import Apollo data
Apollo_gauges = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\mesh_strategy\80mm_validation\*gauges",1)
Apollo_gauges2 = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Near Field Sims\Sims\Latest\80mm_with_afterburn\*gauges",1)


#Quick test to see if enough time in simulation



# for i in range(len(Apollo_gauges)):
#     fig00, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
#     ax1.set_xlabel(i)
#     #theta = 0
#     ax1.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,1]) #OP
#     ax2.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,201])#imp
    
#     #theta = 80
#     ax3.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,200])#OP
#     ax4.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,400])#imp

lim = 0.6
fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
ax1.plot(Apollo_gauges[0][:,0]*1000,Apollo_gauges[0][:,1]) #OP
ax2.plot(Apollo_gauges[0][:,0]*1000,Apollo_gauges[0][:,201])#imp
ax1.set_xlim(0,lim)
ax2.set_xlim(0,lim)
ax1.plot(Apollo_gauges2[0][:,0]*1000,Apollo_gauges2[0][:,1])
ax2.plot(Apollo_gauges2[0][:,0]*1000,Apollo_gauges2[0][:,101])    
#theta = 80
ax3.plot(Apollo_gauges[0][:,0]*1000,Apollo_gauges[0][:,200])#OP
ax4.plot(Apollo_gauges[0][:,0]*1000,Apollo_gauges[0][:,385])#imp
ax3.set_xlim(0,lim)
ax4.set_xlim(0,lim)