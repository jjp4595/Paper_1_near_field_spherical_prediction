"""
Quick checking script to plot gauge history vs time
"""
import preamble_functions as pre
import matplotlib.pyplot as plt #3.0.2


#Import Apollo data
# Apollo_FileList = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\incident\*.txt")
# Apollo_gtable = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\incident\*gtable",1)
# Apollo_gauges = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\incident\*gauges",1)
# Apollo_FileList = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\reflected\*.txt")
# Apollo_gtable = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\reflected\*gtable",1)
# Apollo_gauges = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\reflected\*gauges",1)
Apollo_gauges = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\NewFolder1\*gauges",1)


#Quick test to see if enough time in simulation



for i in range(len(Apollo_gauges)):
    fig00, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2,2)
    ax1.set_xlabel(i)
    #theta = 0
    ax1.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,1]) #OP
    ax2.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,201])#imp
    
    #theta = 80
    ax3.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,200])#OP
    ax4.plot(Apollo_gauges[i][:,0],Apollo_gauges[i][:,400])#imp