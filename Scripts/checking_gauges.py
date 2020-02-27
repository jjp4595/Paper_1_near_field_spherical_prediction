"""
Quick checking script to plot gauge history vs time
"""
import preamble_functions as pre
import matplotlib.pyplot as plt #3.0.2


#Import Apollo data
Apollo_FileList = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\incident\*.txt")
Apollo_gtable = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\incident\*gtable",1)
Apollo_gauges = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\incident\*gauges",1)
# Apollo_FileList = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\reflected\*.txt")
# Apollo_gtable = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\reflected\*gtable",1)
# Apollo_gauges = pre.FileAddressList(r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\reflected\*gauges",1)


#Quick test to see if enough time in simulation
fig00, [ax_test1, ax_test2] = plt.subplots(1,2)
#Check first and last Pre/Imp gauge
ax_test1.plot(Apollo_gauges[0][:,0],Apollo_gauges[0][:,1])

ax_test2.plot(Apollo_gauges[0][:,0],Apollo_gauges[0][:,2])
