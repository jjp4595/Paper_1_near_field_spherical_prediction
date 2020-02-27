import numpy as np


R = np.linspace(0.05, 1.5, 500)
R = np.add(R, 0.0246)
x= np.multiply(R, np.cos(np.deg2rad(45)))

with open("IncidentGauges.txt", "w+") as f:
    for i in range(500):   
        f.write("' " + str(i) + " ' " + str(round(x[i], 4)) + " " + str(round(x[i], 4)) + " " + str(0.001) +"\n")
f.close()    