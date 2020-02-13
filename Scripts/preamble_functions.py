import glob
import numpy as np #1.15.4
import re
import scipy
#Importing Data ----------------------------------------------------------------------------------------------------------
def FileAddressList(fileIN, GaugeData=None):
    #This functions reads the fileIN and then creates the filelist as a list of strings from that recursive search
    #If a second argument is entered, then it will import the gauge data. 
    #In other words, one input for file lists, add a number for gauge data import.
    OutputList=[]
    if GaugeData is None:
        for filename in glob.glob(fileIN, recursive=True):
            OutputList.append(filename)
    else: 
        for filename in glob.glob(fileIN, recursive=True):
            OutputList.append(np.loadtxt(filename))
    return OutputList

#Defining Additional Functions ------------------------------------------------------------------------------------------------------
def ElementSizes(fileIN): #This calculates Element Size to plot on scatter graphs.
    with open(fileIN) as f:
        content = f.readlines()
        content = [x.strip() for x in content] #remove whitespace at end of each line. 
    model_info = np.core.defchararray.find(content, "model resolution level, zone length")
    index = np.argwhere(model_info > 0)
    a=[]
    for word in content[index[0,0]].split():
        try:
            a.append(float(word))
        except ValueError:
            pass
    ElementSize = a[1] / (2 ** a[0])
    return ElementSize

def CPUFinder(fileIN):#This is the CPUFinder function. It returns an integer value for wall time. 
    with open(fileIN) as f:
        content = f.readlines()
        content = [x.strip() for x in content] #remove whitespace at end of each line. 
        index_wall = np.core.defchararray.find(content, "Wall")
        newcontent=[]
        wall_times=0
        for newindex in enumerate(np.argwhere(index_wall > 0)):
            newcontent.append(content[newindex[1][0]])
            wall_times = wall_times + (int(re.search(r'\d+', newcontent[newindex[0]]).group()))
    return wall_times

def MaxOP(fileIN):
    #This function works for imported Apollo Data.
    maxOPlist=[]
    for i in range(len(fileIN)):#Loop to calculate maxOP
        maxOPlist.append((np.amax(fileIN[i][:,1])-101300)/1000)#Gauge of interest always the first one.
    return maxOPlist

def MaxI(fileIN):
    #This function works for imported Apollo Data.
    maxIlist=[]
    for i in range(len(fileIN)):#Loop to calculate maxOP
        maxIlist.append((np.amax(fileIN[i][:,3])))#Gauge of interest always the first one.
    return maxIlist

def PeakImpDist(fileIN, startCol, endCol):
    peakImpDistList=[]
    numCols  = endCol - startCol
    for i in range(numCols):
        peakImpDistList.append(np.amax(fileIN[:,i + startCol]))
    return peakImpDistList

#Calculating scaled distance
def scaled_dist(fileIN): #This calculates scaled distance by subtracting the charge center from y domain max, so this may have to change. 
    with open(fileIN) as f:
        content = f.readlines()
        content = [x.strip() for x in content] #remove whitespace at end of each line. 
    #Getting charge info
    charge_info = np.core.defchararray.find(content, "use charge, main")
    index_charge = np.argwhere(charge_info > 0)
    charge_mass=[]
    for word in content[index_charge[0,0]].split():
        try:
            charge_mass.append(float(word))
        except ValueError:
            pass
   #Getting Standoff info
    chargecenter_info = np.core.defchararray.find(content, "charge center")
    index_chargecenter = np.argwhere(chargecenter_info > 0)
    chargecenter=[]
    for word in content[index_chargecenter[0,0]].split():
        try:
            chargecenter.append(float(word))
        except ValueError:
            pass
   #Getting edge boundary info
    bound_index = 24 #Change this if row that model info is on changes!
    bound_info=[]
    for word in content[bound_index].split():
        try:
            bound_info.append(float(word))
        except ValueError:
            pass
    bound = bound_info[4] #This time we are interested in y domain. 
    standoff = bound - chargecenter[1] 
    z = np.divide(standoff, np.power(charge_mass[0], (1/3))) #Scaled Distance
    return z

#Calculating scaled distance
def standoff_func(fileIN): #This calculates scaled distance by subtracting the charge center from y domain max, so this may have to change. 
    with open(fileIN) as f:
        content = f.readlines()
        content = [x.strip() for x in content] #remove whitespace at end of each line. 
    #Getting charge mass
    charge_info = np.core.defchararray.find(content, "use charge, main")
    index_charge = np.argwhere(charge_info > 0)
    charge_mass=[]
    for word in content[index_charge[0,0]].split():
        try:
            charge_mass.append(float(word))
        except ValueError:
            pass
   #Getting Standoff info
    chargecenter_info = np.core.defchararray.find(content, "charge center")
    index_chargecenter = np.argwhere(chargecenter_info > 0)
    chargecenter=[]
    for word in content[index_chargecenter[0,0]].split():
        try:
            chargecenter.append(float(word))
        except ValueError:
            pass
   #Getting edge boundary info
    bound_index = 24 #Change this if row that model info is on changes!
    bound_info=[]
    for word in content[bound_index].split():
        try:
            bound_info.append(float(word))
        except ValueError:
            pass
    bound = bound_info[4] #This time we are interested in y domain. 
    standoff = bound - chargecenter[1] 
    return standoff

#Guassian Function of multiple curves--------------------------------------------------
def Gaussian_func(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)**2 / (wid**2) ) )
    return y

#Gaussian-2 Curve
def gauss_two_curve(x, *params): #no asterix here as function supplied as a guess anyway
    y1 = np.zeros_like(x)
    y2 = np.zeros_like(x)
    y1 = params[1] * np.exp( -((x - params[0])**2 / (params[2]**2) ) )
    y2 = params[4] * np.exp( -((x - params[3])**2 / (params[5]**2) ) )
    y = y1 + y2
    return y

#Generalized Bell Function ---------------------------------------------------------
def bell_gen(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 4): # http://researchhubs.com/post/maths/fundamentals/bell-shaped-function.html 
        a = params[i] #width
        b = params[i+1] #flatness
        c = params[i+2] #centre
        d = params[i+3] #Scaling Factor
        y = d * (1 / (1 + abs( (x - c) / a)**(2*b)))
    return y

#Trig Function from Jang Paper --------------------------------------------------------------------
def Jang_func(x, z, I, *params):
    #z is scaled distance
    #https://www.un.org/disarmament/un-saferguard/kingery-bulmash/ is the calculator for initial conditions
    theta = np.zeros_like(x)
    beta = np.zeros_like(x)
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        a = params[i]
        b = params[i]
        c = params[i]
        d = 0.2 #Don't know if this is right. 20% of the limiting value of beta which can be 1.0 at its max.
        theta = np.arctan((0.5*x) / z)
        beta = a * (1 + b*(z-c)) * np.exp(-b*(z-c)) + d # the weighting function
        y = I * np.exp(-beta * np.tan(theta) * np.tan(theta))
    return y

#Performance Indicator Function -------------------------------------------------------------------
#RMS Error Chosen
def rmse(predictions, actual):
    return np.sqrt(((predictions-actual)**2).mean())

#Normalising function
def normalise_func(X):
    z = ( X[:] - np.min(X) ) / (np.max(X) - np.min(X))
    return z

#Amplitude function
def exp_func(x, a, tau, c):
    return a * np.exp(-x / tau) + c 
