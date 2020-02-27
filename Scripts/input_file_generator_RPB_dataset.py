#Some Apollo prep functions and eqs.
#Importing packages
import numpy as np
import os

#------------------------------------------------------------------------------
#meta variables 
no_exps = 2


#Charge info
mass = 0.1
charge_rad = 0.0246
tnt_eq = 1
shape_ratio = 0 #0 for sphere

#Stages Variables
term_time = 0.03

#Model
res_level = 5
zone_length = 0.05



#Standoffs/Scaled Distances
ymax = np.linspace(0.055, 1.5, no_exps)
ymax = np.round(ymax, decimals = 4)

#Output and theta range
no_gauges = 1
no_plot_files = 20

            
#file location information
batch_name = "new_batch.bat"
batch_path = r"D:\PhD projects\Pannell_Jordan\RPB\reflected"
template = "RPB_reflected_template.txt"
local_path = r"C:\Users\jorda\Google Drive\Apollo Sims\Impulse Distribution Curve Modelling\Paper_1\RPB\reflected"


def myround(x, base):
    """
    Rounding function
    """
    return base * round(x/base)


#------------------------------------------------------------------------------   
def reflected_file_creator(mass, tnt_eq, shape_ratio, term_time, res_level, zone_length, ymax, no_gauges, no_plot_files, batch_name, batch_path, template, local_path):
    """
    The function below creates the input files and batch file from a given template file. More paramaters can be added that need to be changed. 
    """    
    startfile = 100
    file_list = [str(num)+".txt" for num in list(range(startfile, startfile + no_exps))]
    
    batch_lines = []
    batch_lines.append("cd " + batch_path)
   
    #Sets local path to create files
    os.chdir(local_path)

    
    #Start looping through input files here
    for i in range(len(ymax)):
        #Create batch file lines
        batch_lines.append("blastsimulator_win.exe " + file_list[i])
               
        
        x_max = ymax[i]
        z_max = ymax[i]
                         

        
        #Open template file and strip content
        with open(template) as f1:
            content = f1.readlines()
            content = [x.strip() for x in content]
        
        #Create new input file to write to
        with open(file_list[i], "w") as f2:
            
            #Model----------------
            str_to_search = "\t\t\t\t\t...model resolution level, zone length"
            str_index = np.argwhere(np.core.defchararray.find(content, str_to_search) > 0)
            content[int(str_index + 3)] = "Block 'Part1' 0 0 0  " + str(x_max) + " " + str(ymax[i]) + " " + str(z_max) + " " + " 1 1 1  3 1 3"                  
            
            #Output---------------
            str_to_search = "\t\t\t\t\t...number of gauges"
            gauges_index = np.argwhere(np.core.defchararray.find(content, str_to_search) > 0)
            content[int(gauges_index)] = str(no_gauges) + str_to_search 
                        
            reflected_gauge_string = "' " + str(0) + " ' " + str(0.0001) + " " + str(round(ymax[i]- 0.001, 4)) + " " + str(0.0001)
            reflected_index = int(gauges_index + 1)
            content.insert(reflected_index, reflected_gauge_string)

            #Final step -- write lines of amended content to new file
            f2.write('\n'.join(content) + '\n')
            
        # add lines to delete misc files 
        batch_lines.append("DEL *.mod")
        batch_lines.append("DEL *.vtk")
        batch_lines.append("DEL *_info")
        batch_lines.append("DEL *_plots")
        batch_lines.append("DEL *_status")
        batch_lines.append("DEL *_parts")
        batch_lines.append("DEL *_forces")
        batch_lines.append("DEL *_map") 
        batch_lines.append("DEL *_b1d.dat")
        batch_lines.append("DEL *_effects")
        #batch_lines.append("DEL *_gauges")
    
    #write batch file
    with open(batch_name, "w") as b:
        b.write('\n'.join(batch_lines) + '\n')

#------------------------------------------------------------------------------
reflected_file_creator(mass, tnt_eq, shape_ratio, term_time, res_level, zone_length, ymax, no_gauges, no_plot_files, batch_name, batch_path, template, local_path)

