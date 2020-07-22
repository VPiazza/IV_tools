####################################################################################################################
##                                 ##                                     ##                                      ##
##  This Script is a v2.0 of the   ##          FUNCTIONING:               ##            ELECTRICAL MODEL:         ## 
##  SingleDiode.py. Smoother cod-  ##       - Enter Diode Parameters      ##      ---------------|\|\|\-----O     ## 
##  -ing and commenting, improved  ##       - Enter Voltage Settings      ##  I0 _|_      _|_      Rs             ##
##  functionality (plotting, fit-  ##       - Perform IV calculation      ##   n \ /       / Rsh                  ##
##  -ting) and more user-friendly. ##       - Save data in a file         ##      V       /__                     ##
##                                 ##                                     ##      |        |                      ##
##                                 ##                                     ##      --------------------------O     ##
####################################################################################################################          

#####################
# Importing Modules #
#####################

import math                         # For the IV calculation
import numpy as np                  # For the IV calculation
import matplotlib                   # For plotting
import matplotlib.pyplot as plt     # For plotting
import os                           # For navigate in the folders
import csvkit                       # For opening/writing .csv
import datetime                     # For recording the date            #NOT IMPLEMENTED YET


####################################
#         Script Functions         #
####################################
# Here below the  functions needed #
# for running the script.          #
####################################

def start_program():                ## While Loop: it defines the path of the exp data and it enables to optimize the result  
    answer = ""
    while answer != "yes" and answer != "no":
        ans = input("Do you want to fit data? (yes/no) ")
        answer = ans.strip().lower()
        if answer == "no":
            print("Bye bye")    
            return False
        elif answer == "yes":
            global data_path
            data_path = input("Enter the path where your data are stored: ")          
            # data_path = "C:\\Users\\valee\\VPScripts\\Work\\IV\\Tests"                  #This block can be used to speed up during testing
            return True
        else:
            print("Sorry, I didn't get it.")

#def wanna_numsteps():              ## !!OUTDATED!! It sets how to calculate the voltage range Now is set to: STEPS. 
    return True                      # >>Only for testing purpose. When test is over, re-activate the function<<
    # answer = ""
    # while answer != "steps" and answer != "points":
    #     ans = input("Do you want to set the voltage step or the number of acquisition points? (steps/points) ")
    #     answer = ans.strip().lower()
    #     if answer == "points":   
    #         return False
    #     elif answer == "steps":
    #         return True
    #     else:
    #         print("Sorry, I didn't get it.")

def modif_param(param, old_value):  ## Offers the possibility to skip parameter modification by inserting ""
    output = old_value
    msg = "Enter new " + param + ": "
    try:
        new_value = input(msg)
        if new_value != "":
            output = float(new_value)
    except:
        print("I didn't understand this entry")
        print("The old values of " + param + " is kept.")
    return output
        

####################################
#           Data Functions         #
####################################
# Here below the  functions needed #
# for:                             #
#   - extract data                 #
#   - estimate initial param.      #
####################################

def extract_data():                 ## Extracting form txt or csv. it returns "exp_V" and "exp_I" lists
    exp_V = []
    exp_I = []
    line_num = 0
    run_path = os.getcwd()
    namefile = input("Enter the name of the file: ")
    typefile = input("Which type of data? [txt/csv]")     
    header = int(input("How many lines need to be removed from the top? "))   
    # namefile = "test"                                     #This block can be used to speed up during testing
    # typefile = "csv"
    # header = 0
    os.chdir(data_path)
    if typefile == "csv":
        with open(namefile + ".csv") as file:
            rowreader = csvkit.reader(file, delimiter = ";")
            for row in rowreader:
                line_num+=1
                try:                                        #Verify if the this can act as a text removal.
                    exp_V.append(float(row[0]))
                    exp_I.append(float(row[1]))
                except:
                    print("Line " + str(line_num) + " contains NaNs")
            for i in (0, header):
                exp_V.pop(0)
                exp_I.pop(0)
        file.close()
    elif typefile == "txt":
        my_file = open(namefile + ".txt")
        data = my_file.read()
        my_file.close
        lines = data.split("\n")
        for line in lines:
            row = line.split("\t")
            line_num+=1
            try:
                exp_V.append(float(row[0]))
                exp_I.append(float(row[1]))
            except:
                print("Line " + str(line_num) + " contains NaNs")

    # voltage = np.array([volt])
    # current = np.array([curr])        # Might be necessary to perform a fitting. To be investigated.

    os.chdir(run_path)    
    return exp_V, exp_I

def param_estim(volt, curr):        ## Estimating starting param from the exp data. Might be fundamental for the fitting.
    v_min = min(volt)
    i_min = curr[volt.index(v_min)]
    v_max = max(volt)
    i_max = curr[volt.index(v_max)]
    return -i_min, 2.01, v_min/i_min, v_max/i_max  #I0, eta, Rsh, Rs


####################################
#    Math Functions & Constants    #
####################################
# Here below the  functions needed #
# for:                             #
#   - calculate IV curve           #
#   - (fit experimental data)      # TO BE IMPLEMENTED
####################################

def current_calc(voltage_range):    ## IV calculator: it requires the circuit elements and a V list; it returns the corresponding V and I lists.
    volt_fit = []
    curr_fit = []
    for Vd in voltage_range:
        Idf = diode_forw.I0 * (np.exp(Vd/(diode_forw.eta*Vt))-1)
        # Idr = -diode_rev.I0 * (np.exp(-Vd/(diode_rev.eta*Vt))-1)
        Ileak = Vd/resistance_circ.Rsh
        Itot = Idf + Ileak # + Idr
        Vtot = Vd + Itot*resistance_circ.Rs
        volt_fit.append(Vtot)
        curr_fit.append(Itot)
    return volt_fit, curr_fit

q = 1.60217662e-19  # C
Kb = 1.38e-23  # JK^-1
T = 298 #K (25°C)
Vt =Kb*T/q # V

####################################
#               Objects            #
####################################
#  Here below the objects which    #
#  define the circuit elements.    #
####################################

                                    ## All the properties in these classes must be float.
class Diode():                      ## Saturation Current and ideality factor
    def __init__(self, I0, eta):
        self.I0 = I0
        self.eta = eta

class Resistances():                ## Shunt and Series resistances
    def __init__(self, Rsh, Rs):
        self.Rsh = Rsh
        self.Rs = Rs

class VoltageSettings():            ## voltage range definition. It returns voltage_range list.                                     
    def __init__(self, V_low, V_high, V_step):         
        self.V_low = V_low
        self.V_high = V_high
        self.V_step = V_step
        self.V_points = ((V_high-V_low)/V_step) + 1   # n° points includes also V_start and V_end
    
    def voltage_range(self):
        voltage_range = []
        volt = self.V_low
        for i in range(0, int(self.V_points)):
            volt = self.V_low + self.V_step * i
            voltage_range.append(volt)
        return voltage_range

#######################################################################################################
#####                                        RUNNING SCRIPT                                       #####
while start_program():
    # Using experimental data to define the initial value of the circuit parameters
    exp_V, exp_I = extract_data()
    I0,eta,Rsh,Rs = param_estim(exp_V, exp_I)
    my_voltage = VoltageSettings(min(exp_V),max(exp_V), 0.001)
    diode_forw=Diode(I0,eta)   
    resistance_circ=Resistances(Rsh, Rs)  
    # my_voltage = VoltageSettings(-2,2,0.001)                                 #This block can be used to speed up during testing
    # diode_forw=Diode(1e-12,2)    
    # resistance_circ=Resistances(1e12, 1e6)                          

    # Defining the voltage settings of the simulation
    conf1 = ""                                    
    while conf1 != "y":
        print("The simulation will be performed from " + str(my_voltage.V_low) + " V to " + str(my_voltage.V_high) 
              + " V for every " + str(my_voltage.V_step) + " V.")
        conf1 = input("Are the voltage settings correct? [y/n]: ").strip().lower()
        if conf1 != "y":                                                        
            print("Previous parameters: ")
            print("Start = " + str(my_voltage.V_low) + " V")
            print("End = " + str(my_voltage.V_high) + " V")
            print("Step = " + str(my_voltage.V_step) + " V")
            my_voltage = VoltageSettings(modif_param("start", my_voltage.V_low),modif_param("end", my_voltage.V_high), 
                                         modif_param("step", my_voltage.V_step))
           
    voltage_range = my_voltage.voltage_range()

    # Visual optimization of the fit (Iterate plotting) 
    conf2 = ""
    while conf2 != "y":

        # Simulate the IV curve from the circuit model
        volt_fit, curr_fit = current_calc(voltage_range)
        ################################## FITTING METHOD 2 #######################################
        ##### Another option can be to use mesh grids to minimize the equation. To be studied.#####
        # voltage_space = np.linspace(exp_V[0], exp_V[len(exp_V)-1],len(exp_V))
        # current_space = np.linspace(exp_I[0], exp_I[len(exp_I)-1], len(exp_I))
        # voltage_space, current_space = np.meshgrid(voltage_space, current_space)
        # I_eq = lambda V, i: i - (diode_forw.I0 * ((np.exp(V-i*resistance_circ.Rs) / diode_forw.eta * Vt))-1) - ((V-i*resistance_circ.Rs) / resistance_circ.Rsh)
        
        # Plotting experimental and simulated curves
        fig, ax = plt.subplots(figsize=(8,8))
        # CS1 = ax.contour(voltage_space, current_space, I_eq(voltage_space, current_space), [0])
        CS1 = ax.plot(volt_fit, curr_fit, "r-")
        CS2 = ax.plot(exp_V, exp_I, "bo")
        ax.set_xlabel('Voltage [V]')
        ax.set_ylabel('Current [A] ')
        ax.set_xlim(min(exp_V)*1.1, max(exp_V)*1.1)
        ax.set_ylim(min(exp_I)*10, max(exp_I)*1.1)
        ax.set_title('IV Curve')
        # ax.clabel(CS1, fontsize=10, fmt='Fitting')
        # ax.clabel(CS2, fontsize=10, fmt='Experimental')
        # plt.legend('Experimental')
        plt.show()
        conf2 = input("Is the fitting good enough? [y/n]: ").strip().lower()
        if conf2 != "y":                                                        
            print("Previous parameters: ")
            print("I0 = " + str(diode_forw.I0) + " A")
            print("eta = " + str(diode_forw.eta))
            print("Rsh = " + str(resistance_circ.Rsh) + " Ohm")
            print("Rs = " + str(resistance_circ.Rs) + " Ohm")
            diode_forw = Diode(modif_param("I0", diode_forw.I0),modif_param("eta", diode_forw.eta))
            resistance_circ = Resistances(modif_param("Rsh", resistance_circ.Rsh),modif_param("Rs", resistance_circ.Rs))
    
    # Saving results:it creates a folder in a given path and write the parameters and simulated IV in a csv file
    namefile = input("Enter the name of the file: ")              
    namefolder = "IV_Simulation"
    path = input("Enter the path to save the result: ") 
    # namefile = "TestResults_"                                 #This block can be used to speed up during testing
    # namefolder = "IV_Simulation"
    # path = data_path
    os.chdir(path)
    try:
        os.mkdir(namefolder)
        os.chdir(namefolder)
    except:
        os.chdir(namefolder)
    with open(namefile + ".csv","w") as file:
        writer = csvkit.writer(file, delimiter= ";")
        writer.writerows([["I0 [A]", "eta", "Rs [Ohm]", "Rsh [Ohm]"],[diode_forw.I0, diode_forw.eta, resistance_circ.Rs, resistance_circ.Rsh]])
        writer.writerows([["Voltage","Current"],["V", "A"]])
        for volt in volt_fit:
            writer.writerow([str(volt), str(curr_fit[volt_fit.index(volt)])])
    file.close()









