import matplotlib.pyplot as plt
import numpy as np
import math

def cutoff_frequency(T_debeye):
    #Boltzmann constant in J/K
    kb = 1.381 * 10**(-23)
    #Reduced Plank's constant in J*s
    hbar = 1.055 * 10**(-34)
    #Formula for debeye cutoff frequency
    debeye_cutoff = (kb*T_debeye)/(hbar)
    return debeye_cutoff

def g(w,debeye_cutoff):
    #Equation 2.5 in the book for one mole of atoms 
    return (6.022*10**23)* (9*(w**2))/(debeye_cutoff**3)

def bose_occupation_factor(w, T):
    #Boltzmann constant in J/K
    kb = 1.381 * 10**(-23)
    #Reduced Plank's constant in J*s
    hbar = 1.055 * 10**(-34)
    #Definition of Bose occupation factor
    try:
        return 1/(math.e**((hbar*w)/(kb*T))-1)
    except:
        return 0

def energy_expectation_value(w, T, debeye_cutoff):
    #Reduced Plank's constant in J*s
    hbar = 1.055 * 10**(-34)
    #Equation 2.8 in the book 7
    return hbar*w*g(w,debeye_cutoff)*bose_occupation_factor(w, T)

#Implement Simpson's Rule to calculate the integral
#We update xi by h for each iteration
#If even, multiply h/3*f(xi) by 2. If odd, multiply by 4
def simpsons(b,a,n, debeye_cutoff, T):
    h = (b-a)/n 
    I = (h/3)*(energy_expectation_value(a, T, debeye_cutoff) + energy_expectation_value(b, T, debeye_cutoff))
    k = 1
    xi = a + h
    n = n-1 
    while(n > 0):
        if(k%2 == 0):
            I += 2*(h/3)*energy_expectation_value(xi, T, debeye_cutoff)
        elif(k%2 == 1):
            I += 4*(h/3)*energy_expectation_value(xi, T, debeye_cutoff)
        n -= 1
        k += 1
        xi += h
    return I

#Function to calculate a derivative (For a central derivative, dividing by h instead of 2h as x2 is to the right and x1 is to the left of our data point so x2-x1 = 2h
def Derivative(E1, T1, E2, T2):
#Step value between T2 and T1
    h = T2 - T1
#Formula for derivative
    df = (E2-E1)/(h)
    return df

#Finding Debye Heat Capacity by first finding <E> and then d<E?/DT
def Debye_Heat_Capacity(T_debye):
    #Max Frequency
    frequency_max = cutoff_frequency(T_debye)

    #Storing energy expectation value energy
    E_Data = []
    #100 equally spaced temperature values between 0K and the debye temperature 
    T_Data = range(0, T_debye, 1)
    length = len(T_Data)
    #Calculating expectation value of energy at different temperature values
    for ii in range(0,length):
        #Equation 2.8 for expectation value of energy (Multiplying by 10^23 to do it in units per mol)
        E_Data.append(simpsons(frequency_max, 0, 1000, frequency_max, T_Data[ii]))

    #Creating variable to store heat capacity
    C_Data = []

    #Heat capacity is DE/DT so calculating this
    for ii in range(0,length):
    #If this is the first element, we need to take a forward derivative on this element (Push in the first element and second element into derivative function)
        if(ii == 0):
            C_Data.append(Derivative(E_Data[0], T_Data[0],E_Data[1], T_Data[1]))
    #If this is a central element, we take a central derivative (Push in last element and next element into derivative function)
        elif(ii == length-1):
            C_Data.append(Derivative(E_Data[length-2], T_Data[length - 2],E_Data[length-1], T_Data[length-1]))
    #If this is the last element, we take a backward derivative (Push in last element and second-last element into derivative function)
        else:
            C_Data.append(Derivative(E_Data[ii-1], T_Data[ii-1],E_Data[ii+1], T_Data[ii+1]))

    return C_Data

T_debye_nickel = 345
T_debye_strontium = 148
T_debye_cadmium = 221
T_debye_cerium = 138
T_debye_erbium = 163
T_debye_platinum = 225

T_Data_Nickel = list(range(0, T_debye_nickel, 1))
T_Data_Strontium = list(range(0, T_debye_strontium, 1))
T_Data_Cadmium = list(range(0, T_debye_cadmium, 1))
T_Data_Cerium = list(range(0, T_debye_cerium, 1))
T_Data_Erbium = list(range(0, T_debye_erbium, 1))
T_Data_Platinum = list(range(0, T_debye_platinum, 1))

C_Data_Nickel = Debye_Heat_Capacity(T_debye_nickel)
C_Data_Strontium = Debye_Heat_Capacity(T_debye_strontium)
C_Data_Cadmium = Debye_Heat_Capacity(T_debye_cadmium)
C_Data_Cerium = Debye_Heat_Capacity(T_debye_cerium)
C_Data_Erbium = Debye_Heat_Capacity(T_debye_erbium)
C_Data_Platinum = Debye_Heat_Capacity(T_debye_platinum)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(T_Data_Nickel, C_Data_Nickel, "-r")
ax.set_title("Heat Capacity of Nickel (Debye Model)")
ax.set_xlabel("T (K)")
ax.set_ylabel("Heat Capacity (J/(Mol*K))")
fig.savefig("Heat_Capacities\\Nickel.png")

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(T_Data_Strontium, C_Data_Strontium, "-b")
ax2.set_title("Heat Capacity of Strontium (Debye Model)")
ax2.set_xlabel("T (K)")
ax2.set_ylabel("Heat Capacity (J/(Mol*K))")
fig2.savefig("Heat_Capacities\\Strontium.png")

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(T_Data_Cadmium, C_Data_Cadmium, "-g")
ax3.set_title("Heat Capacity of Cadmium (Debye Model)")
ax3.set_xlabel("T (K)")
ax3.set_ylabel("Heat Capacity (J/(Mol*K))")
fig3.savefig("Heat_Capacities\\Cadmium.png")

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.plot(T_Data_Cerium, C_Data_Cerium, "-c")
ax4.set_title("Heat Capacity of Cerium (Debye Model)")
ax4.set_xlabel("T (K)")
ax4.set_ylabel("Heat Capacity (J/(Mol*K))")
fig4.savefig("Heat_Capacities\\Cerium.png")

fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
ax5.plot(T_Data_Erbium, C_Data_Erbium, "-y")
ax5.set_title("Heat Capacity of Erbium (Debye Model)")
ax5.set_xlabel("T (K)")
ax5.set_ylabel("Heat Capacity (J/(Mol*K))")
fig5.savefig("Heat_Capacities\\Erbium.png")

fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
ax6.plot(T_Data_Platinum, C_Data_Platinum, "-m")
ax6.set_title("Heat Capacity of Platinum (Debye Model)")
ax6.set_xlabel("T (K)")
ax6.set_ylabel("Heat Capacity (J/(Mol*K))")
fig6.savefig("Heat_Capacities\\Platnium.png")

dulong_petit = []
length = len(T_Data_Nickel)
for ii in range(0, length):
    dulong_petit.append(24.9)

fig7 = plt.figure()
ax7 = fig7.add_subplot(111)
ax7.plot(T_Data_Nickel, C_Data_Nickel, "-r", label = "Nickel")
ax7.plot(T_Data_Strontium, C_Data_Strontium, "-b", label = "Strontium")
ax7.plot(T_Data_Cadmium, C_Data_Cadmium, "-g", label = "Cadmium")
ax7.plot(T_Data_Cerium, C_Data_Cerium, "-c", label = "Cerium")
ax7.plot(T_Data_Erbium, C_Data_Erbium, "-y", label = "Erbium")
ax7.plot(T_Data_Platinum, C_Data_Platinum, "-m", label = "Platinum")
ax7.plot(T_Data_Nickel, dulong_petit, "-k", label = "Dulong-Petit")
ax7.set_title("Heat Capacity of Varius Elements (Debye Model)")
ax7.set_xlabel("T (K)")
ax7.set_ylabel("Heat Capacity (J/(Mol*K))")
ax7.legend()
fig7.savefig("Heat_Capacities\\Combined.png")


def scaled_temp(T_debye, T_data):
    length = len(T_data)
    for ii in range(0, length):
        T_data[ii] = T_data[ii]/T_debye
    return T_data

T_Data_Nickel = scaled_temp(T_debye_nickel, T_Data_Nickel)
T_Data_Strontium = scaled_temp(T_debye_strontium, T_Data_Strontium)
T_Data_Cadmium = scaled_temp(T_debye_cadmium, T_Data_Cadmium)
T_Data_Cerium = scaled_temp(T_debye_cerium, T_Data_Cerium)
T_Data_Erbium = scaled_temp(T_debye_erbium, T_Data_Erbium)
T_Data_Platinum = scaled_temp(T_debye_platinum, T_Data_Platinum)

fig8 = plt.figure()
ax8 = fig8.add_subplot(111)
ax8.plot(T_Data_Nickel, C_Data_Nickel, "-r", label = "Nickel")
ax8.plot(T_Data_Strontium, C_Data_Strontium, "-b", label = "Strontium")
ax8.plot(T_Data_Cadmium, C_Data_Cadmium, "-g", label = "Cadmium")
ax8.plot(T_Data_Cerium, C_Data_Cerium, "-c", label = "Cerium")
ax8.plot(T_Data_Erbium, C_Data_Erbium, "-y", label = "Erbium")
ax8.plot(T_Data_Platinum, C_Data_Platinum, "-m", label = "Platinum")
ax8.plot(T_Data_Nickel, dulong_petit, "-k", label = "Dulong-Petit")
ax8.set_title("Heat Capacity of Varius Elements (Debye Model)")
ax8.set_xlabel("T/Td")
ax8.set_ylabel("Heat Capacity (J/(Mol*K))")
ax8.legend()
fig8.savefig("Heat_Capacities\\Combined_and_Scaled.png")









