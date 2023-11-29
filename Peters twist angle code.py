import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import interpolate


b = 66.45
c_r = 10.22960376
taper = 0.28
G = 26*10**9


#chord
def chord(y):
    chord = c_r - c_r * (1-taper) / (b/2) * y
    return chord

#Spar Area
area_spar = float(input("Enter the area of spar point areas: "))

#spar height
def left_sparheight(y):
    left_sparheight = 0.1266068 * chord(y)
    return left_sparheight

def right_sparheight(y):
    right_sparheight = 0.11369081 * chord(y)
    return right_sparheight

#Get spar thickness
def left_spart (y) :
    left_spart = 2*area_spar / left_sparheight(y)
    return left_spart

def right_spart (y) :
    right_spart = 2*area_spar / right_sparheight(y)
    return right_spart


#Get top and bottom thickness
top_bottom_thickness = float(input("what is the thickness of the top and bottom surfaces? "))

#Polar Moment of Inertia Function 
def J(y):
    sum_s_t = (left_sparheight(y) / left_spart(y)) + (right_sparheight(y) / right_spart(y)) +  (0.40019 + 0.4) * chord(y)/top_bottom_thickness
    internal_A = (right_sparheight(y) + left_sparheight(y))/2 * 0.4 * chord(y)
    J_value = 4*internal_A**2/sum_s_t
    return J_value



#Stefan's code
shear_lst = []
torque_lst = []
moment_lst = []

"""
sw = 0 for positive load factor, maximum shear case
sw = 1 for positive load factor, maximum bending stress case
sw = 2 for negative load factor, maximum shear and maximum bending stress case
"""
sw = 0


if sw == 0:
    with open('positive_OEW+Payload.txt', 'r') as fin:
        for line in fin:
            line = line.strip('\n')
            if line.find('0.0') == 0 or line.find('0.0') == 1:
                break
            shear_lst.append(float(line))
        #header = fin.readline()
        for line in fin:
            line = line.strip('\n')
            if line.find('0.0') == 0 or line.find('0.0') == 1:
                break
            torque_lst.append(float(line))
        #header = fin.readline()
        for line in fin:
            line = line.strip('\n')
            if line.find('0.0') == 0 or line.find('0.0') == 1:
                break
            moment_lst.append(float(line))

elif sw == 1:
    with open('positive_OEW+Payload+fuel.txt', 'r') as fin:
        for line in fin:
            line = line.strip('\n')
            if line.find('0.0') == 0 or line.find('0.0') == 1:
                break
            shear_lst.append(float(line))
        #header = fin.readline()
        for line in fin:
            line = line.strip('\n')
            if line.find('0.0') == 0 or line.find('0.0') == 1:
                break
            torque_lst.append(float(line))
        #header = fin.readline()
        for line in fin:
            line = line.strip('\n')
            if line.find('0.0') == 0 or line.find('0.0') == 1:
                break
            moment_lst.append(float(line))

elif sw == 2:
    with open('negative_OEW+Payload+fuel.txt', 'r') as fin:
        for line in fin:
            line = line.strip('\n')
            if line.find('0.0') == 0 or line.find('0.0') == 1:
                break
            shear_lst.append(float(line))
        #header = fin.readline()
        for line in fin:
            line = line.strip('\n')
            if line.find('0.0') == 0 or line.find('0.0') == 1:
                break
            torque_lst.append(float(line))
        #header = fin.readline()
        for line in fin:
            line = line.strip('\n')
            if line.find('0.0') == 0 or line.find('0.0') == 1:
                break
            moment_lst.append(float(line))

shear_lst.append(0)
torque_lst.append(0)
moment_lst.append(0)

Y = np.linspace(0, 33.226, 50)

shear_function = sp.interpolate.interp1d(Y,shear_lst,kind='cubic',fill_value="extrapolate")
torque_function = sp.interpolate.interp1d(Y,torque_lst,kind='cubic',fill_value="extrapolate")
moment_function = sp.interpolate.interp1d(Y,moment_lst,kind='cubic',fill_value="extrapolate")

#twist angle
twist_angle = integrate.quad(lambda y:  torque_function(y) / (J(y) * G), 0, b/2 )
print("The twist angle is: ", twist_angle[0]*180/np.pi)
print("the estimated error is: ", twist_angle[1]*180/np.pi)
print("The right spar thickness is: ", right_spart(0))

#twist angle vs spanwise location
y_values = np.linspace(0, b/2, 100)
twist_values = [integrate.quad(lambda y: torque_function(y) *180/ (J(y) * G* np.pi), 0, y)[0] for y in y_values]

J_values = [J(y) for y in y_values]

#J Value plot
plt.plot(y_values, J_values, label='J(y)')
plt.xlabel('y')
plt.ylabel('J(y)')
plt.title('J(y) vs. y')
plt.legend()
plt.grid(True)
plt.show()

#twist angle plot
   
plt.plot(y_values, twist_values, label='Twist (y)')
plt.xlabel('y')
plt.ylabel('angle')
plt.title('twist angle vs. y')
plt.legend()
plt.grid(True)
plt.show()




