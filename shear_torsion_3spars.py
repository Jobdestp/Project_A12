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


#spar height
def h_front(y):
    left_sparheight = 0.1266068 * chord(y)
    return left_sparheight

def h_rear(y):
    right_sparheight = 0.11369081 * chord(y)
    return right_sparheight

#Get spar thickness
t_front = float(input("what is the front spar web thickness? "))
t_middle = float(input("what is the middle web thickness? "))
t_rear = float(input("what is the rear web thickness? "))


#Get top and bottom thickness
t_skin = float(input("what is the thickness of the top and bottom surfaces? "))

#Polar Moment of Inertia Function 





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

Y = np.linspace(0, 33.225, 50)

shear_function = sp.interpolate.interp1d(Y,shear_lst,kind='cubic',fill_value="extrapolate")
torque_function = sp.interpolate.interp1d(Y,torque_lst,kind='cubic',fill_value="extrapolate")
moment_function = sp.interpolate.interp1d(Y,moment_lst,kind='cubic',fill_value="extrapolate")

f = 0.5

def h_middle(y):
    return (1-f)*h_front(y)+f*h_rear(y)

def A_front(y):
    return f*0.2*chord(y)*(h_front(y)+h_middle(y))

def A_rear(y):
    return (1-f)*0.2*chord(y)*(h_rear(y)+h_middle(y))

def solver(y):
    c1 = 2*f*0.4*chord(y)/t_skin+h_front(y)/t_front+h_middle(y)/t_middle
    c2 = -h_middle(y)/t_middle
    c3 = -2*A_front(y)*G
    c4 = -h_middle(y)/t_middle
    c5 = 2*(1-f)*0.4*chord(y)/t_skin+h_rear(y)/t_rear+h_middle(y)/t_middle
    c6 = -2*A_rear(y)*G
    c7 = 2*A_front(y)
    c8 = 2*A_rear(y)
    c9 = 0
    matrix = np.array([[c1, c2, c3],
                   [c4, c5, c6],
                   [c7, c8, c9]])
    righthandside = np.array([0., 0., torque_function(y)])
    solution = np.linalg.solve(matrix, righthandside)
    return solution[0], solution[1]



front_stress = np.zeros(50)
rear_stress = np.zeros(50)
middle_stress = np.zeros(50)
skin_stress_q1 = np.zeros(50)
skin_stress_q2 = np.zeros(50)

for i in range(50):
    q1, q2 = solver(Y[i])
    print (q1, q2)
    front_stress[i] = q1/t_front
    rear_stress[i] = q2/t_rear
    middle_stress[i] = (q1-q2)/t_middle
    skin_stress_q1[i] = q1/t_skin
    skin_stress_q2[i] = q2/t_skin

print (min(front_stress/10**6))
if min(middle_stress) == 0:
    print (max(middle_stress/10**6))
else:
    print (min(middle_stress/10**6))
print (min(rear_stress/10**6))
print (min(skin_stress_q1/10**6))
print (min(skin_stress_q2/10**6))




    





