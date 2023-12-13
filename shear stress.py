import numpy as np
import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt



"""
sw = 0 for positive load factor, maximum shear case
sw = 1 for positive load factor, maximum bending stress case
sw = 2 for negative load factor, maximum shear and maximum bending stress case
"""

b = 66.45
c_r = 10.22960376
taper = 0.28
G = 26*10**9

f = float(input("where is the middle spar? "))
t_front = float(input("what is the front spar web thickness? "))
t_middle = float(input("What is the middle spar web thickness? "))
t_rear = float(input("What is the rear spar web thickness? "))
t_skin = float(input("What is the skin thickness? "))


#chord
def chord(y):
    chord = c_r - c_r * (1-taper) / (b/2) * y
    return chord

#spar height
def h_front(y):
    front_sparheight = 0.1266068 * chord(y)
    return front_sparheight

def h_rear(y):
    rear_sparheight = 0.11369081 * chord(y)
    return rear_sparheight

def h_middle(y):
    middle_sparheight = h_front(y) * (1-f) + h_rear(y) * f
    return middle_sparheight



for j in range(3):
    shear_lst = []
    torque_lst = []
    moment_lst = []

    sw = j


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


    y_values = np.linspace(0, b/2, 50)

    

   

    #Average Shear

    def V_shear(y):
        average_shear = shear_function(y) / (h_front(y) * t_front + h_middle(y) * t_middle + h_rear(y) * t_rear)
        return average_shear

    #max shear is average shear multiplied by 1.5 + safety factor of 1.5

    V_shearvalues = []
    

    for y in y_values:
        V_shearvalues.append(V_shear(y) * 2.25/10**6)


    #torsion shear
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

    #create lists of total shear stress in Mpa with safety factor 1.5
    for i in range(50):
        q1, q2 = solver(y_values[i])
        front_stress[i] = 1.5*q1/(t_front *10**6) + V_shearvalues[i]
        rear_stress[i] = -1.5*q2/(t_rear*10**6) + V_shearvalues[i]
        middle_stress[i] = -1.5*(q1-q2)/(t_middle*10**6) + V_shearvalues[i]
        skin_stress_q1[i] = 1.5*q1/(t_skin*10**6)
        skin_stress_q2[i] = 1.5*q2/(t_skin*10**6)

    

    #interpolate the shear stress list to create a function
    front_stressfunction = sp.interpolate.interp1d(y_values,front_stress,kind='cubic',fill_value="extrapolate")
    middle_stressfunction = sp.interpolate.interp1d(y_values,middle_stress,kind='cubic',fill_value="extrapolate")
    rear_stressfunction = sp.interpolate.interp1d(y_values,rear_stress,kind='cubic',fill_value="extrapolate")

    
    #print information
    if sw == 0:
        print("The maximum shear stress due to sw", j, "is" , max(front_stress, key = abs), max(middle_stress, key = abs), max(rear_stress, key = abs), "MPa")

        plt.plot(y_values, front_stress, label='Shear Stress')
        plt.xlabel('spanwise location [m]')
        plt.ylabel('Shear Stress [MPa]')
        plt.title('Shear Stress vs y Values')
        plt.legend()
        plt.grid(True)
        plt.show()


    elif sw == 2:
        print("The minimum shear stress due to sw", j, "is" , max(front_stress, key = abs), max(middle_stress, key = abs), max(rear_stress, key = abs), "MPa")
        plt.plot(y_values, front_stress, label='Shear Stress')
        plt.xlabel('spanwise location [m]')
        plt.ylabel('Shear Stress [MPa]')
        plt.title('Shear Stress vs y Values')
        plt.legend()
        plt.grid(True)
        plt.show()

    else:
        print("The maximum shear stress due to shear in sw", j, "is" , max(front_stress, key = abs), max(middle_stress, key = abs), max(rear_stress, key = abs), "MPa")
        plt.plot(y_values, front_stress, label='Shear Stress')
        plt.xlabel('spanwise location [m]')
        plt.ylabel('Shear Stress [MPa]')
        plt.title('Shear Stress vs y Values')
        plt.legend()
        plt.grid(True)
        plt.show()


