import numpy as np
import scipy as sp
from scipy import interpolate

shear_lst = []
torque_lst = []
moment_lst = []

"""
sw = 0 for positive load factor, maximum shear case
sw = 1 for positive load factor, maximum bending stress case
sw = 2 for negative load factor, maximum shear and maximum bending stress case
"""
sw = 2


if sw == 2:
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
