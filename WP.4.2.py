import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy import integrate
import sympy as symp

halfspan=33.23
y_tab=np.arange(0, halfspan+0.01, 0.01) #y coordinate along the wing. y_max=halfspan        
v_tab=[]                                #Deflection list.
I_tab=[]                                #Area moment of inertia with changing y. List obtained from other Python codes.
E=1                                     #Young's modulus. Constant.

#M_calculations
for y in y_tab:
    I_tab.append(1)
    y=y+0.01



def I(y):
    return I_tab[int(y*100)]    #most probably will be replaced




for i in range(len(y_tab)):
    y_i=i/100
    ii=I(y_i)#will be replaced with I formula 
    f=lambda y,x: -(5521.165-10*halfspan*y+10*y**2/2)/E/ii #the stuff inside barckets is internal moment equation for cantaliever beam with distributed load of 10N/m 
    estimatedeflection, errordeflection=sp.integrate.dblquad(f, 0, y_i, lambda y:0, lambda y: y)
    v_tab.append(estimatedeflection)   

plt.plot(y_tab, v_tab)
plt.show()

print(f'Tip deflection is', round(v_tab[-1], 2) ,f'm.')
