import numpy as np
import scipy as sp
from scipy import integrate
import matplotlib.pyplot as plt


halfspan=33.23
y_tab=np.arange(0, halfspan+0.01, 0.01) #y coordinate along the wing. y_max=halfspan        
v_tab=[]                                #Deflection list.
I_tab=[]                                #Area moment of inertia with changing y. List obtained from other Python codes.
E=68.9*10**9                                     #Young's modulus. Constant.

print(f'Max deflection should be:', 0.15*2*halfspan, f'm.')

#M_calculations
for y in y_tab:
    I_tab.append(0.002494)
    y=y+0.01



def I(y):
    return I_tab[int(y*100)]    #most probably will be replaced




for i in range(len(y_tab)):
    y_i=i/100
    ii=I(y_i)#will be replaced with I formula 
    f=lambda y,x: -(74.36*y**3+1617*y**2 - 3.161*10**5*y + 6.036*10**6)/E/ii #the stuff inside barckets is internal moment equation for cantaliever beam with distributed load of 10N/m 
    estimatedeflection, errordeflection=sp.integrate.dblquad(f, 0, y_i, lambda y:0, lambda y: y)
    v_tab.append(estimatedeflection)   

plt.plot(y_tab, v_tab)
plt.show()

print(f'Tip deflection is', abs(round(v_tab[-1], 3)) ,f'm.')
print(f'Corresponding I average value of:', I_tab[1], f'm^4.')
