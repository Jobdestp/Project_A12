# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:36:40 2023

@author: CKW
"""
import numpy as np

#skin buckling
import math

b_w = 66.45150601     #span
c_r  = 10.229604    #root chord
taper = 0.28        #taper ratio
t = float(input('skin thickness:'))
E = 68.9 * 10**9
v = 0.33
k_c = 4
lst_skin_critical_stress = []

y = np.array([0.0, 0.6768367346938775, 1.353673469387755, 2.0305102040816325, 2.70734693877551,
              3.384183673469388, 4.061020408163266, 4.737857142857144, 5.414693877551021, 6.091530612244898,
              6.768367346938776, 7.445204081632653, 8.12204081632653, 8.798877551020408, 9.475714285714286,
              10.152551020408163, 10.829387755102041, 11.506224489795918, 12.183061224489796, 12.859897959183673,
              13.53673469387755, 14.213571428571428, 14.890408163265305, 15.567244897959183, 16.24408163265306,
              16.920918367346937, 17.597755102040817, 18.274591836734694, 18.95142857142857, 19.628265306122447,
              20.305102040816325, 20.981938775510204, 21.65877551020408, 22.335612244897958, 23.012448979591837,
              23.689285714285715, 24.366122448979594, 25.04295918367347, 25.71979591836735, 26.396632653061224,
              27.073469387755105, 27.75030612244898, 28.42714285714286, 29.10397959183674, 29.780816326530615,
              30.457653061224494, 31.13448979591837, 31.81132653061225, 32.48816326530613, 33.225])


def chord(y):
    chord = c_r - c_r * (1-taper) / (b_w/2) * y
    return chord 

numberstringers = (float(input("# of stringers")))

for i in range(0, len(y)):
    # Length of a and b at different sections
    b_stringerspacing = 0.4 * (c_r - c_r * ((1 - taper) / (b_w / 2)) * y[i]) / (numberstringers + 1)



#lst_aspect_ratio=[4]
#lst_web_buckling_coefficient=[]
#lst_skin_buckling_coefficient=[]
#lst_web_critical_stress=[]
#lst_skin_critical_stress=[]

    


    skin_critical_stress = ((math.pi ** 2 * k_c * E) / (12 * (1 - v ** 2))) * ( t / b_stringerspacing ) ** 2
    lst_skin_critical_stress.append(skin_critical_stress)

print(lst_skin_critical_stress)



