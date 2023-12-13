# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 17:00:21 2023

@author: CKW
"""

import numpy as np
import math

b_w = 66.45150601  # span
c_r = 10.229604  # root chord
taper = 0.28  # taper ratio

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
    return c_r - c_r * (1 - taper) / (b_w / 2) * y

b_rib = float(input('Rib spacing: '))
t_spar = float(input('Thickness sparweb: '))

E = 68.9 * 10**9
v = 0.33
lst_web_buckling_coefficient_left = []
lst_web_buckling_coefficient_right = []
lst_web_buckling_coefficient_mid = []
lst_web_critical_stress_left = []
lst_web_critical_stress_right = []
lst_web_critical_stress_mid = []

for i in range(0, len(y)):
    # Length of a and b at different sections
    a_right = 0.11369 * (c_r - c_r * ((1 - taper) / (b_w / 2)) * y[i])
    a_left = 0.126607 * (c_r - c_r * ((1 - taper) / (b_w / 2)) * y[i])
    a_mid = a_right + a_left / 2
    b = b_rib

    # Calculate aspect ratio and buckling coefficient
    k_left = a_left / b
    k_right = a_right / b
    k_mid = a_mid / b
    k_s_left = 5.34 + 4 / (k_left ** 2)
    k_s_right = 5.34 + 4 / (k_right ** 2)
    k_s_mid = 5.34 + 4 / (k_mid ** 2)

    # Calculate critical web buckling stress and skin buckling stress
    web_critical_buckling_stress_left = ((math.pi ** 2 * k_s_left * E) / (12 * (1 - v ** 2))) * (t_spar / b) ** 2
    web_critical_buckling_stress_right = ((math.pi ** 2 * k_s_right * E) / (12 * (1 - v ** 2))) * (t_spar / b) ** 2
    web_critical_buckling_stress_mid = ((math.pi ** 2 * k_s_mid * E) / (12 * (1 - v ** 2))) * (t_spar / b) ** 2
    lst_web_critical_stress_left.append(web_critical_buckling_stress_left)
    lst_web_critical_stress_right.append(web_critical_buckling_stress_right)
    lst_web_critical_stress_mid.append(web_critical_buckling_stress_mid)
    lst_web_buckling_coefficient_left.append(k_s_left)
    lst_web_buckling_coefficient_right.append(k_s_right)
    lst_web_buckling_coefficient_mid.append(k_s_mid)

print('The web critical buckling stress is: ', lst_web_critical_stress_mid, lst_web_critical_stress_right,
      lst_web_critical_stress_left)
print('The web buckling coefficient is: ', lst_web_buckling_coefficient_mid, lst_web_buckling_coefficient_left,
      lst_web_buckling_coefficient_right)
