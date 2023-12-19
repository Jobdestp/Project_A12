# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import numpy as np
import  sympy as sp
import math

thickness = float(input('Enter the thickness of the stringer: '))
ribs_spacing = float(input('Enter the ribs spacing in meters: '))
area_of_cross_section = float(input("Enter the area of the L-shaped cross-section: "))
wingspan = 33.225

K = 1  # loading factor when the stringer is pinned at both ends
t = thickness
l = ribs_spacing
A = area_of_cross_section
K_skin = 4
v = 0.33
E = 68.9 * 10**9
sigma_tensile = 310 * 10 ** 6
c_r  = 10.229604
taper = 0.28
b_w = 66.45150601
Critical_Buckling_Stress_New = []



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

for i in range(0, len(y)):
    area_crossection_new = area_of_cross_section * (c_r - c_r * (1-taper) / (b_w/2) * y[i])

    t_new = t * (c_r - c_r * (1-taper) / (b_w/2) * y[i])

    stringer_length = (area_crossection_new + t_new ** 2) / (2 * t_new)

    I_new = (1/(12) * stringer_length * t_new**3 + area_crossection_new * (t_new/2)**2 + 1/12 * t_new * (stringer_length - t_new)**3 + area_crossection_new * ((stringer_length - t_new)/2 + t_new)**2)

    sigma_column = (K * math.pi ** 2 * E * I_new) / (l ** 2 * area_crossection_new)

    Critical_Buckling_Stress_New.append(sigma_column)


print('The critical column buckling stress is:', Critical_Buckling_Stress_New)



