# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import  sympy as sp
import math

thickness = float(input('Enter the thickness of the stringer: '))
ribs_spacing = float(input('Enter the ribs spacing in meters: '))
area_of_cross_section = float(input("Enter the area of the L-shaped cross-section (from WP4): "))
wingspan = 33.225

K = 1  # loading factor when the stringer is pinned at both ends
t = thickness
l = ribs_spacing
A = area_of_cross_section
K_skin = 4
v = 0.33
E = 68.9 * 10**9
sigma_tensile = 310 * 10 ** 6

lst_column_critical_stress = []
lst_I = []

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

for i in range(0, len(y)):  
    def stringer_length(t):
        return (area_of_cross_section * (c_r - c_r * ((1 - taper) / (b_w / 2)) * y[i])) + (t * (c_r - c_r * ((1 - taper) / (b_w / 2)) * y[i]) ** 2) / (2 * t*(c_r - c_r * ((1 - taper) / (b_w / 2)) * y[i]))

    h = stringer_length(t)

    def rib_length(r):
        return wingspan / ribs_spacing

    result_rib = rib_length(ribs_spacing)

    def calculate_second_moment_of_inertia(t, h, A):
        return (1/(12) * h * t**3 + A * (t/2)**2 + 1/12 * t * (h - t)**3 + A * ((h - t)/2 + t)**2)

    result_I = calculate_second_moment_of_inertia(t, h, A)

    I = result_I

def calculate_column_buckling_stress(I, K, E, l, A):
    return (K * math.pi ** 2 * E * I) / (l ** 2 * A)


result_sigma = calculate_column_buckling_stress(I, K, E, l, A)

print('The second moment of inertia of the stringer is:', I)
print('There are ', result_rib, 'ribs on one side of the wing')
print('The critical column buckling stress is:', result_sigma)
print('The height of stringer is: ', h)




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
