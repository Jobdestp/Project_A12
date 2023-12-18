# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
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

def stringer_length(t):
    return (area_of_cross_section + t ** 2) / (2 * t)

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



