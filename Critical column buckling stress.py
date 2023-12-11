# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import  sympy as sp
import math

length_1 = float(input('Enter the length of stringer segment: '))
area = float(input('Enter the area of a single stringer: '))
length_2 = float(input('Enter the length of the stringer: '))
thickness = float(input('Enter the thickness of the stringer: '))
area_of_cross_section = float(input("Enter the area of the L-shaped cross-section: "))
user_specified_thickness = float(input("Enter the thickness of each side (t): "))

K = 1  # loading factor when the stringer is pinned at both ends
A = area
h = length_2
l = length_1
t = thickness
E = 68.9 * 10**9



def calculate_dimensions(area, user_thickness):
    # Define symbols
    s, t = sp.symbols('s t', positive=True, real=True)

    # Define the equation for the area with user-specified thickness
    equation = sp.Eq((2 * s + user_thickness) * user_thickness, area)

    # Solve the equation for s
    solution = sp.solve(equation, s)

    # Extract the value for s from the first tuple in the solution list
    length = solution[0]

    return length, user_thickness

length, thickness = calculate_dimensions(area_of_cross_section, user_specified_thickness)

def calculate_second_moment_of_inertia(t, h, A):
    I = 1/(12) * h * t**3 + A * (t/2)**2 + 1/12 * t * (h - t)**3 + A * ((h - t)/2 + t)**2
    return I

result_I = calculate_second_moment_of_inertia(t, h, A)

print('The second moment of inertia of the stringer is:', result_I)

def calculate_column_buckling_stress(I, K, E, l, A):
    sigma_buckling = (K * math.pi ** 2 * E * I) / (l ** 2 * A)
    return sigma_buckling

result_sigma = calculate_column_buckling_stress(result_I, K, E, l, A)

print('The critical column buckling stress is :', result_sigma)
print(f"Length of each side: {length.evalf()}")
print(f"Thickness: {thickness}")


