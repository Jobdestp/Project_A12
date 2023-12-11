import sympy as sp

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

# Example usage:
area_of_cross_section = float(input("Enter the area of the L-shaped cross-section: "))
user_specified_thickness = float(input("Enter the thickness of each side (t): "))

length, thickness = calculate_dimensions(area_of_cross_section, user_specified_thickness)

print(f"Length of each side: {length.evalf()}")
print(f"Thickness: {thickness}")
