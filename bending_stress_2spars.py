import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad
import scipy as sp
from scipy import interpolate

b = 66.45150601     #span
c_r  = 10.229604    #root chord
taper = 0.28        #taper ratio

y = np.linspace(0, 33.225, 50)

def chord(y):
    chord = c_r - c_r * (1-taper) / (b/2) * y
    return chord 





# Function to calculate weighted centroid and moment of inertia
def calculate_centroid_and_moment_of_inertia(spar_areas, spar_coordinates, point_coordinates, point_areas, centroid):
    total_area = sum(spar_areas) + sum(point_areas)

    centroid_sum = np.zeros(2)
    moment_of_inertia_sum = 0

    # Lists to store individual moment of inertia values for spar and point areas
    spar_moment_of_inertia = []
    point_moment_of_inertia = []

    # Lists to store differences between y-coordinates and y-value of the centroid
    spar_y_diff = []
    point_y_diff = []

    # Calculate centroid and moment of inertia contribution from spar areas
    for i, area in enumerate(spar_areas):
        centroid_sum += area * np.array(spar_coordinates[i])
        spar_y_diff_i = spar_coordinates[i][1] - centroid[1]
        spar_moment_of_inertia_i = area * spar_y_diff_i ** 2
        moment_of_inertia_sum += spar_moment_of_inertia_i
        spar_moment_of_inertia.append(spar_moment_of_inertia_i)
        spar_y_diff.append(spar_y_diff_i)

    # Calculate centroid and moment of inertia contribution from point areas
    for i in range(len(point_areas)):
        centroid_sum += point_areas[i] * np.array(point_coordinates[i])
        point_y_diff_i = point_coordinates[i][1] - centroid[1]
        point_moment_of_inertia_i = point_areas[i] * point_y_diff_i ** 2
        moment_of_inertia_sum += point_moment_of_inertia_i
        point_moment_of_inertia.append(point_moment_of_inertia_i)
        point_y_diff.append(point_y_diff_i)

    # Calculate the weighted centroid
    weighted_centroid = centroid_sum / total_area

    return weighted_centroid, moment_of_inertia_sum, spar_moment_of_inertia, point_moment_of_inertia, spar_y_diff, point_y_diff

# Function to convert coordinates to the new axis system with the weighted centroid as the origin
def convert_to_new_axis_system(coordinates, weighted_centroid):
    return [(coord[0] - weighted_centroid[0], coord[1] - weighted_centroid[1]) for coord in coordinates]

# Coordinates of spar points
spar_coordinates = [(0, 0), (0.40, 0.01238), (0.40, 0.12608), (0, 0.1266)]

# Prompt user for the area of spar point areas
area_spar = float(input("Enter the area of spar point areas: "))
spar_areas = [area_spar] * len(spar_coordinates)

# Number of point areas between spar 1 and 2
num_point_areas_1_2 = int(input("Enter the number of point areas between point 1 and 2: "))
# Number of point areas between spar 3 and 4
num_point_areas_3_4 = int(input("Enter the number of point areas between point 3 and 4: "))

# Generate evenly spaced point areas between coordinates of spar point areas at point 1 and 2
point_coordinates_1_2 = [(spar_coordinates[0][0] + i * (spar_coordinates[1][0] - spar_coordinates[0][0]) / (num_point_areas_1_2 + 1),
                          spar_coordinates[0][1] + i * (spar_coordinates[1][1] - spar_coordinates[0][1]) / (num_point_areas_1_2 + 1))
                         for i in range(1, num_point_areas_1_2 + 1)]

# Generate evenly spaced point areas between coordinates of spar point areas at point 3 and 4
point_coordinates_3_4 = [(spar_coordinates[2][0] + i * (spar_coordinates[3][0] - spar_coordinates[2][0]) / (num_point_areas_3_4 + 1),
                          spar_coordinates[2][1] + i * (spar_coordinates[3][1] - spar_coordinates[2][1]) / (num_point_areas_3_4 + 1))
                         for i in range(1, num_point_areas_3_4 + 1)]

# Prompt user for the area of point areas
area_point = float(input("Enter the area of point areas: "))
point_areas_3_4 = [area_point] * len(point_coordinates_3_4)
point_areas_1_2 = [area_point] * len(point_coordinates_1_2)

# Calculate the centroid and moment of inertia
weighted_centroid, moment_of_inertia, spar_moment_of_inertia, point_moment_of_inertia, spar_y_diff, point_y_diff = calculate_centroid_and_moment_of_inertia(
    spar_areas, spar_coordinates, point_coordinates_3_4 + point_coordinates_1_2,
    point_areas_3_4 + point_areas_1_2, centroid=np.array([0, 0])
)

# Convert spar and point coordinates to the new axis system
spar_coordinates_new_axis = convert_to_new_axis_system(spar_coordinates, weighted_centroid)
point_coordinates_1_2_new_axis = convert_to_new_axis_system(point_coordinates_1_2, weighted_centroid)
point_coordinates_3_4_new_axis = convert_to_new_axis_system(point_coordinates_3_4, weighted_centroid)

# Calculate the moment of inertia in the new axis system
moment_of_inertia_new_axis = 0

for i, area in enumerate(spar_areas):
    y_value_new_axis = spar_coordinates_new_axis[i][1]
    spar_moment_of_inertia_i = area * y_value_new_axis ** 2
    moment_of_inertia_new_axis += spar_moment_of_inertia_i

for i in range(len(point_areas_1_2)):
    y_value_new_axis = point_coordinates_1_2_new_axis[i][1]
    point_moment_of_inertia_i = point_areas_1_2[i] * y_value_new_axis ** 2
    moment_of_inertia_new_axis += point_moment_of_inertia_i

for i in range(len(point_areas_3_4)):
    y_value_new_axis = point_coordinates_3_4_new_axis[i][1]
    point_moment_of_inertia_i = point_areas_3_4[i] * y_value_new_axis ** 2
    moment_of_inertia_new_axis += point_moment_of_inertia_i


chord_values = chord(y)
multiplied_moment_of_inertia = chord(y)**3 * moment_of_inertia_new_axis



# Multiply the total moment of inertia by the array [10.22960, -0.2216748821]
#multiplied_moment_of_inertia = np.square(chord(np.array([0, 0]))) * moment_of_inertia_new_axis

#print(f"Multiplied Moment of Inertia Array: {multiplied_moment_of_inertia}")

# Interpolate to get a function
moi = interp1d(y, multiplied_moment_of_inertia, kind='cubic', fill_value="extrapolate")

# Plotting
plt.figure(figsize=(10, 6))

# Plot spar points in the new axis system
plt.scatter(*zip(*spar_coordinates_new_axis), s=np.array(spar_areas) * 500, label='Spar Areas (New Axis System)', color='red', alpha=0.7)

# Plot point areas between point 3 and 4 in the new axis system
plt.scatter(*zip(*point_coordinates_3_4_new_axis), s=np.array(point_areas_3_4) * 500, label='Point Areas (3 to 4, New Axis System)', color='blue', alpha=0.7)

# Plot point areas between point 1 and 2 in the new axis system
plt.scatter(*zip(*point_coordinates_1_2_new_axis), s=np.array(point_areas_1_2) * 500, label='Point Areas (1 to 2, New Axis System)', color='green', alpha=0.7)

# Plot the new y-axis (centroid as the origin) in the new axis system
plt.axvline(0, color='purple', linestyle='--', label='New Y-axis (centroid as origin, New Axis System)')

# Plot the new x-axis (centroid as the origin) in the new axis system
plt.axhline(0, color='orange', linestyle='--', label='New X-axis (centroid as origin, New Axis System)')

# Add labels and legend
plt.xlabel('X-axis (New Axis System)')
plt.ylabel('Y-axis (New Axis System)')
plt.axhline(0, color='black', linewidth=0.5)  # Zero line
plt.legend()
plt.title('Weighted Centroid Calculation with New Axes (Centroid as Origin)')

# Show the plot
plt.show()

# Generate y values for smooth plotting
moi_smooth = moi(y)

# Plotting
plt.figure(figsize=(10, 6))

# Plot the function moi
plt.plot(y, moi_smooth, label='Multiplied Moment of Inertia (Function)', color='purple')

# Add labels and legend
plt.xlabel('y')
plt.ylabel('Multiplied Moment of Inertia')
plt.axhline(0, color='black', linewidth=0.5)  # Zero line
plt.legend()
plt.title('Multiplied Moment of Inertia')
plt.show()


# Given data
e = 68.9 * 10**9

m = np.array([[25029385.656322654, 23987436.860227615, 22972151.460157588, 21984303.595079064,
               21024422.74875815, 20093439.382383946, 19190356.42432302, 18320966.458186615,
               17452779.77179053, 16586035.222988505, 15753703.804096347, 14950502.511911029,
               14177940.417561684, 13435756.449944522, 12723723.66405924, 12043015.999765215,
               11388956.881350175, 10778673.798561746, 10107751.31967802, 9376032.40002238,
               8687629.705944244, 8025110.133390517, 7392862.166318535, 6789398.610918606,
               6214779.229427664, 5668622.664489855, 5150640.046910202, 4660492.647419574,
               4197830.6947142305, 3762283.4882114544, 3353462.533430459, 2970961.0235873763,
               2614354.123328785, 2283198.082389323, 1977034.6066072197, 1695381.4925145633,
               1437741.235206663, 1203596.1377194359, 992408.603715985, 803620.3326889806,
               636649.5767562416, 490887.3529702796, 365695.2372573205, 260405.48664167273,
               174312.91168651686, 106641.5328194136, 56544.60472525939, 23083.79864952872,
               5067.459384203866, 0]])

# Calculate the y-coordinate of the centroid
centroid_y = weighted_centroid[1]

# Calculate the y-coordinate of each point mass
point_mass_y_coordinates = [coord[1] for coord in point_coordinates_1_2 + point_coordinates_3_4]

# Find the absolute difference between the y-coordinate of the centroid and each point mass
distances_from_centroid = [abs(centroid_y - y) for y in point_mass_y_coordinates]

# Find the maximum distance
max_distance = max(distances_from_centroid)

#print("Maximum distance on the y-axis between centroid and any point mass:", max_distance)

shear_lst = []
torque_lst = []
moment_lst = []

"""
sw = 0 for positive load factor, maximum shear case
sw = 1 for positive load factor, maximum bending stress case
sw = 2 for negative load factor, maximum shear and maximum bending stress case
"""
sw = 1


if sw == 0:
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



shear_function = sp.interpolate.interp1d(y,shear_lst,kind='cubic',fill_value="extrapolate")
torque_function = sp.interpolate.interp1d(y,torque_lst,kind='cubic',fill_value="extrapolate")
moment_function = sp.interpolate.interp1d(y,moment_lst,kind='cubic',fill_value="extrapolate")


# Interpolate the given data
#g = interp1d(y, moment_lst, kind='cubic', fill_value="extrapolate")

# Interpolate to get a function
moi = interp1d(y, multiplied_moment_of_inertia, kind='cubic', fill_value="extrapolate")

# Define the moment of inertia function
#def moi_function(y):
#    return 2.65301785e+00 -1.14981463e-01*y + 1.24582058e-03*y**2
#    return 5.68905490e-01 -2.46562931e-02*y + 2.67150169e-04*y**2
#    return 3.79270326e-01 -1.64375287e-02*y + 1.78100112e-04*y**2
#    return 1.89635163e+00 -8.21876437e-02*y + 8.90500562e-04*y**2
#    return 3.79270326e+00 -1.64375287e-01*y + 1.78100112e-03*y**2 #this satisfies requirements
#     return 9.48175816e+00 -4.10938219e-01*y + 4.45250281e-03*y**2
#    return 1.89635163e+01 -8.21876437e-01*y + 8.90500562e-03*y**2
#1.89635163e2 -8.21876437e1*y + 8.90500562e0*y**2

# Define your functions here (moment_function, chord, moi)

# Define the integrand function
def bending_stress_function(y):
    return moment_function(y) * max_distance * chord(y) / moi(y)

# Allowable stress (you need to set this based on your design requirements)
allowable_stress = 276000000  # Replace with the actual allowable stress value

# Calculate the margin of safety
margin_of_safety = allowable_stress / (bending_stress_function(y) * 1.5)  # With Safety factor

# Find the minimum margin of safety (most critical location)
min_margin_of_safety = np.min(margin_of_safety)

# Plot the margin of safety
plt.figure(figsize=(10, 6))

# Plot margin of safety
plt.plot(y, margin_of_safety, label='Margin of Safety', color='red')

# Highlight the most critical location
plt.scatter(y[np.argmin(margin_of_safety)], min_margin_of_safety, color='red', marker='o', label='Critical Point')

# Add labels and legend
plt.xlabel('Spanwise Location (m)')
plt.ylabel('Margin of Safety')
plt.xlim(-1, 35)
plt.ylim(0, 5)
plt.axhline(0, color='black', linewidth=0.5)  # Zero line
plt.legend()
plt.title('Margin of Safety along the Wing')

# Show the plot
plt.show()

# Print the minimum margin of safety and its corresponding location
print(f"Minimum Margin of Safety: {min_margin_of_safety}")
print(f"Corresponding Location (y): {y[np.argmin(margin_of_safety)]}")

# Initialize maxim before the loop
maxim = 0

# Second plot for lowest margin of safety every four meters
plt.figure(figsize=(10, 6))

# Calculate the lowest margin of safety for every four meters
y_intervals = np.arange(0, 35, 4)
min_margin_of_safety_intervals = []

# Iterate over intervals
for interval in y_intervals:
    # Get the margin of safety within the interval
    margin_in_interval = margin_of_safety[(y >= interval) & (y < interval + 4)]

    # Find the maximum value within the interval
    max_margin_in_interval = np.max(margin_in_interval)

    # Update maxim if needed
    if max_margin_in_interval > maxim:
        maxim = max_margin_in_interval

    # Append the minimum margin of safety within the interval
    min_margin_of_safety_intervals.append(np.min(margin_in_interval))

# Plot the lowest margin of safety for every four meters
plt.plot(y_intervals, min_margin_of_safety_intervals, marker='o', linestyle='-', color='blue', label='Lowest Margin of Safety')

# Add labels and legend
plt.xlabel('Spanwise Location (m)')
plt.ylabel('Lowest Margin of Safety')
plt.xlim(-1, 35)
plt.ylim(0, 5)
plt.axhline(0, color='black', linewidth=0.5)  # Zero line
plt.legend()
plt.title('Lowest Margin of Safety every Four Meters')

# Show the plot
plt.show()

# Print the minimum margin of safety and its corresponding location
print(f"Minimum Margin of Safety: {min_margin_of_safety}")
print(f"Corresponding Location (y): {y[np.argmin(margin_of_safety)]}")




for i in range(0, 50):
    print (bending_stress_function(y[i])/10**6)
    if bending_stress_function(y[i])/10**6>maxim:
        maxim = bending_stress_function(i)/10**6

print (maxim)
