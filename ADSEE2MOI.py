import numpy as np
import matplotlib.pyplot as plt

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
spar_coordinates = [(0, 0), (0.4, 0.01238), (0.4, 0.12608), (0, 0.1266)]

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
spar_moment_of_inertia_new_axis = []
point_moment_of_inertia_new_axis = []

for i, area in enumerate(spar_areas):
    y_value_new_axis = spar_coordinates_new_axis[i][1]
    spar_moment_of_inertia_i = area * y_value_new_axis ** 2
    moment_of_inertia_new_axis += spar_moment_of_inertia_i
    spar_moment_of_inertia_new_axis.append(spar_moment_of_inertia_i)

for i in range(len(point_areas_1_2)):
    y_value_new_axis = point_coordinates_1_2_new_axis[i][1]
    point_moment_of_inertia_i = point_areas_1_2[i] * y_value_new_axis ** 2
    moment_of_inertia_new_axis += point_moment_of_inertia_i
    point_moment_of_inertia_new_axis.append(point_moment_of_inertia_i)

for i in range(len(point_areas_3_4)):
    y_value_new_axis = point_coordinates_3_4_new_axis[i][1]
    point_moment_of_inertia_i = point_areas_3_4[i] * y_value_new_axis ** 2
    moment_of_inertia_new_axis += point_moment_of_inertia_i
    point_moment_of_inertia_new_axis.append(point_moment_of_inertia_i)

# Print the coordinates of spar point areas in the new axis system
print("\nSpar Point Areas in the New Axis System:")
for i, coord in enumerate(spar_coordinates_new_axis, start=1):
    print(f"{i}. {coord}")

# Print the coordinates of point areas between Spar 1 and Spar 2 in the new axis system
print("\nPoint Areas (between Spar 1 and Spar 2) in the New Axis System:")
for i, coord in enumerate(point_coordinates_1_2_new_axis, start=1):
    print(f"{i}. {coord}")

# Print the coordinates of point areas between Spar 3 and Spar 4 in the new axis system
print("\nPoint Areas (between Spar 3 and Spar 4) in the New Axis System:")
for i, coord in enumerate(point_coordinates_3_4_new_axis, start=1):
    print(f"{i}. {coord}")

# Print the moment of inertia values in the new axis system
print("\nMoment of Inertia Values in the New Axis System:")
print(f"Total Moment of Inertia in the New Axis System: {moment_of_inertia_new_axis}")
print(f"Moment of Inertia for Spar Areas in the New Axis System: {spar_moment_of_inertia_new_axis}")
print(f"Moment of Inertia for Point Areas in the New Axis System: {point_moment_of_inertia_new_axis}")

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

