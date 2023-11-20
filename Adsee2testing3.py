import numpy as np
import matplotlib.pyplot as plt

# Function to calculate weighted centroid and moment of inertia
def calculate_centroid_and_moment_of_inertia(spar_areas, spar_coordinates, point_coordinates, point_areas, centroid):
    total_area = sum(spar_areas) + sum(point_areas)

    centroid_sum = np.zeros(2)
    moment_of_inertia_sum = 0

    # Calculate centroid and moment of inertia contribution from spar areas
    for i, area in enumerate(spar_areas):
        centroid_sum += area * np.array(spar_coordinates[i])
        moment_of_inertia_sum += area * np.square(spar_coordinates[i][1] - centroid[1])

    # Calculate centroid and moment of inertia contribution from point areas
    for i in range(len(point_areas)):
        centroid_sum += point_areas[i] * np.array(point_coordinates[i])
        moment_of_inertia_sum += point_areas[i] * np.square(point_coordinates[i][1] - centroid[1])

    # Calculate the weighted centroid
    centroid = centroid_sum / total_area

    return centroid, moment_of_inertia_sum

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
centroid, moment_of_inertia = calculate_centroid_and_moment_of_inertia(
    spar_areas, spar_coordinates, point_coordinates_3_4 + point_coordinates_1_2,
    point_areas_3_4 + point_areas_1_2, centroid=np.array([0, 0])
)

# Print the coordinates of spar point areas
print("Spar Point Areas:")
for i, coord in enumerate(spar_coordinates, start=1):
    print(f"{i}. {coord}")

# Print the coordinates of point areas between Spar 1 and Spar 2
print("\nPoint Areas (between Spar 1 and Spar 2):")
for i, coord in enumerate(point_coordinates_1_2, start=1):
    print(f"{i}. {coord}")

# Print the coordinates of point areas between Spar 3 and Spar 4
print("\nPoint Areas (between Spar 3 and Spar 4):")
for i, coord in enumerate(point_coordinates_3_4, start=1):
    print(f"{i}. {coord}")

# Print the centroid and moment of inertia
print("\nWeighted Centroid:", centroid)
print("Summation of Moment of Inertia about Y-axis:", moment_of_inertia)

# Plotting
plt.figure(figsize=(10, 6))

# Plot spar points
plt.scatter(*zip(*spar_coordinates), s=np.array(spar_areas) * 500, label='Spar Areas', color='red', alpha=0.7)

# Plot point areas between point 3 and 4
plt.scatter(*zip(*point_coordinates_3_4), s=np.array(point_areas_3_4) * 500, label='Point Areas (3 to 4)', color='blue', alpha=0.7)

# Plot point areas between point 1 and 2
plt.scatter(*zip(*point_coordinates_1_2), s=np.array(point_areas_1_2) * 500, label='Point Areas (1 to 2)', color='green', alpha=0.7)

# Plot the centroid
plt.scatter(*centroid, marker='x', color='black', label='Weighted Centroid', s=200)

# Add labels and legend
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.axhline(0, color='black', linewidth=0.5)  # Zero line
plt.legend()
plt.title('Weighted Centroid Calculation')

# Show the plot
plt.show()


