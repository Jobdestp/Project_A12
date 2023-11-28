import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad


b = 66.45150601     # span
c_r = 10.229604      # root chord
taper = 0.28         # taper ratio

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
    return c_r - c_r * (1 - taper) / (b / 2) * y

def calculate_centroid_and_moment_of_inertia(spar_areas, spar_coordinates, point_coordinates_top, point_coordinates_bot, point_areas_top, point_areas_bot):
    total_area = sum(spar_areas) + sum(point_areas_top) + sum(point_areas_bot)

    centroid_sum = np.zeros(2)
    moment_of_inertia_sum = 0

    spar_moment_of_inertia = []
    point_moment_of_inertia = []
    spar_y_diff = []
    point_y_diff = []

    for i, area in enumerate(spar_areas):
        centroid_sum += area * np.array(spar_coordinates[i])
        spar_y_diff_i = spar_coordinates[i][1] - centroid_sum[1]
        spar_moment_of_inertia_i = area * spar_y_diff_i ** 2
        moment_of_inertia_sum += spar_moment_of_inertia_i
        spar_moment_of_inertia.append(spar_moment_of_inertia_i)
        spar_y_diff.append(spar_y_diff_i)

    for i in range(len(point_areas_top)):
        centroid_sum += point_areas_top[i] * np.array(point_coordinates_top[i])
        point_y_diff_i = point_coordinates_top[i][1] - centroid_sum[1]
        point_moment_of_inertia_i = point_areas_top[i] * point_y_diff_i ** 2
        moment_of_inertia_sum += point_moment_of_inertia_i
        point_moment_of_inertia.append(point_moment_of_inertia_i)
        point_y_diff.append(point_y_diff_i)

    for i in range(len(point_areas_bot)):
        centroid_sum += point_areas_bot[i] * np.array(point_coordinates_bot[i])
        point_y_diff_i = point_coordinates_bot[i][1] - centroid_sum[1]
        point_moment_of_inertia_i = point_areas_bot[i] * point_y_diff_i ** 2
        moment_of_inertia_sum += point_moment_of_inertia_i
        point_moment_of_inertia.append(point_moment_of_inertia_i)
        point_y_diff.append(point_y_diff_i)

    weighted_centroid = centroid_sum / total_area

    return weighted_centroid, moment_of_inertia_sum, spar_moment_of_inertia, point_moment_of_inertia, spar_y_diff, point_y_diff

def convert_to_new_axis_system(coordinates, weighted_centroid):
    return [(coord[0] - weighted_centroid[0], coord[1] - weighted_centroid[1]) for coord in coordinates]

spar_coordinates = [(0, 0), (0.20, 0.00619), (0.40, 0.01238), (0.40, 0.12608), (0.20, 0.12634), (0, 0.1266)]

area_spar = float(input("Enter the area of spar point areas: "))
spar_areas = [area_spar] * len(spar_coordinates)

num_point_areas_bot_01 = int(input("Enter the number of point areas bottom between spar 0 and 1: "))
num_point_areas_bot_12 = int(input("Enter the number of point areas bottom between spar 1 and 2: "))
num_point_areas_top_34 = int(input("Enter the number of point areas top between spar 3 and 4: "))
num_point_areas_top_45 = int(input("Enter the number of point areas top between spar 4 and 5: "))

point_coordinates_bot_01 = [
    (spar_coordinates[0][0] + i * (spar_coordinates[1][0] - spar_coordinates[0][0]) / (num_point_areas_bot_01 + 1),
     spar_coordinates[0][1] + i * (spar_coordinates[1][1] - spar_coordinates[0][1]) / (num_point_areas_bot_01 + 1))
    for i in range(1, num_point_areas_bot_01 + 1)]

point_coordinates_bot_12 = [
    (spar_coordinates[1][0] + i * (spar_coordinates[2][0] - spar_coordinates[1][0]) / (num_point_areas_bot_12 + 1),
     spar_coordinates[1][1] + i * (spar_coordinates[2][1] - spar_coordinates[1][1]) / (num_point_areas_bot_12 + 1))
    for i in range(1, num_point_areas_bot_12 + 1)]

point_coordinates_top_34 = [
    (spar_coordinates[3][0] + i * (spar_coordinates[4][0] - spar_coordinates[3][0]) / (num_point_areas_top_34 + 1),
     spar_coordinates[3][1] + i * (spar_coordinates[4][1] - spar_coordinates[3][1]) / (num_point_areas_top_34 + 1))
    for i in range(1, num_point_areas_top_34 + 1)]

point_coordinates_top_45 = [
    (spar_coordinates[4][0] + i * (spar_coordinates[5][0] - spar_coordinates[4][0]) / (num_point_areas_top_45 + 1),
     spar_coordinates[4][1] + i * (spar_coordinates[5][1] - spar_coordinates[4][1]) / (num_point_areas_top_45 + 1))
    for i in range(1, num_point_areas_top_45 + 1)]

point_coordinates_bot = point_coordinates_bot_01 + point_coordinates_bot_12
point_coordinates_top = point_coordinates_top_34 + point_coordinates_top_45

area_point = float(input("Enter the area of point areas: "))
point_areas_top = [area_point] * len(point_coordinates_top)
point_areas_bot = [area_point] * len(point_coordinates_bot)

weighted_centroid, moment_of_inertia, spar_moment_of_inertia, point_moment_of_inertia, spar_y_diff, point_y_diff = calculate_centroid_and_moment_of_inertia(
    spar_areas, spar_coordinates,
    point_coordinates_top, point_coordinates_bot, point_areas_top, point_areas_bot)

spar_coordinates_new_axis = convert_to_new_axis_system(spar_coordinates, weighted_centroid)
point_coordinates_bot_new_axis = convert_to_new_axis_system(point_coordinates_bot, weighted_centroid)
point_coordinates_top_new_axis = convert_to_new_axis_system(point_coordinates_top, weighted_centroid)

moment_of_inertia_new_axis = 0

for i, area in enumerate(spar_areas):
    y_value_new_axis = spar_coordinates_new_axis[i][1]
    spar_moment_of_inertia_i = area * y_value_new_axis ** 2
    moment_of_inertia_new_axis += spar_moment_of_inertia_i

for i in range(len(point_areas_bot)):
    y_value_new_axis = point_coordinates_bot_new_axis[i][1]
    point_moment_of_inertia_i = point_areas_bot[i] * y_value_new_axis ** 2
    moment_of_inertia_new_axis += point_moment_of_inertia_i

for i in range(len(point_areas_top)):
    y_value_new_axis = point_coordinates_top_new_axis[i][1]
    point_moment_of_inertia_i = point_areas_top[i] * y_value_new_axis ** 2
    moment_of_inertia_new_axis += point_moment_of_inertia_i

chord_values = chord(y)
multiplied_moment_of_inertia = chord(y)**3 * moment_of_inertia_new_axis

print(f"Multiplied Moment of Inertia Array: {multiplied_moment_of_inertia}")

moi = interp1d(y, multiplied_moment_of_inertia, kind='cubic', fill_value="extrapolate")

plt.figure(figsize=(10, 6))

plt.scatter(*zip(*spar_coordinates_new_axis), s=np.array(spar_areas) * 500, label='Spar Areas (New Axis System)', color='red', alpha=0.7)
plt.scatter(*zip(*point_coordinates_top_new_axis), s=np.array(point_areas_top) * 500, label='Point Areas (Top, New Axis System)', color='blue', alpha=0.7)
plt.scatter(*zip(*point_coordinates_bot_new_axis), s=np.array(point_areas_bot) * 500, label='Point Areas (Bottom, New Axis System)', color='green', alpha=0.7)
plt.axvline(0, color='purple', linestyle='--', label='New Y-axis (centroid as origin, New Axis System)')
plt.axhline(0, color='orange', linestyle='--', label='New X-axis (centroid as origin, New Axis System)')

plt.xlabel('X-axis (New Axis System)')
plt.ylabel('Y-axis (New Axis System)')
plt.axhline(0, color='black', linewidth=0.5)
plt.legend()
plt.title('Weighted Centroid Calculation with New Axes (Centroid as Origin)')
plt.show()

moi_smooth = moi(y)

plt.figure(figsize=(10, 6))

plt.plot(y, moi_smooth, label='Multiplied Moment of Inertia (Function)', color='purple')
plt.xlabel('y')
plt.ylabel('Multiplied Moment of Inertia')
plt.axhline(0, color='black', linewidth=0.5)
plt.legend()
plt.title('Multiplied Moment of Inertia')
plt.show()

# Given data
e = 68.9 * 10**9
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

# Interpolate the given data
g = interp1d(y, m, kind='cubic', fill_value="extrapolate")

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

# Define the integrand function
def integrand_function(y):
    return g(y) / (e * moi(y) )

# Array to store integration results for n
n = []

# Loop through each y value and perform integration for n
for y_value in y:
    result, error = quad(integrand_function, 0, y_value)  # Increase the limit if needed
    n.append(result)

# Interpolate the integration results for n
h = interp1d(y, n, kind='cubic', fill_value="extrapolate")

# Array to store integration results for p
p = []

# Loop through each y value and perform integration for p
for y_value in y:
    result2, error2 = quad(h, 0, y_value)  # Increase the limit if needed
    p.append(result2)

# Plot the functions n(y) and p(y)
plt.plot(y, n, label='dv/dy(y)')
plt.plot(y, p, label='deflection(y)')
plt.title('Plot of dv/dy(y) and deflection(y)')
plt.xlabel('y')
plt.ylabel('Function Values')
plt.legend()
plt.grid(True)
plt.show()
