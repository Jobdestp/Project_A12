import matplotlib.pyplot as plt
import math

def trapezoid_centroid(x1, y1, x2, y2, x3, y3, x4, y4, thickness_left, thickness_right, thickness_top, thickness_bottom):
    # Calculate the areas and centroids of each section
    area_left = (y4 - y1) * thickness_left
    area_right = (y3 - y2) * thickness_right
    area_top = (x3 - x4) * thickness_top
    area_bottom = ((((x2 - x1)**2) + ((y2-y1)**2))**0.5) * thickness_bottom

    total_area = area_left + area_right + area_top + area_bottom

    # Calculate the centroids of each section
    x_centroid_left = (x1 + x4) / 2
    y_centroid_left = (y1 + y4) / 2

    x_centroid_right = (x2 + x3) / 2
    y_centroid_right = (y2 + y3) / 2

    x_centroid_top = (x3 + x4) / 2
    y_centroid_top = (y3 + y4) / 2

    x_centroid_bottom = (x1 + x2) / 2
    y_centroid_bottom = (y1 + y2) / 2

    # Calculate the weighted average of the centroids based on the areas
    x_centroid = (area_left * x_centroid_left + area_right * x_centroid_right + area_top * x_centroid_top + area_bottom * x_centroid_bottom) / total_area
    y_centroid = (area_left * y_centroid_left + area_right * y_centroid_right + area_top * y_centroid_top + area_bottom * y_centroid_bottom) / total_area

    return x_centroid, y_centroid

def plot_trapezoid(x1, y1, x2, y2, x3, y3, x4, y4, thickness_left, thickness_right, thickness_top, thickness_bottom):
    plt.plot([x1, x2, x3, x4, x1], [y1, y2, y3, y4, y1], marker='o', linestyle='-', color='b', label='Trapezoid')

    # Plotting the sides
    plt.plot([x1, x2], [y1, y2], linestyle='--', color='g', label='Bottom Side')
    plt.plot([x3, x4], [y3, y4], linestyle='--', color='m', label='Top Side')
    plt.plot([x1, x4], [y1, y4], linestyle='--', color='y', label='Left Side')
    plt.plot([x2, x3], [y2, y3], linestyle='--', color='c', label='Right Side')

    # Plotting the centroids of the sides
    centroid_bottom = ((x1 + x2) / 2, (y1 + y2) / 2)
    centroid_top = ((x3 + x4) / 2, (y3 + y4) / 2)
    centroid_left = ((x1 + x4) / 2, (y1 + y4) / 2)
    centroid_right = ((x2 + x3) / 2, (y2 + y3) / 2)

    plt.scatter(*centroid_bottom, color='g', marker='x', label='Bottom Centroid')
    plt.scatter(*centroid_top, color='m', marker='x', label='Top Centroid')
    plt.scatter(*centroid_left, color='y', marker='x', label='Left Centroid')
    plt.scatter(*centroid_right, color='c', marker='x', label='Right Centroid')

    centroid = trapezoid_centroid(x1, y1, x2, y2, x3, y3, x4, y4, thickness_left, thickness_right, thickness_top, thickness_bottom)
    plt.scatter(*centroid, color='r', label='Centroid')

    # Plotting the coordinates
    plt.text(x1, y1, '({}, {})'.format(x1, y1), fontsize=8, verticalalignment='bottom')
    plt.text(x2, y2, '({}, {})'.format(x2, y2), fontsize=8, verticalalignment='bottom')
    plt.text(x3, y3, '({}, {})'.format(x3, y3), fontsize=8, verticalalignment='top')
    plt.text(x4, y4, '({}, {})'.format(x4, y4), fontsize=8, verticalalignment='top')

    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Trapezoid and Centroid')
    plt.legend()
    plt.grid(True)
    plt.show()

# Example usage:
x1, y1 = map(float, input("Enter the coordinates of point 1 (x y): ").split())
x2, y2 = map(float, input("Enter the coordinates of point 2 (x y): ").split())
x3, y3 = map(float, input("Enter the coordinates of point 3 (x y): ").split())
x4, y4 = map(float, input("Enter the coordinates of point 4 (x y): ").split())

thickness_left = float(input("Enter the thickness of the left side: "))
thickness_right = float(input("Enter the thickness of the right side: "))
thickness_top = float(input("Enter the thickness of the top side: "))
thickness_bottom = float(input("Enter the thickness of the bottom side: "))

plot_trapezoid(x1, y1, x2, y2, x3, y3, x4, y4, thickness_left, thickness_right, thickness_top, thickness_bottom)
centroid = trapezoid_centroid(x1, y1, x2, y2, x3, y3, x4, y4, thickness_left, thickness_right, thickness_top, thickness_bottom)
print("Calculated Centroid:", centroid)

