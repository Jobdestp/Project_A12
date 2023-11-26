import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

# Given data
e = 68.9 * 10**9
y = np.array([0.0, 0.6768367346938775, 1.353673469387755, 2.0305102040816325, 2.70734693877551, 3.384183673469388, 4.061020408163266, 4.737857142857144, 5.414693877551021, 6.091530612244898, 6.768367346938776, 7.445204081632653, 8.12204081632653, 8.798877551020408, 9.475714285714286, 10.152551020408163, 10.829387755102041, 11.506224489795918, 12.183061224489796, 12.859897959183673, 13.53673469387755, 14.213571428571428, 14.890408163265305, 15.567244897959183, 16.24408163265306, 16.920918367346937, 17.597755102040817, 18.274591836734694, 18.95142857142857, 19.628265306122447, 20.305102040816325, 20.981938775510204, 21.65877551020408, 22.335612244897958, 23.012448979591837, 23.689285714285715, 24.366122448979594, 25.04295918367347, 25.71979591836735, 26.396632653061224, 27.073469387755105, 27.75030612244898, 28.42714285714286, 29.10397959183674, 29.780816326530615, 30.457653061224494, 31.13448979591837, 31.81132653061225, 32.48816326530613, 33.225])

m = np.array([4.47194376e+07, 4.27872239e+07, 4.08950227e+07, 3.90441701e+07,
 3.72357910e+07, 3.54710653e+07, 3.37503060e+07, 3.20764547e+07,
 3.04370474e+07, 2.88326461e+07, 2.72768474e+07, 2.57678728e+07,
 2.43066073e+07, 2.28931718e+07, 2.15276605e+07, 2.02106449e+07,
 1.89404458e+07, 1.77235603e+07, 1.66290946e+07, 1.54399973e+07,
 1.43039203e+07, 1.32141744e+07, 1.21722627e+07, 1.11774384e+07,
 1.02295043e+07, 9.32806429e+06, 8.47272515e+06, 7.66304344e+06,
 6.89854075e+06, 6.17870036e+06, 5.50296881e+06, 4.87075583e+06,
 4.28143445e+06, 3.73434093e+06, 3.22877430e+06, 2.76399552e+06,
 2.33922645e+06, 1.95364860e+06, 1.60640106e+06, 1.29657810e+06,
 1.02322310e+06, 7.85320600e+05, 5.81788415e+05, 4.11470291e+05,
 2.73113424e+05, 1.65305672e+05, 8.64423845e+04, 3.46362639e+04,
 7.39957469e+03, 0.00000000e+00])

# Interpolate the given data
g = interp1d(y, m, kind='cubic', fill_value="extrapolate")

# Define the moment of inertia function
def moi_function(y):
#    return 5.68905490e-01 -2.46562931e-02*y + 2.67150169e-04*y**2
#    return 3.79270326e-01 -1.64375287e-02*y + 1.78100112e-04*y**2
#    return 1.89635163e+00 -8.21876437e-02*y + 8.90500562e-04*y**2
    return 3.79270326e+00 -1.64375287e-01*y + 1.78100112e-03*y**2
#     return 9.48175816e+00 -4.10938219e-01*y + 4.45250281e-03*y**2
#    return 1.89635163e+01 -8.21876437e-01*y + 8.90500562e-03*y**2
#1.89635163e2 -8.21876437e1*y + 8.90500562e0*y**2

# Define the integrand function
def integrand_function(y):
    return g(y) / (e * moi_function(y))

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
plt.plot(y, n, label='n(y)')
plt.plot(y, p, label='p(y)')
plt.title('Plot of n(y) and p(y)')
plt.xlabel('y')
plt.ylabel('Function Values')
plt.legend()
plt.grid(True)
plt.show()
