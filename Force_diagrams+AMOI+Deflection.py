import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
#import initializations

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

# Multiply the total moment of inertia by the array [10.22960, -0.2216748821]
multiplied_moment_of_inertia = np.array([1.09259970e+03, -7.05468577e+01, -6.84947548e+00,
        -1.08929900e-02]) * moment_of_inertia_new_axis
print(f"Multiplied Moment of Inertia Array: {multiplied_moment_of_inertia}")

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


y0lst = np.zeros(0)
a_i0lst = np.zeros(0)
Cl0lst = np.zeros(0)
Cd0lst = np.zeros(0)
Cm0lst = np.zeros(0)
sw = 0
fuel = 1

MTOW = 325100*9.81

if fuel == 0:
    MTOW -= 126800*9.81 #initially 139170

plf = -1
nlf = -1
CL_d = -0.5
L_plf = plf*MTOW
L_nlf = nlf*MTOW
rho = 1.225225
S = 435
#V = 350
thrust = 0
moment_arm = 2.4 #from catia

y_second_lst = [0, 33.226]
clst = [10.23, 2.864]

CLlst = []
lst1 = []
lst2 = []



with open('a0.txt', 'r') as fin:
    for line in fin:
        line = line.strip()
        if line != '':
            if line.find('CL') != -1:
                CLlst = line.strip()
                CLlst = CLlst.split('=')
                CL0 = float(CLlst[1].strip())
            if line.find('Coefficients') != -1:
                sw = 0
            if sw == 1:
                line = line.strip()
                lst1 = line.split(' ')
                if lst1[0][0] != '-':
                    sw = 2
            if sw == 2:
                line = line.strip()
                lst1 = line.split(' ')
                for elmnt in lst1:
                    if elmnt != '':
                        lst2.append(elmnt)
                y0lst = np.append(y0lst, float(lst2[0]))
                a_i0lst = np.append(a_i0lst, float(lst2[2]))
                Cl0lst = np.append(Cl0lst, float(lst2[3]))
                Cd0lst = np.append(Cd0lst, float(lst2[5]))
                Cm0lst = np.append(Cm0lst, float(lst2[7]))
                lst2 = []
            if line.find('y-span') != -1:
                sw = 1


Cl0_y = sp.interpolate.interp1d(y0lst,Cl0lst,kind='cubic',fill_value="extrapolate")
a_i0_y = sp.interpolate.interp1d(y0lst,a_i0lst,kind='cubic',fill_value="extrapolate")
Cd0_y = sp.interpolate.interp1d(y0lst,Cd0lst,kind='cubic',fill_value="extrapolate")
Cm0_y = sp.interpolate.interp1d(y0lst,Cm0lst,kind='cubic',fill_value="extrapolate")

c_y = sp.interpolate.interp1d(y_second_lst,clst,kind='linear',fill_value="extrapolate")


y10lst = np.zeros(0)
a_i10lst = np.zeros(0)
Cl10lst = np.zeros(0)
Cd10lst = np.zeros(0)
Cm10lst = np.zeros(0)
sw = 0

with open('a10.txt', 'r') as fin:
    for line in fin:
        line = line.strip()
        if line != '':
            if line.find('CL') != -1:
                CLlst = line.strip()
                CLlst = CLlst.split('=')
                CL10 = float(CLlst[1].strip())
            if line.find('Coefficients') != -1:
                sw = 0
            if sw == 1:
                line = line.strip()
                lst1 = line.split(' ')
                if lst1[0][0] != '-':
                    sw = 2
            if sw == 2:
                line = line.strip()
                lst1 = line.split(' ')
                for elmnt in lst1:
                    if elmnt != '':
                        lst2.append(elmnt)
                y10lst = np.append (y10lst, float(lst2[0]))
                a_i10lst = np.append(a_i10lst, float(lst2[2]))
                Cl10lst = np.append (Cl10lst, float(lst2[3]))
                Cd10lst = np.append (Cd10lst, float(lst2[5]))
                Cm10lst = np.append (Cm10lst, float(lst2[7]))
                lst2 = []
            if line.find('y-span') != -1:
                sw = 1


Cl10_y = sp.interpolate.interp1d(y10lst,Cl10lst,kind='cubic',fill_value="extrapolate")
a_i10_y = sp.interpolate.interp1d(y0lst,a_i10lst,kind='cubic',fill_value="extrapolate")
Cd10_y = sp.interpolate.interp1d(y10lst,Cd10lst,kind='cubic',fill_value="extrapolate")
Cm10_y = sp.interpolate.interp1d(y10lst,Cm10lst,kind='cubic',fill_value="extrapolate")

a_lst = [0, 10]
a_i_minilst = [0, 0]



    

"""def interpolate_lists(list1, list2, t):
    
    Linear interpolation between two lists.

    Parameters:
    - list1, list2: Lists containing the same number of elements.
    - t: Interpolation parameter (between 0 and 1).

    Returns:
    - Interpolated list.
    
    interpolated_list = [a + t * (b - a) for a, b in zip(list1, list2)]
    return interpolated_list"""

#a_d = 180*np.arcsin((CL_d-CL0)*np.sin(10*np.pi/180)/(CL10-CL0))/np.pi #a_d stored in degrees

"""def spanwise_distr_graphs(y):
    Cl_d_y = Cl0_y(y) + (CL_d-CL0)*(Cl10_y(y)-Cl0_y(y))/(CL10-CL0)
    L_y = Cl_d_y*0.5*(V**2)*rho*c_y(y)
    a_i_a = interpolate_lists(a_i0_y(y), a_i10_y(y), a_d/10)
    Cd_d_y = Cl_d_y*np.pi*a_i_a/180
    D_y = Cd_d_y*0.5*(V**2)*rho*c_y(y)
    Cm_qc_d_y = Cm0_y(y) + (CL_d-CL0)*(Cm10_y(y)-Cm0_y(y))/(CL10-CL0)
    M_qc_y = Cm_qc_d_y*0.5*(V**2)*rho*(c_y(y))**2
    N_y = L_y*np.cos(np.pi*a_d/180)+D_y*np.sin(np.pi*a_d/180)
    T_y = D_y*np.cos(np.pi*a_d/180)-L_y*np.sin(np.pi*a_d/180)
    M_04c_y = M_qc_y+0.15*N_y
    W_y = 37.6*c_y(y)**1.5
    W_f_y = 433.56*(c_y(y))**2
    return N_y, T_y, M_04c_y, W_y"""

def spanwise_distr_root_forces(y):
    Cl_d_y = Cl0_y(y) + (CL_d-CL0)*(Cl10_y(y)-Cl0_y(y))/(CL10-CL0)
    #a_d = 180*np.arcsin((CL_d-CL0)*np.sin(10*np.pi/180)/(CL10-CL0))/np.pi #a_d stored in degrees
    L_y = Cl_d_y*0.5*(V**2)*rho*c_y(y)
    a_i_minilst[0] = a_i0_y(y)
    a_i_minilst[1] = a_i10_y(y)
    a_i_a_y = sp.interpolate.interp1d(a_lst,a_i_minilst,kind='linear',fill_value="extrapolate")
    Cd_d_y = Cl_d_y*np.pi*a_i_a_y(a_d)/180
    D_y = Cd_d_y*0.5*(V**2)*rho*c_y(y)
    Cm_qc_d_y = Cm0_y(y) + (CL_d-CL0)*(Cm10_y(y)-Cm0_y(y))/(CL10-CL0)
    #Cm_qc_d_y = -0.1
    M_qc_y = Cm_qc_d_y*0.5*(V**2)*rho*(c_y(y))**2
    #print (M_qc_y)
    N_y = L_y*np.cos(np.pi*a_d/180)+D_y*np.sin(np.pi*a_d/180)
    T_y = D_y*np.cos(np.pi*a_d/180)-L_y*np.sin(np.pi*a_d/180)
    M_04c_y = M_qc_y+0.15*N_y*c_y(y)
    #print (M_04c_y)
    W_y = 9.81*37.6*c_y(y)**1.5
    W_f_y = 425.75204*(c_y(y))**2
    return L_y, D_y, M_04c_y, W_y, W_f_y


def L_distr(x):
    L, D, M_04c, W, W_f = spanwise_distr_root_forces(x)
    return L

def D_distr(x):
    L, D, M_04c, W, W_f = spanwise_distr_root_forces(x)
    return D

def W_distr(x):
    L, D, M_04c, W, W_f = spanwise_distr_root_forces(x)
    return W

def M_04c_distr(x):
    L, D, M_04c, W, W_f = spanwise_distr_root_forces(x)
    return M_04c

def W_f_distr(x):
    L, D, M_04c, W, W_f = spanwise_distr_root_forces(x)
    return W_f

    

Y = np.linspace (0, 33.226, 50)
shear_normal = np.zeros(50)
shear_axial = np.zeros(50)
moment = np.zeros(50)
torque = np.zeros(50)

#CL_d = 2*L_plf/(rho*V**2*S)
V = np.sqrt(abs(2*L_nlf/(rho*CL_d*S)))
a_d = 180*np.arcsin((CL_d-CL0)*np.sin(10*np.pi/180)/(CL10-CL0))/np.pi
#print (a_d, V)

for i in range (49, -1, -1):
    S_lift, err_S_lift = sp.integrate.quad(L_distr, Y[i], 33.226)
    S_weight, err_S_weight = sp.integrate.quad(W_distr, Y[i], 33.226)
    S_weight_fuel, err_S_weight_fuel = sp.integrate.quad(W_f_distr, Y[i], 33.226)
    S_drag, err_S_drag = sp.integrate.quad(D_distr, Y[i], 33.226)
    Torque, err_Torque = sp.integrate.quad(M_04c_distr, Y[i], 33.226)

    S_weight *= plf
    S_weight_fuel *= plf
    
    if fuel == 1:
        S_normal = (S_lift-S_weight-S_weight_fuel)*np.cos(np.pi*a_d/180)+S_drag*np.sin(np.pi*a_d/180)
        S_axial = -(S_lift-S_weight-S_weight_fuel)*np.sin(np.pi*a_d/180)+S_drag*np.cos(np.pi*a_d/180)
    else:
        S_normal = (S_lift-S_weight)*np.cos(np.pi*a_d/180)+S_drag*np.sin(np.pi*a_d/180)
        S_axial = -(S_lift-S_weight)*np.sin(np.pi*a_d/180)+S_drag*np.cos(np.pi*a_d/180)
    
    shear_normal[i] = S_normal
    shear_axial[i] = S_axial
    torque[i] = Torque
    if Y[i] < 5.343:#L/G
        shear_normal[i] -= plf*3028*9.81*np.cos(np.pi*a_d/180)#L/G
    if Y[i] < 11.63:#engine
        shear_normal[i] -= plf*9630*9.81*np.cos(np.pi*a_d/180)#engine
        shear_axial[i] += thrust*np.cos(np.pi*34.4/180)#thrust vector
        torque[i] -= moment_arm*thrust*np.cos(np.pi*34.4/180)
    

def internal_normal_distr(x):
    internal_normal = sp.interpolate.interp1d(Y,shear_normal,kind='cubic',fill_value="extrapolate")
    return internal_normal(x)

for i in range (49, -1, -1):
    Moment, err_Moment = sp.integrate.quad(internal_normal_distr, Y[i], 33.226)
    moment[i] = Moment
    if Y[i] < 11.63:
        moment[i] = moment[i] - thrust*moment_arm*np.sin(np.pi*34.4/180)
        
#print (shear_normal[0])
#print (torque[0])
#print (moment[0])
#print (shear_axial[0])

#moment_function = sp.interpolate.interp1d(Y,moment,kind='cubic',fill_value="extrapolate")
#print (moment)

with open('negative_OEW+Payload+fuel.txt', 'w') as fout:
    for i in range(0, 50):
        fout.write(f'{shear_normal[i]}\n')
    for i in range(0, 50):
        fout.write(f'{torque[i]}\n')
    for i in range(0, 50):
        fout.write(f'{moment[i]}\n')

s=0
for i in range(33):
    print 
    print (M_04c_distr(i), c_y(i))
    s+= M_04c_distr(i)

#print (s)

plt.subplot(221)
plt.plot(Y, shear_normal, color="black")
plt.title("Internal shear diagram perpendicular to chord")

plt.subplot(222)
plt.plot(Y, shear_axial, color="black")
plt.title("Internal shear diagram parallel to chord")

plt.subplot(223)
plt.plot(Y, moment, color="black")
plt.title("Bending moment diagram")

plt.subplot(224)
plt.plot(Y, torque, color="black")
plt.title("Torque diagram at 0.4c")

plt.show()



#Olliver's part
E=68.9*10**9
I_tab=[]

def I(y):
    I=multiplied_moment_of_inertia[0]+multiplied_moment_of_inertia [1]*y+multiplied_moment_of_inertia [2]*y**2+multiplied_moment_of_inertia [3]*y**3
    return I

slope_function=np.zeros(50)
deflection=np.zeros(50)

def slope_integrand(y):
    slope_integrand=-moment_function(y)/E/I(y)
    return slope_integrand

for i in range(0,50,1):
    slope, err_slope = sp.integrate.quad(slope_integrand, 0,Y[i])
    slope_function[i] = slope

def deflection_integrand(x):
    deflection_integrand = sp.interpolate.interp1d(Y,slope_function,kind='cubic',fill_value="extrapolate")
    return deflection_integrand(x)

for i in range(0,50,1):
    Deflection, err_Deflection=sp.integrate.quad(deflection_integrand, 0, Y[i])
    deflection[i]=Deflection

plt.plot(Y, deflection)
plt.title("Deflcetion")
plt.show()

I_tab=[]
for i in range(50):
    I_tab.append(I(Y[i]))

plt.plot(Y, I_tab)
plt.show()
