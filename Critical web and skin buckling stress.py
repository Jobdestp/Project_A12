import math

a_0 = float(input('Initial longer side of the plate length: '))
b_0 = float(input('Initial thickness of the plate: '))
t_spar = float(input('Spar thickness: '))
l = float(input('Length of wingspan(one side): '))
gradient_a = float(input('Change of longer side of the plate along wingspan: '))
gradient_b = float(input('Change of shorter side of the plate along wingspan: '))
E = 68.9 * 10**9
v = 0.33
lst_aspect_ratio=[]
lst_web_buckling_coefficient=[]
lst_skin_buckling_coefficient=[]
lst_web_critical_stress=[]
lst_skin_critical_stress=[]


num_sections = 1000
wingspan_length = l
section_steps = l / num_sections

a = a_0
b = b_0


for i in range(0, num_sections):
    # Length of a and b at different sections
    a -= gradient_a * section_steps * i
    b -= gradient_b * section_steps * i

    # Calculate aspect ratio and buckling coefficient
    k = a / b
    k_s = 5.34 + 4 / ( k ** 2 )
    k_c =

    # Calculate critical web buckling stress and skin buckling stress
    web_critical_buckling_stress = ((math.pi ** 2 * k_s * E) / (12 * (1 - v ** 2))) * ( t_spar / b ) ** 2
    skin_critical_buckling_stress = ((math.pi ** 2 * k_c * E) / (12 * (1 - v ** 2))) * ( t_spar / b ) ** 2
    lst_web_critical_stress.append(web_critical_buckling_stress)
    lst_skin_critical_stress.append(skin_critical_buckling_stress)
    lst_web_buckling_coefficient.append( k_s )
    lst_skin_buckling_coefficient.append( k_c)


print('The web critical buckling stress is: ',lst_web_critical_stress)
print('The skin critical buckling stress is: ',lst_skin_critical_stress)




