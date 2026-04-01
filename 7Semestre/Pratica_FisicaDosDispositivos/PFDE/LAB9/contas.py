import numpy as np

R1 = 1015
R2 = 989e3

VDC = np.array([0, 5, 10, 15, 20, 25])

V1 = np.array([0, 3.18, 8.04, 12.96, 17.85, 22.5])

V2 = np.array([0, 0.422, 1.076, 1.632, 2.08, 2.37])

Vled1 = VDC-V1

VDC2 = 12

print('======= Parte 1 ======')

print('Vled1=', Vled1)

Iled1 = V1/R1
print('Iled1=', Iled1)

Vled2 = VDC2-V2
print('Vled2=', Vled2)

Iled2 = V2/R2
print('Iled2=', Iled2)


print('======= Parte 2 =======')


V1 = np.array([0, 3.2, 8.0, 12.95, 17.8, 22.9])

V2 = np.array([0, -0.372, -0.902, -1.291, -1.404, -1.435])

Vled1 = VDC-V1

print('Vled1=', Vled1)

Iled1 = V1/R1
print('Iled1=', Iled1)

Vled2 = V2
print('Vled2=', Vled2)

Iled2 = V2/R2
print('Iled2=', Iled2)


print('========== Parte 3 ==========')

R2 = np.array([19.99e3, 0.498e6, 1.001e6, 1.499e6, 2e6, 2.5e6, 3e6, 3.5e6, 4.01e6, 4.51e6, 5e6, 9.96e6])

Vled2 = np.array([-38.8e-3, -885e-3, -1.404, -1.460, -1.473, -1.477, -1.490, -1.493, -1.494, -1.494, -1.499, -1.504])


Iled2 = Vled2/R2

print('Iled2=', Iled2)
