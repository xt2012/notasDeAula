import numpy as np


Va = [ 2, 4, 6, 8, 10, 12]
Vb = [1.88, 3.85, 5.81, 7.79, 9.77, 11.74] #tensão sobre o resistor

Vab = np.subtract(Va,Vb)

print(Vab)

I = Vab

R1 = Vb/I

print(R1)

P1 = Vb*I

print(P1)



# analise da potencia irradiada pelo sol 

A = 6.0877e12
T = 5800
sigma = 5.67e-8

Pr = sigma*A*T**4

print(Pr)
