import numpy as np

Va = np.array([0, 0.503, 0.987, 1.506, 2, 2.5, 3, 3.49, 4.01])

Vb = np.array([0, 0.325, 0.614, 0.976, 1.374, 1.807, 2.26, 2.77, 3.29])

R2=9.6

I = Vb/R2

print("I=",I)

Vab = np.subtract(Va,Vb)

print("Vab", Vab)

R1 = Vab/I

print("R1=", R1)

termal_dis=0.011

power=Vab*I

temp = power/termal_dis

print("temp=", temp)
