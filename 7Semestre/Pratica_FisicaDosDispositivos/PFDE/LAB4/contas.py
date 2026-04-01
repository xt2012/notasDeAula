import numpy as np

V2 = 12

VR1 = np.array([45.3e-3, 0.985, 1.424, 1.725, 1.951, 2.14, 2.30])

VR2 = np.array([0, 0.5, 0.999, 1.501, 2.00, 2.50, 3.00])

R1 = 101.3
R2 = 997

Iled = VR2/R2

print("Iled=" , Iled)

Vldr = V2-VR1

print("Vldr=" ,Vldr)

Ildr = VR1/R1

print("Ildr=" , Ildr)

Rldr = Vldr/Ildr

print("Rldr=", Rldr)

Gldr = 1/Rldr

print("Gldr=", Gldr)
