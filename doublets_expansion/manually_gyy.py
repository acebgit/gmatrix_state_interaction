import numpy as np

a = 0.02404144
b = -0.99970303
c = 3.98140665e-03j
d = -9.57472014e-05j
ge = 2.002319304363
gyy_reference = -0.064

A = np.conj(b) * b - np.conj(-a) * a + np.conj(-d) * d - np.conj(c) * c
B = np.conj(a) * b - np.conj(b) * a + np.conj(c) * d - np.conj(d) * c

Gyy = ge**2 * ( (A.real)**2  + (A.imag)**2 + B**2 )
gyy = (np.sqrt(Gyy) - ge) *1000
print(gyy)