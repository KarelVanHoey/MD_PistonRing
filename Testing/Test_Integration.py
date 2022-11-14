import scipy.integrate as spint
import numpy as np

# def F1(input):
#     return input ** 2 + 1

# int1 = spint.quad(F1, 0, 4)
# print(int1[0])

# F2 = lambda x: x**2 + 1

# int2 = spint.quad(F2, 0, 4)
# print(int2[0])


y = np.array([0, 1, 2, 3, 4])
x = y
print(spint.cumtrapz(y, x)) # Contains integral up to each timestep; element n is the integral until x[n+1]
print((y**2) / 2)