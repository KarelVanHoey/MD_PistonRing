import numpy as np
import matplotlib.pyplot as plt

# a = np.linspace(0,0.001,21)
# print("a =", a)
# b = np.round(a*1000, 5)
# print("b =", b)
# for i in range(len(b)):
#     # print('{0:.2f}'.format(i))
#     print("Figures/PT@Time_"+"{0:.2f}".format(round(a[i]*1000,5))+"ms.png")

# a = np.arange(10)
# print(a)
# a += 1
# print(a)

CrownHeight = 1e-5
Thickness = 1.5e-3
x = np.linspace(-Thickness/2, Thickness/2, 100)
h0 = 5e-5

# a = np.array([min(.75*CrownHeight, i) for i in 4 * CrownHeight * (x**2) / (Thickness**2)]) + h0

h_orig = 4 * CrownHeight * (x**2) / (Thickness**2)
h_transl = h_orig - 0.25 * CrownHeight
h_capped = np.maximum(h_transl, 0)

plt.plot(x, h_orig, x, h_transl, x, h_capped)
plt.show()