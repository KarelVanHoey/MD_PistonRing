import numpy as np

a = np.linspace(0,0.001,21)
print("a =", a)
b = np.round(a*1000, 5)
print("b =", b)
for i in range(len(b)):
    # print('{0:.2f}'.format(i))
    print("Figures/PT@Time_"+"{0:.2f}".format(round(a[i]*1000,5))+"ms.png")