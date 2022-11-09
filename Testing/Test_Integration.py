import scipy.integrate as spint

def F1(input):
    return input ** 2 + 1

int1 = spint.quad(F1, 0, 4)
print(int1[0])

F2 = lambda x: x**2 + 1

int2 = spint.quad(F2, 0, 4)
print(int2[0])