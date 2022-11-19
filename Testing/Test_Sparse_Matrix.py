import scipy.sparse as sparse
import numpy as np

# a = np.array([[1.0, 2.0, 3.0]])
# print(np.shape(np.transpose(a)))
# b = np.array([1.0, 2.0, 3.0])

# I = sparse.diags(np.ones(3))

# # print(I @ a)
# print(I @ b)

# Implementation now:
# phi = np.arange(4)
# phi = sparse.csr_matrix(phi)
# # print(phi)
# eye = np.identity(4)

# print(phi.multiply(eye))

# print('----------')
# In assignment:
# phi = np.arange(4)
# phi = sparse.diags(phi)
# # print(phi)
# eye = np.identity(4)

# print(phi.multiply(eye))

#### OK: they give the same result

a = np.array([[1, 2], [3, 4]])
print('a=',a)
b = sparse.diags([2, 2])
print('b=',b.toarray())
c = sparse.csr_matrix(np.matrix([1, 2]).T)
print('c=',c.toarray())

print(b @ c)
# print(sparse.csc_matrix(np.matrix(np.arange(4)).T).toarray())