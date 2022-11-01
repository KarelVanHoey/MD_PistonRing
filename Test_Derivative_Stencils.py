import numpy as np
import scipy.sparse as sparse

dx = 1.0
Nx = 6
CentralStencilD=np.array([-1.0,0.0,1.0])/(2.0*dx)
BackwardStencilD=np.array([-1.0,1.0])/dx
ForwardStencilD=np.array([-1.0,1.0])/dx
CentralStencilDD=np.array([1.0, -2.0, 1.0])/(dx**2)


# Central DDX
M = sparse.diags(CentralStencilD, [-1, 0, 1], shape=(Nx,Nx), dtype='float', format='csr')
# print(M.toarray())

M[-1,-2:]=BackwardStencilD
# print(M.toarray())

M[0,:2]=ForwardStencilD
# print(M.toarray())

# Backward DDX
N = sparse.diags(BackwardStencilD, [-1, 0], shape=(Nx,Nx), dtype='float', format='csr')
# print(N.toarray())

N[0,:2] = ForwardStencilD
# print(N.toarray())

# Forward DDX
I = sparse.diags(ForwardStencilD, [0, 1], shape=(Nx,Nx), dtype='float', format='csr')
# print(I.toarray())

I[-1,-2:]=BackwardStencilD
# print(I.toarray())

# Central D2DX2
K = sparse.diags(CentralStencilDD, [-1,0,1], shape=[Nx,Nx], dtype='float', format='csr')
# print(K.toarray())

K[0,:3]=CentralStencilDD
# print(K.toarray())

K[-1,-3:]=CentralStencilDD
# print(K.toarray())