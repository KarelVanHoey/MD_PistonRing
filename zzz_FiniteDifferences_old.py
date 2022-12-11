# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 15:03:02 2020
@ProjectTitle : Thermal Mixed-Hydrodynamic Lubrication Solver for simulation of a compression ring - cylinder liner contact in an IC Engine.
@author: Dieter Fauconnier
@version : 1.0
"""

from webbrowser import BackgroundBrowser
import numpy as np # Matrix Definitions and operations.
import scipy.sparse as sparse # Sparse Matrix Definitions and operations


class FiniteDifferences:
    
    
#################
##### TO DO #####
#################
    def __init__(self,Grid):
        
        #Schemes for first-order (D) and second order (DD) derivatives
        self.CentralStencilD=np.array([-1.0,0.0,1.0])/(2.0*Grid.dx)
        self.BackwardStencilD=np.array([-1.0,1.0])/Grid.dx
        self.ForwardStencilD=np.array([-1.0,1.0])/Grid.dx
        self.BackwardStencilD2=np.array([0.0,-1.0,1.0])/Grid.dx
        self.ForwardStencilD2=np.array([-1.0,1.0,0.0])/Grid.dx
        
        self.CentralStencilDD=np.array([1.0, -2.0, 1.0])/(Grid.dx**2)

        
        # Sparse Matrix Operators
        self.Identity=sparse.identity(Grid.Nx, dtype='float', format="csr")

        # Because the domain is not periodic (x begins and ends at piston ring end), the best approximation of the derivatives
        # is given by the forward or backward approach.
        # There is no variation in the circumferential direction.

        # First derivatives (DDX)
        self.DDXCentral=sparse.diags(self.CentralStencilD, [-1, 0, 1], shape=(Grid.Nx, Grid.Nx), dtype='float', format="csr")
        self.DDXCentral[-1,-2:]=self.BackwardStencilD #Define right boundary stencil
        self.DDXCentral[0,:2]=self.ForwardStencilD #Define left boundary stencil

        self.DDXBackward=sparse.diags(self.BackwardStencilD, [-1, 0], shape=(Grid.Nx, Grid.Nx), dtype='float', format='csr')
        self.DDXBackward[0,:2]=self.ForwardStencilD

        self.DDXForward=sparse.diags(self.ForwardStencilD, [0, 1], shape=(Grid.Nx, Grid.Nx), dtype='float', format='csr')
        self.DDXForward[-1,-2:]=self.BackwardStencilD
        
        # Second Derivative (D2DX2)
        self.D2DX2=sparse.diags(self.CentralStencilDD, [-1, 0, 1], shape=(Grid.Nx, Grid.Nx), dtype='float', format='csr')
        self.D2DX2[0,:3]=self.CentralStencilDD #Define right boundary stencil
        self.D2DX2[-1,-3:]=self.CentralStencilDD #Define left boundary stencil
        

    # Efficient Implementation for 1D csr type FD matrix
    #do not Change implementation below!
    def SetDirichletLeft(self,M):
        M.data[[0,1,2]]=[1.0, 0.0, 0.0]
        

    def SetDirichletRight(self,M):
        M.data[[-3,-2,-1]]=[0.0, 0.0, 1.0]
    
    def SetNeumannLeft(self,M):
        M.data[[0,1,2]]=self.ForwardStencilD2
    
    def SetNeumannRight(self,M):
        M.data[[-3,-2,-1]]=self.BackwardStencilD2
