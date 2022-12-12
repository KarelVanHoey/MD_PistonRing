# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 10:59:11 2020
@ProjectTitle : Thermal Mixed-Hydrodynamic Lubrication Solver for simulation of a compression ring - cylinder liner contact in an IC Engine.
@author: Dieter Fauconnier
@version : 1.0
"""


import numpy as np # Matrix Definitions and operations.
import scipy.sparse as sparse # Sparse Matrix Definitions and operations
import scipy.sparse.linalg as linalg # Sparse Matrix Linear Algebra
import matplotlib.pyplot as plt
import VisualLib as vis

class ReynoldsSolver:
    def __init__(self,Grid,Time,Ops,FluidModel,Discretization):
        self.MaxIter=[]
        self.TolP=[]
        self.UnderRelaxP=[]
        self.SetSolver()
        self.VisualFeedbackLevel=0
        
        self.Grid=Grid
        self.Time=Time
        self.Ops=Ops
        self.FluidModel=FluidModel
        self.Discretization=Discretization
        self.VisualFeedbackLevel=0
    
    def SetSolver(self,MaxIter: int=10000,TolP: float=1.0e-5 ,UnderRelaxP: float=0.001,TolT: float=1.0e-5 ,UnderRelaxT: float=0.001, VisualFeedbackLevel: int=0):
        self.MaxIter=MaxIter
        self.TolP=TolP
        self.UnderRelaxP=UnderRelaxP
        self.TolT=TolT
        self.UnderRelaxT=UnderRelaxT
        self.VisualFeedbackLevel=VisualFeedbackLevel

 
#################
##### TO DO #####
#################       
    def SolveReynolds(self,StateVector,time): # StateVector is both in and output
        
        #1. reset convergence
        epsP=np.zeros(self.MaxIter+1)
        epsP[0]=1.0
        epsT=np.zeros(self.MaxIter+1)
        epsT[0]=1.0
        
        
        #2. Predefine variables outside loop for Computational Efficiency
        DensityFunc    =self.FluidModel.Density
        ViscosityFunc  =self.FluidModel.DynamicViscosity
        SpecHeatFunc   =self.FluidModel.SpecificHeatCapacity
        ConducFunc     =self.FluidModel.ThermalConductivity      
        # PreviousDensity    =self.FluidModel.Density(StateVector[time-1]) ## Not used
        Density_prev = DensityFunc(StateVector[time-1])

        PHI = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        DPHIDX = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        DPDX = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        UPLUSDT = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        UMINDT = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        E1 = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        I = sparse.identity(self.Grid.Nx, dtype='float', format="csr")
        
        DDX=self.Discretization.DDXCentral
        DDXBackward=self.Discretization.DDXBackward
        DDXForward=self.Discretization.DDXForward
        D2DX2=self.Discretization.D2DX2
        SetDirichletLeft=self.Discretization.SetDirichletLeft
        SetDirichletRight=self.Discretization.SetDirichletRight
        # SetNeumannLeft=self.Discretization.SetNeumannLeft
        # SetNeumannRight=self.Discretization.SetNeumannRight
        
        #define your own when desired
        
        #3. Iterate

        k=0
        while ((((epsP[k]>self.TolP)) or (epsT[k]>self.TolT)) and (k<self.MaxIter)):
        
     
            #0. Calc Properties
            Density = DensityFunc(StateVector[time])
            SpecHeat = SpecHeatFunc(StateVector[time])
            Viscosity = ViscosityFunc(StateVector[time])
            Conduc = ConducFunc(StateVector[time])
            

            # phi = np.divide(np.multiply(Density, StateVector[time].h**3), 12 * Viscosity) # kan gwn normaal (numpy array)
            phi = Density * StateVector[time].h**3 / (12 * Viscosity)

            PHI.data = phi
            DPHIDX.data = DDX @ phi
        
            #1. LHS Pressure
            M = PHI @ D2DX2 + DPHIDX @ DDX
            
        
            #2. RHS Pressure
            #kan gwn normaal + squeeze term eerst uitzetten bij uittesten 
            # b = self.Ops.SlidingVelocity[time]/2 * DDX @ sparse.csc_matrix(np.matrix(np.multiply(Density, StateVector[time].h)).T) # + sparse.csc_matrix(np.matrix((np.multiply(Density, StateVector[time].h) - np.multiply(Density_prev, StateVector[time-1].h)) / self.Time.dt).T) 
            squeeze = (Density * StateVector[time].h - Density_prev * StateVector[time-1].h) / self.Time.dt
            b = self.Ops.SlidingVelocity[time]/2 * (DDX @ (Density * StateVector[time].h)) + squeeze
   
            #3. Set Boundary Conditions Pressure --> NOTE work with absolute pressure!!

            SetDirichletLeft(M) 
            SetDirichletRight(M)

            b[0] = self.Ops.AtmosphericPressure     
            b[-1] = self.Ops.CylinderPressure[time] 
            
            #4. Solve System for Pressure + Update
            p_star = linalg.spsolve(M,b)
            Delta_p = np.maximum(p_star, 0) - StateVector[time].Pressure #np.zeros(self.Grid.Nx) kan gwn nul zijn

            ## Update pressure
            StateVector[time].Pressure += self.UnderRelaxP * Delta_p
            # print("Pressure")
            # print(StateVector[time].Pressure)

            # ## Victor: Update properties dependent on pressure --> moet dit? Staat niet in algoritme vermeld. Kunnen desnoods eens testen of dit sneller tot convergentie leidt.
            # Density = DensityFunc(StateVector[time])
            # SpecHeat = SpecHeatFunc(StateVector[time])
            # Viscosity = ViscosityFunc(StateVector[time])
            # Conduc = ConducFunc(StateVector[time])
            
            #5. LHS Temperature ---> Absolute temperaturen!
            #Uaveraged = 0

            av_u = - StateVector[time].h**2 / (12 * Viscosity) * (DDX @ StateVector[time].Pressure) + self.Ops.SlidingVelocity[time] / 2
            Uaveraged = av_u

            # u_plus = np.maximum(av_u, 0)
            # u_min = np.minimum(av_u, 0)
            
            # UPLUSDT.data[:] = u_plus * self.Time.dt
            # UMINDT.data[:] = u_min * self.Time.dt
            # E1.data[:] = Conduc / (Density * SpecHeat

            D = UPLUSDT @ DDXBackward + UMINDT @ DDXForward
            E = - E1 @ D2DX2
            M1 = I + D + E

            # #6. RHS Temperature
            Q = StateVector[time].h**2 / (12 * Viscosity) * (DDX @ StateVector[time].Pressure)**2 + Viscosity *( self.Ops.SlidingVelocity[time]**2 / StateVector[time].h**2 )
            RHS = StateVector[time-1].Temperature + self.Time.dt * Q / (Density * SpecHeat)

            # #Boundary conditions
            # if self.Ops.SlidingVelocity[time] <= 0:
            #     M1[0,0:1] = [-1/self.Grid.dx, 1/self.Grid.dx] ### Karel: hier moet het denk ik Engine.CompressionRing.Thickness/self.Grid.Nx zijn.
            #     M1[-1, -1] = 1
            #     M1[0,3:] = 0   ### Karel: niet juist denk ik, alle getallen vanaf 3 in die rij moeten nul zijn --> [0, 3:] (zie voorbeeldje in Test_Sparse_Matrix.py)
            #     M1[-1,1:-1] = 0  ### Karel: hier moet het [-1, 1:-1] zijn denk ik
            #     RHS[0] = 0
            #     RHS[-1] = self.Ops.OilTemperature
            # else:
            #     M1[0,0] = 1     ## Nog eens checken!
            #     M1[2:, 0] = 0
            #     M1[-1,-2:] = [-1/self.Grid.dx, 1/self.Grid.dx]
            #     M1[-1,1:-2] = 0
            #     RHS[0] = self.Ops.OilTemperature
            #     RHS[-1] = 0
            # #7. Solve System for Temperature + Update

            T_star = linalg.spsolve(M1, RHS)
            delta_T = T_star - StateVector[time].Temperature
            StateVector[time].Temperature += delta_T * self.UnderRelaxT

            #8. Calculate other quantities: Hydrodynamic load (eq. 37 in assignment), Wall shear stress, Viscous friction force (store all in StateVector)
            #################################################################################################################
            ### Opmerking: moet dit nu hier of moet dit pas in stap 11? --> lijkt mij nutteloos om dit hier te berekenen. ###
            #################################################################################################################

            
            #9. Residuals & Report
            k += 1

            epsP[k] = np.linalg.norm(np.divide(Delta_p, StateVector[time].Pressure)) / self.Grid.Nx
            # epsT[k] = np.linalg.norm(np.divide(delta_T, StateVector[time].Temperature)) / self.Grid.Nx


           
            #10. Provide a plot of the solution
            if (k % 500 == 0):
                CFL=np.max(Uaveraged)*self.Time.dt/self.Grid.dx
                print("ReynoldsSolver:: CFL", np.round(CFL,2) ,"Residual [P,T] @Time:",round(self.Time.t[time]*1000,5),"ms & Iteration:",k,"-> [",np.round(epsP[k],6),",",np.round(epsT[k],6),"]")
                if self.VisualFeedbackLevel>2:
                    fig=vis.Report_PT(self.Grid,StateVector[time]) 
                    plt.close(fig)

                
            if (epsP[k]<=self.TolP) and (epsT[k]<=self.TolT):
                print("ReynoldsSolver:: Convergence [P,T] to the predefined tolerance @Time:",round(self.Time.t[time]*1000,5),"ms & Iteration:",k,"-> [",np.round(epsP[k],6),",",np.round(epsT[k],6),"]")
                
            if k>=self.MaxIter:
                print("ReynoldsSolver:: Residual [P,T] @Time:",round(self.Time.t[time]*1000,5),"ms & Iteration:",k,"-> [",np.round(epsP[k],6),",",np.round(epsT[k],6),"]")
                print("ReynoldsSolver:: Warning: Maximum Iterations without converging to the predefined tolerance]")


            
        #11. Calculate other quantities (e.g. Wall Shear Stress, Hydrodynamic Load, ViscousFriction)
        StateVector[time].HydrodynamicLoad = np.trapz(StateVector[time].Pressure, dx=self.Grid.dx)
        # print("Help2")
        # print( StateVector[time].HydrodynamicLoad)
        WallShearStress_0 = Viscosity * ((self.Ops.SlidingVelocity[time-1]-self.Ops.SlidingVelocity[time])/StateVector[time].h) - (DDX @ StateVector[time].Pressure)*(StateVector[time].h/2) #Uit cursus gehaald, idk of dit correct is...
        WallShearStress_h = Viscosity * ((self.Ops.SlidingVelocity[time-1]-self.Ops.SlidingVelocity[time])/StateVector[time].h) - (DDX @ StateVector[time].Pressure)*(StateVector[time].h/2)
        StateVector[time].WallShearStress = WallShearStress_h


