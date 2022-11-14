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
        PreviousDensity    =self.FluidModel.Density(StateVector[time-1])
        
        DDX=self.Discretization.DDXCentral
        # DDXBackward=self.Discretization.DDXBackward
        # DDXForward=self.Discretization.DDXForward
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
            Density_prev = PreviousDensity() #Needed to calculate time differentiation in step 2.
            SpecHeat = SpecHeatFunc(StateVector[time])
            Viscosity = ViscosityFunc(StateVector[time])
            Conduc = ConducFunc(StateVector[time])

            phi = Density*StateVector[time].h**3 / (12 * Viscosity)
        
            #1. LHS Pressure

            M = phi * D2DX2 + DDX @ phi @ DDX
        
            #2. RHS Pressure
            
            b = self.Ops.SlidingVelocity/2 * DDX @ (Density * StateVector[time].h) + (Density * StateVector[time].h - Density_prev * StateVector[time-1].h) / self.Time.dt 
                #Note: Squeeze term: backward time differentiation for d(rho*h)/dt!
   
            #3. Set Boundary Conditions Pressure --> NOTE work with absolute pressure!!
            SetDirichletRight(M)
            SetDirichletLeft(M)

            b[0] = self.Ops.AtmosphericPressure     # [psi] --> Moet hier gedefinieerd worden waar in de cycli we zitten?
            b[-1] = self.Ops.CylinderPressure

            #4. Solve System for Pressure + Update
            p_star = linalg.spsolve(M,b)
            Delta_p = max(p_star, 0) - StateVector[time].Pressure

            ## Update pressure
            StateVector[time].Pressure += self.UnderRelaxP * Delta_p

            # ## Update properties dependent on pressure --> moet dit? Staat niet in algoritme vermeld. Kunnen desnoods eens testen of dit sneller tot convergentie leidt.
            # Density = DensityFunc(StateVector[time])
            # SpecHeat = SpecHeatFunc(StateVector[time])
            # Viscosity = ViscosityFunc(StateVector[time])
            # Conduc = ConducFunc(StateVector[time])
            
            #5. LHS Temperature

            av_u = - StateVector[time].h**2 / (12 * Viscosity) * DDX @ StateVector[time].Pressure + self.Ops.SlidingVelocity / 2
            Uaveraged = av_u        # Om if (k % 500 == 0): ... te laten werken
            # u_plus = np.where(av_u < 0, 0, av_u)
            # u_min = np.where(av_u > 0, 0, av_u)
            # D = u_plus * self.Time.dt @ DDXForward + u_min * self.Time.dt @ DDXBackward
            # E = - Conduc / (Density * SpecHeat) * self.Time.dt @ D2DX2
            # M = np.identity(Grid.Nx) + D + E


            #6. RHS Temperature

            # Q = StateVector[time].h**2 / (12 * Viscosity) * (DDX @ StateVector[time].Pressure)**2 + Viscosity * self.Ops.SlidingVelocity**2 / StateVector[time].h**2
            # RHS = (StateVector[time].Temperature - StateVector[time-1].Temperature) + (self.Time.dt * Q) / (Density * SpecHeat)

            #Boundary conditions

            #7. Solve System for Temperature + Update

            # StateVector[time].Temperature = linalg.spsolve(M, RHS)

            # Density = DensityFunc(StateVector[time])
            # SpecHeat = SpecHeatFunc(StateVector[time])
            # Viscosity = ViscosityFunc(StateVector[time])
            # Conduc = ConducFunc(StateVector[time])
            
            #8. Calculate other quantities: Hydrodynamic load (eq. 37 in assignment), Wall shear stress, Viscous friction force (store all in StateVector)
            #################################################################################################################
            ### Opmerking: moet dit nu hier of moet dit pas in stap 11? --> lijkt mij hier nutteloos om dit te berekenen. ###
            #################################################################################################################

            
            #9. Residuals & Report
            k += 1

            epsP[k] = np.linalg.norm(np.divide(Delta_p, StateVector[time].Pressure)) / self.Grid.Nx


           
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
        # WallShearStress = Viscosity *()         #Uit cursus gehaald, idk of dit correct is...


