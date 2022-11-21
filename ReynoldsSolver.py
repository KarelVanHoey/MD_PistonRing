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
            Density_prev = DensityFunc(StateVector[time-1]) #Needed to calculate time differentiation in step 2. --> or use PreviousDensity()?
            SpecHeat = SpecHeatFunc(StateVector[time])
            Viscosity = ViscosityFunc(StateVector[time])
            Conduc = ConducFunc(StateVector[time])

            phi = np.divide(np.multiply(Density, StateVector[time].h**3), 12 * Viscosity)
            ### OLD: phi_sparse = sparse.csr_matrix(phi) #Nodig om sparse matrix te krijgen als resultaat van phi.multiply(sparse)
            ### Karel: Ik vind deze methode niet zo overzichtelijk. Ik zou van phi_sparse een diagonaalmatrix maken zoals in de opgave
            ### Karel: phi_sparse = sparse.diags(phi)
            ### Karel: phi_column = sparse.csc_matrix(np.matrix(phi).T)
            ### Victor: akkoord, had dat gemist in de opgave.
            phi_diag = sparse.diags(phi)
            phi_column = sparse.csc_matrix(np.matrix(phi).T)
        
            #1. LHS Pressure

            ### OLD: M = phi_sparse.multiply(D2DX2) + DDX @ phi_sparse.multiply(DDX) #Checken of dit ok is!!
            ### Karel: suggestie: sparse.diags(phi) @ D2DX2 + sparse.diags(DDX @ phi_column) @ DDX
            ### Victor: Akkoord --> geeft nu wel foutmelding: hij wilt DDX @ phi_column niet in sparse.diags, wss omdat hij dit als sparse column ziet + diags heeft een rij nodig...
            
            #################################################################################################################################
            # Oere lelijk, maar enige dat gelijk werkt. sparse.diags aanvaard enkel niet-sparse rij array. Dus via .T gaat sparse kolom     #
            # naar sparse rij. .toarray() maakt er dan een normale array van, die om een of andere reden genest is: [[data1, data2,...]]    #
            # de [0] haalt er dan de juiste array uit...                                                                                    #
            #################################################################################################################################
            M = phi_diag @ D2DX2 + sparse.diags((DDX @ phi_column).T.toarray()[0]) @ DDX
            
        
            #2. RHS Pressure
            ### OLD: b = self.Ops.SlidingVelocity[time]/2 * DDX * np.multiply(Density, StateVector[time].h) + (np.multiply(Density, StateVector[time].h) - np.multiply(Density_prev, StateVector[time-1].h)) / self.Time.dt 
                #Note: Squeeze term: backward time differentiation for d(rho*h)/dt!
                ### Karel: Ik denk dat deze niet correct is.
                ### Karel: suggestie: U/2 * DDX @ sparse.csc_matrix(np.matrix(np.multiply(Density, StateVector[time].h).T)) + "kolomvector van die hele zooi hierboven"
                ### Victor: Akkoord
            b = self.Ops.SlidingVelocity[time]/2 * DDX @ sparse.csc_matrix(np.matrix(np.multiply(Density, StateVector[time].h)).T) + sparse.csc_matrix(np.matrix((np.multiply(Density, StateVector[time].h) - np.multiply(Density_prev, StateVector[time-1].h)) / self.Time.dt).T)

   
            #3. Set Boundary Conditions Pressure --> NOTE work with absolute pressure!!

            SetDirichletLeft(M) # OLD: --> geeft om een of andere reden foutmelding ... "TypeError: memoryview: invalid slice key"
            # Victor: met aanpassing Karel kunnen we deze functies toch gebruiken!
            # OLD: M.data[0] = 1
            # OLD: M.data[1] = 0
            # OLD: M.data[2] = 0
    
            SetDirichletRight(M)
            # OLD: M.data[-1] = 1
            # OLD: M.data[-2] = 0
            # OLD: M.data[-3] = 0 
            # print(M)

            b[0] = self.Ops.AtmosphericPressure     
            b[-1] = self.Ops.CylinderPressure[time] # Victor: [psi] --> Moet hier gedefinieerd worden waar in de cycli we zitten?
            
            #4. Solve System for Pressure + Update
            p_star = linalg.spsolve(M,b)
            ### OLD: Delta_p = np.zeros(p_star.shape)
            ### for i in range(len(Delta_p)):
            ###    Delta_p[i] = max(p_star[i], 0) - StateVector[time].Pressure[i]
            ### Karel: Suggestie om if-loop te vermijden en snelheid te winnen:
            ### Karel: Delta_p = np.maximum(p_star, np.zeros(len(Delta_p))) - StateVector[time].Pressure
            ### Victor: Akkoord
            Delta_p = np.maximum(p_star, np.zeros(len(p_star))) - StateVector[time].Pressure

            ## Update pressure
            StateVector[time].Pressure += self.UnderRelaxP * Delta_p

            # ## Victor: Update properties dependent on pressure --> moet dit? Staat niet in algoritme vermeld. Kunnen desnoods eens testen of dit sneller tot convergentie leidt.
            # Density = DensityFunc(StateVector[time])
            # SpecHeat = SpecHeatFunc(StateVector[time])
            # Viscosity = ViscosityFunc(StateVector[time])
            # Conduc = ConducFunc(StateVector[time])
            
            #5. LHS Temperature
            print(sparse.csr_matrix(StateVector[time].Pressure))
            print(DDX)
            av_u = - np.multiply(np.divide(StateVector[time].h**2 , (12 * Viscosity)), (DDX @ sparse.csr_matrix(StateVector[time].Pressure))) + self.Ops.SlidingVelocity / 2 #hier zit nog een fout: hoe dp/dx krijgen?
            ### Karel: Ik denk dat hier best alles omgezet wordt in kolomvectoren, nu zijn het gewone arrays (1 rij). Dan gaat DDX correct werken.
            Uaveraged = av_u        # Om if (k % 500 == 0): ... te laten werken
            print(av_u)
            u_plus = np.where(av_u < 0, 0, av_u)
            ### Hier kan ook np.maximum gebruikt worden zoals in lijn 128
            print(u_plus)
            u_min = np.where(av_u > 0, 0, av_u)
            print(u_min)
            u_plus_sparse = sparse.csr_matrix(u_plus) * self.Time.dt
            print('1', u_plus_sparse)
            u_min_sparse = sparse.csr_matrix(u_min) * self.Time.dt
            print('2', u_min_sparse)
            D = u_plus_sparse.multiply(DDXForward)  + u_min_sparse.multiply(DDXBackward)
            E1 = - Conduc / (Density * SpecHeat) * self.Time.dt
            E = sparse.csr_matrix(E1) @ D2DX2
            M1 = np.identity(self.Grid.Nx) + D + E


            #6. RHS Temperature

            Q = StateVector[time].h**2 / (12 * Viscosity) * (DDX @ StateVector[time].Pressure)**2 + Viscosity * self.Ops.SlidingVelocity**2 / StateVector[time].h**2
            RHS = (StateVector[time].Temperature - StateVector[time-1].Temperature) + (self.Time.dt * Q) / (Density * SpecHeat)

            #Boundary conditions
            if self.Ops.SlidingVelocity <= 0:
                M1[0,0:1] = [-1/self.Grid.Nx, 1/self.Grid.Nx] ### Karel: hier moet het denk ik Engine.CompressionRing.Thickness/self.Grid.Nx zijn.
                M1[-1, -1] = 1
                M1[0,3:-1] = 0   ### Karel: niet juist denk ik, alle getallen vanaf 3 in die rij moeten nul zijn --> [0, 3:] (zie voorbeeldje in Test_Sparse_Matrix.py)
                M1[-1,1:-2] = 0  ### Karel: hier moet het [-1, 1:-1] zijn denk ik
                RHS[0] = 0
                RHS[-1] = self.Ops.OilTemperature
            else:
                M1[0,0] = 1
                M1[2:-1, 0] = 0
                M1[-1,-2:-1] = [-1/self.Grid.Nx, 1/self.Grid.Nx]
                M1[-1,1:-3] = 0
                RHS[0] = self.Ops.OilTemperature
                RHS[-1] = 0
            #7. Solve System for Temperature + Update

            T_star = linalg.spsolve(M1, RHS)
            delta_T = T_star - StateVector[time].Temperature
            StateVector[time].Temperature += delta_T * self.UnderRelaxT

            # Density = DensityFunc(StateVector[time])
            # SpecHeat = SpecHeatFunc(StateVector[time])
            # Viscosity = ViscosityFunc(StateVector[time])
            # Conduc = ConducFunc(StateVector[time])
            
            #8. Calculate other quantities: Hydrodynamic load (eq. 37 in assignment), Wall shear stress, Viscous friction force (store all in StateVector)
            #################################################################################################################
            ### Opmerking: moet dit nu hier of moet dit pas in stap 11? --> lijkt mij nutteloos om dit hier te berekenen. ###
            #################################################################################################################

            
            #9. Residuals & Report
            k += 1

            epsP[k] = np.linalg.norm(np.divide(Delta_p, StateVector[time].Pressure)) / self.Grid.Nx
            ### Victor: @Niels: ook een residual voor de temp. maken. Voorstel:
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
        # WallShearStress = Viscosity *()         #Uit cursus gehaald, idk of dit correct is...


