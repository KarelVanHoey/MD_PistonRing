# -*- coding: utf-8 -*-
"""

Created on Tue Aug 25 17:37:40 2020
@ProjectTitle : Thermal Mixed-Hydrodynamic Lubrication Solver for simulation of a compression ring - cylinder liner contact in an IC Engine.
@author: Dieter Fauconnier
@version : 1.0
"""


"""
Import libraries
"""
import socket
import sys

hostname=socket.gethostname()
Tier2List=['swallot', 'skitty', 'victini', 'slaking', 'kirlia', 'doduo']
if any(list(i in hostname for i in Tier2List)):
    import matplotlib
    matplotlib.use('Agg')
    print('HPC detected: Matplotlib used as backend')

import matplotlib.pyplot as plt
import numpy as np 
import copy as copy

from EngineParts import Engine  #import all classes from file
from TriboContact_worn import TriboContact #import all classes from file
from Grid import Grid #import all classes from file
from Time import Time #import all classes from file
from Ops import Ops #import all classes from file
from FluidLibrary import Liquid,Gas #import all classes from file
from TwoPhaseModel import CavitationModel 
from SolutionState import State #import all classes from file
from FiniteDifferences import FiniteDifferences
from ReynoldsSolver import ReynoldsSolver
from IOHDF5 import IOHDF5
import time as TimeKeeper
import VisualLib as vis


"""General Settings for Input and Output """
VisualFeedbackLevel=0 # [0,1,2,3] = [none, per time step, per load iteration, per # reynolds iterations]
SaveFig2File=False # Save figures to file? True/False
LoadInitialState=False # Load The InitialState? True/False
InitTime=0.0 #Initial Time to Load?
SaveStates=True # Save States to File? True/False

"""I/O Operator"""
IO=IOHDF5()

""" Input Parameters"""
EngineType='VW 2.0 R4 16v TDI CR 103kW'
OilTemperature=95.0 #C
EngineRPM=2400.0 #rpm
EngineAcceleration=0.0 #;


""" Define Engine Geometry"""
Engine=Engine(EngineType)


"""Define Dry Contact parameters"""
Contact=TriboContact(Engine)


"""1D Computational Grid"""
Nodes=256 #;
Grid=Grid(Contact,Nodes)

"""Temporal Discretization"""
TimeStep=5e-5 # Choose Temperal Resolution (original 1e-5)
EndTime=4.0*np.pi/(EngineRPM*(2.0*np.pi/60.0))
Time=Time(EndTime,TimeStep)

"""Define Operational Conditions""" 
Ops=Ops(Time,Engine,EngineRPM,EngineAcceleration,OilTemperature)


"""Define Two-Phase Lubricant-Vapour flow"""
Oil=Liquid('SAE5W40')
Vapour=Gas('SAE5W40')
Mixture=CavitationModel(Oil,Vapour)


"""Define the State Vector = List of All States over time"""
StateVector=[]
for t in range(Time.nt):
    StateVector.append(State(Grid))


""" Spatial Discretization by Finite Differences """
Discretization=FiniteDifferences(Grid)


""" Initialize Reynolds Solver"""
MaxIterReynolds=5000 
TolP=1e-4 
UnderRelaxP=0.001 
TolT=1e-4 
UnderRelaxT=0.01
Reynolds=ReynoldsSolver(Grid,Time,Ops,Mixture,Discretization)
Reynolds.SetSolver(MaxIterReynolds,TolP,UnderRelaxP,TolT,UnderRelaxT,VisualFeedbackLevel)

""" Set Load Balance loop"""
MaxIterLoad= 41 #originally 40
Tolh0=1e-3 
UnderRelaxh0=0.4
# Delta_Load = 0.0

"""Start from Initial guess or Load Initial State"""

time=(np.abs(Time.t - InitTime)).argmin()

if LoadInitialState:
    
    """Start from previous solution: Load Data at t=0"""
    FileName='Data_worn_v3/Time_'+str(round(Time.t[time]*1000,5))+'ms.h5'
    Data=IO.ReadData(FileName)
    StateVector[time].h0=float(Data['State']['h0'])
    StateVector[time].Hersey=Data['State']['Hersey']
    StateVector[time].Lambda=float(Data['State']['Lambda'])
    StateVector[time].HydrodynamicLoad=float(Data['State']['HydrodynamicLoad'])
    StateVector[time].ViscousFriction=float(Data['State']['ViscousFriction'])
    StateVector[time].AsperityLoad=float(Data['State']['AsperityLoad'])
    StateVector[time].AsperityFriction=float(Data['State']['AsperityFriction'])
    StateVector[time].AsperityContactArea=float(Data['State']['AsperityContactArea'])
    StateVector[time].AsperityContactPressure=float(Data['State']['AsperityContactPressure'])
    StateVector[time].HertzianContactPressure=float(Data['State']['HertzianContactPressure'])
    StateVector[time].COF=float(Data['State']['COF'])
    StateVector[time].WearDepthRing=float(Data['State']['WearDepthRing'])
    
    StateVector[time].h= Data['State']['h']
    StateVector[time].Pressure=Data['State']['Pressure']
    StateVector[time].Temperature=Data['State']['Temperature']
    StateVector[time].WallShearStress=Data['State']['WallShearStress']
    StateVector[time].WearLocationsCylinder=Data['State']['WearLocationsCylinder']
    StateVector[time].WearDepthCylinder=Data['State']['WearDepthCylinder']


else:    
    
    """Start from Scratch: Set State Initial Conditions"""
    StateVector[time].Lambda=1.758; 
    StateVector[time].h0=StateVector[time].Lambda*Contact.Roughness
    StateVector[time].h= StateVector[time].h0 + (4.0*Engine.CompressionRing.CrownHeight/Engine.CompressionRing.Thickness**2.0)*Grid.x**2.0
    StateVector[time].Hersey=0.001*np.abs(Ops.PistonVelocity[time])/np.abs(Ops.CompressionRingLoad[time])
    StateVector[time].Pressure=Ops.AtmosphericPressure+0.0*Grid.x
    StateVector[time].Temperature=Ops.OilTemperature+0.0*Grid.x
    StateVector[time].ViscousFriction=0.0
    
    Contact.AsperityContact(StateVector,time)
    StateVector[time].COF=0.0
    StateVector[time].WearDepthRing=0.0
    StateVector[time].WearLocationsCylinder=np.unique(np.round(Ops.PistonPosition,8))
    StateVector[time].WearDepthCylinder=0.0*StateVector[time].WearLocationsCylinder
    
    if SaveStates:
        FileName='Data_worn/Time_'+str(round(Time.t[time]*1000,5))+'ms.h5'
        #Data2File={'Grid': Grid,'Time': Time,'State': StateVector[time]}
        Data2File={'State': StateVector[time]}
        IO.SaveData(FileName,Data2File)

h_template = 4 * Engine.CompressionRing.CrownHeight * (Grid.x**2) / (Engine.CompressionRing.Thickness**2)
h_transl = h_template - 0.25 * Engine.CompressionRing.CrownHeight
h_capped = np.maximum(h_transl, 0)

# print(Ops.SlidingVelocity[time])
# vis.Report_Ops(Time,Ops,time)
# plt.show()
# print(Ops.PistonPosition)
"""Start Time Loop"""
start_time = TimeKeeper.time()
while time<Time.nt:
    
    
    """Initialize State"""
    time += 1
    StateVector[time]=copy.deepcopy(StateVector[time-1])
    print("Time Loop:: Start Calculation @ Time:",round(Time.t[time]*1000,5),"ms \n")

    eps_h0 = np.ones(MaxIterLoad+1)
    Delta_Load = np.ones(MaxIterLoad)
    h0_k = np.zeros(MaxIterLoad + 2)
    h0_k[0] = StateVector[time-1].h0
    h0_k[1] = h0_k[0] * 1.01
    k_load = 1

    # """Test"""
    # StateVector[time].h_worn = np.array([max(.25*Engine.CompressionRing.CrownHeight, i) for i in 4 * Engine.CompressionRing.CrownHeight * (Grid.x**2) / (Engine.CompressionRing.Thickness**2)]) + h0_k[k_load]
    # StateVector[time].h = 4 * Engine.CompressionRing.CrownHeight * (Grid.x**2) / (Engine.CompressionRing.Thickness**2) + h0_k[k_load]

    # plt.plot(StateVector[time].h,label='original')
    # plt.plot(StateVector[time].h_worn,label='worn')
    # plt.legend()
    # plt.show()


    """Start Load Balance Loop"""
    #TODO
    while (k_load < MaxIterLoad) and abs(Delta_Load[k_load-1]) >= 1:# and (eps_h0[k_load] > Tolh0): #
        """a. Calculate Film Thickness Profile"""
        StateVector[time].h = h_capped + h0_k[k_load]
        
        
        """b. Calculate Asperity Load"""
        StateVector[time].Lambda = h0_k[k_load] / Contact.Roughness     # lambda hier gwn berekenen en min naar asperity contact verplaatsen
        Contact.AsperityContact(StateVector,time)
        
        """c. Solve Reynolds""" 
        Reynolds.SolveReynolds(StateVector,time)
        
        """d. Newton Raphson Iteration to find the h0"""

        Delta_Load[k_load] = StateVector[time].HydrodynamicLoad + StateVector[time].AsperityLoad - Ops.CompressionRingLoad[time]
        h0_k[k_load + 1] = max(h0_k[k_load] - UnderRelaxh0 * ((Delta_Load[k_load]) / (Delta_Load[k_load] - Delta_Load[k_load - 1])) * (h0_k[k_load] - h0_k[k_load - 1]), 0.1 * Contact.Roughness)
        
        
        """e. Update & Calculate Residual"""      
        k_load += 1 

        eps_h0[k_load] = abs(h0_k[k_load] / h0_k[k_load - 1] - 1) 
       
        """Load Balance Output""" 
        print("Load Balance:: Residuals [h0] @Time:",round(Time.t[time]*1000,5),"ms & Iteration:",k_load,"-> [",np.round(eps_h0[k_load],2+int(np.abs(np.log10(Tolh0)))),"]")
        print("Load Balance:: Residuals [DeltaLoad] @Time:",round(Time.t[time]*1000,5),"ms & Iteration:",k_load,"-> [",np.round(Delta_Load[k_load - 1],2),"N]\n")

        if VisualFeedbackLevel>1:
           fig=vis.Report_PT(Grid,StateVector[time])                       
           if SaveFig2File:
               figname="Figures/PT@Time_"+str(round(Time.t[time]*1000,5))+"ms_LoadIteration_"+str(k_load)+".png" 
               fig.savefig(figname, dpi=300)  
           plt.close(fig)

        StateVector[time].h0 = h0_k[k_load]
        # StateVector[time].h = StateVector[time].h0 + 4 * Engine.CompressionRing.CrownHeight * (Grid.x**2) / (Engine.CompressionRing.Thickness**2)    
    
    # plt.plot(StateVector[time].h)
    # plt.show()
    # plt.plot(Delta_Load)
    # plt.show()
    
    """Visual Output per time step""" 
    if VisualFeedbackLevel>0:
        vis.Report_Ops(Time,Ops,time)
        fig=vis.Report_PT(Grid,StateVector[time])
        # fig = vis.Report_Ops_PT(Time,Ops,time, Grid,StateVector[time])
        if SaveFig2File:# and round(Time.t[time]*1000,5)*100 % 10 == 0:
            # figname="Figures/PT@Time_"+str(round(Time.t[time]*1000,5))+"ms.png" 
            figname="Figures_worn/PT@Time_"+"{0:.2f}".format(round(Time.t[time]*1000,5))+"ms.png"  
            fig.savefig(figname, dpi=300)
        plt.close(fig)
    
    
       
    
    
    """ Calculate Other Variables of Interest, e.g. COF wear"""
    #TODO
    StateVector[time].Hersey = Ops.SlidingVelocity[time] * StateVector[time].Viscosity / StateVector[time].HydrodynamicLoad
    StateVector[time].COF= (StateVector[time].ViscousFriction + StateVector[time].AsperityFriction) / StateVector[time].HydrodynamicLoad  
    Contact.Wear(Ops,Time,StateVector,time)
 
    
    """Save Output""" 
    if SaveStates:
        FileName='Data_worn/Time_'+str(round(Time.t[time]*1000,5))+'ms.h5'
        Data2File={'State': StateVector[time]}
        IO.SaveData(FileName,Data2File)

    """Close all open Figures"""
    plt.close('all')


print("\n Main Program Completed in %s seconds" % round(TimeKeeper.time() - start_time,0))