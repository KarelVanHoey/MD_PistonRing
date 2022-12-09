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
from TriboContact import TriboContact #import all classes from file
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
VisualFeedbackLevel=1 # [0,1,2,3] = [none, per time step, per load iteration, per # reynolds iterations]
SaveFig2File=False # Save figures to file? True/False
LoadInitialState=False # Load The InitialState? True/False
InitTime=0.0 #Initial Time to Load?
SaveStates=False # Save States to File? True/False

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
TimeStep=1e-5 # Choose Temperal Resolution 
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
MaxIterReynolds=6000 # originally 5000
TolP=1e-4 
UnderRelaxP=0.001 
TolT=1e-4 
UnderRelaxT=0.01 
Reynolds=ReynoldsSolver(Grid,Time,Ops,Mixture,Discretization)
Reynolds.SetSolver(MaxIterReynolds,TolP,UnderRelaxP,TolT,UnderRelaxT,VisualFeedbackLevel)

""" Set Load Balance loop"""
MaxIterLoad= 40 #originally 40
Tolh0=1e-3 
UnderRelaxh0=0.2
Delta_Load = 0.0

"""Start from Initial guess or Load Initial State"""

"""Start from Scratch: Set State Initial Conditions"""
time=120
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
    

####################################################################################################################
## Test of ReynoldsSolver for a constant film thickness at a given t, without squeeze term or temperature effects ##
####################################################################################################################
eps_h0 = np.ones(MaxIterLoad+1)
Delta_Load = np.zeros(MaxIterLoad)
h0_k = np.zeros(MaxIterLoad + 1)
h0_k[0] = StateVector[time].h0
h0_k[1] = h0_k[0] * 1.01
k_load = 1


while (k_load < MaxIterLoad) and (eps_h0[k_load] > Tolh0): 

    StateVector[time].h0=h0_k[k_load]
    StateVector[time].h = StateVector[time].h0+ 4.0 * Engine.CompressionRing.CrownHeight * (Grid.x**2) / (Engine.CompressionRing.Thickness**2)
    Reynolds.SolveReynolds(StateVector,time)
    
    Delta_Load[k_load] = StateVector[time].HydrodynamicLoad - Ops.CompressionRingLoad[time]
    h0_k[k_load + 1] = max(h0_k[k_load] - UnderRelaxh0 * ((Delta_Load[k_load]) / (Delta_Load[k_load] - Delta_Load[k_load - 1])) * (h0_k[k_load] - h0_k[k_load - 1]), 0.1 * Contact.Roughness)
    k_load += 1 
    
    eps_h0[k_load] = abs(h0_k[k_load] / h0_k[k_load - 1] - 1) 
    
    """Load Balance Output""" 
    print("Load Balance:: Residuals [h0] @Time:",round(Time.t[time]*1000,5),"ms & Iteration:",k_load,"-> [",np.round(eps_h0[k_load],2+int(np.abs(np.log10(Tolh0)))),"]\n")
    if VisualFeedbackLevel>1:
        fig=vis.Report_PT(Grid,StateVector[time])                       
        if SaveFig2File:
            figname="Figures/PT@Time_"+str(round(Time.t[time]*1000,5))+"ms_LoadIteration_"+str(k_load)+".png" 
            fig.savefig(figname, dpi=300)  
        plt.close(fig)

plt.plot(Grid.x,StateVector[time].h)
plt.title("h")
plt.show()

plt.plot(Grid.x,StateVector[time].Pressure)
plt.title("Pressure")
plt.show()    

plt.semilogy(np.arange(0,MaxIterLoad,1),Delta_Load)
plt.title("Delta Load")
plt.show()

plt.plot(np.arange(0,MaxIterLoad+1,1),h0_k)
plt.title("h0")
plt.show()
