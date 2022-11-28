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
# from IOHDF5 import IOHDF5
import time as TimeKeeper
import VisualLib as vis



# """General Settings for Input and Output """
# VisualFeedbackLevel=1 # [0,1,2,3] = [none, per time step, per load iteration, per # reynolds iterations]
# SaveFig2File=False # Save figures to file? True/False
# LoadInitialState=False # Load The InitialState? True/False
# InitTime=0.0 #Initial Time to Load?
# SaveStates=False # Save States to File? True/False

# """I/O Operator"""
# IO=IOHDF5()

""" Input Parameters"""
EngineType='VW 2.0 R4 16v TDI CR 103kW'
OilTemperature=95.0 #C
EngineRPM=2400.0 #rpm
EngineAcceleration=0.0 #;


""" Define Engine Geometry"""
Engine=Engine(EngineType)


"""Define Dry Contact parameters"""
Contact=TriboContact(Engine)


# """1D Computational Grid"""
Nodes=10 #;
Grid=Grid(Contact,Nodes)

# """Temporal Discretization"""
# TimeStep=1e-5 # Choose Temperal Resolution 
# EndTime=4.0*np.pi/(EngineRPM*(2.0*np.pi/60.0))
# Time=Time(EndTime,TimeStep)

# """Define Operational Conditions""" 
# Ops=Ops(Time,Engine,EngineRPM,EngineAcceleration,OilTemperature)


# """Define Two-Phase Lubricant-Vapour flow"""
# Oil=Liquid('SAE5W40')
# Vapour=Gas('SAE5W40')
# Mixture=CavitationModel(Oil,Vapour)


# """Define the State Vector = List of All States over time"""
# StateVector=[]
# for t in range(Time.nt):
#     StateVector.append(State(Grid))


""" Spatial Discretization by Finite Differences """
Discretization=FiniteDifferences(Grid)

# print(Grid.dx)

DDX = Discretization.DDXCentral
# print("DDX =", np.round(DDX.toarray(),2))
D2DX2 = Discretization.D2DX2
# print("D2DX2 =", np.round(D2DX2.toarray(), 2))

x = Grid.x
print("x =", x)

fx = np.sin(x)
dfx_an = np.cos(x)
ddfx_an = - np.sin(x)

dfx = DDX @ fx
print("dfx =", dfx)
print("dfx_an =", dfx_an)
ddfx = D2DX2 @ fx
print("ddxf =", ddfx)
print("ddxf_an =", ddfx_an)