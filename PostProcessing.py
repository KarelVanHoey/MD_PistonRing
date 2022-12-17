#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 10:49:17 2020
@ProjectTitle : Thermal Mixed-Hydrodynamic Lubrication Solver for simulation of a compression ring - cylinder liner contact in an IC Engine.
@author: Dieter Fauconnier
@version : 1.0
"""



"""
Import libraries
"""
import numpy as np # Matrix Definitions and operations.
import scipy.integrate as integral # Sparse Matrix Definitions and operations
import matplotlib.pyplot as plt


from EngineParts import Engine  #import all classes from file
from TriboContact import TriboContact #import all classes from file
from Grid import Grid #import all classes from file
from Time import Time #import all classes from file
from Ops import Ops #import all classes from file
from FluidLibrary import Liquid,Gas #import all classes from file
from TwoPhaseModel import CavitationModel 
from SolutionState import State #import all classes from file
from FiniteDifferences import FiniteDifferences
from IOHDF5 import IOHDF5

# from IPython import get_ipython
# get_ipython().magic('reset -sf')


"""I/O Operator"""
IO=IOHDF5()

""" Input Parameters"""
EngineType='VW 2.0 R4 16v TDI CR 103kW'
OilTemperature=95.0 #C
EngineRPM=2400.0 #rpm
EngineAcceleration=0.0

""" Define Engine Geometry"""
Engine=Engine(EngineType)

"""Define Dry Contact parameters"""
Contact=TriboContact(Engine)


"""1D Computational Grid"""
Nodes=256
Grid=Grid(Contact,Nodes)

"""Temporal Discretization"""
TimeStep=5e-5 # Choose Temperal Resolution 
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

""" "Spatial "Discretization by Finite Differences """
Discretization=FiniteDifferences(Grid)



"""Read Data"""
time=0
for time in range(Time.nt-1):
    FileName='Data/Solution_Time_'+str(round(Time.t[time]*1000,4))+'ms.h5' 

    Data=IO.ReadData(FileName)
    StateVector[time].h0=float(Data['State']['h0']) # ok
    StateVector[time].Hersey=float(Data['State']['Hersey']) # ok
    StateVector[time].Lambda=float(Data['State']['Lambda']) # ok
    StateVector[time].HydrodynamicLoad=float(Data['State']['HydrodynamicLoad']) # ok
    StateVector[time].ViscousFriction=float(Data['State']['ViscousFriction'])   # ok
    StateVector[time].AsperityLoad=float(Data['State']['AsperityLoad']) # ok
    StateVector[time].AsperityFriction=float(Data['State']['AsperityFriction']) # ok
    StateVector[time].AsperityContactArea=float(Data['State']['AsperityContactArea']) # ok
    StateVector[time].AsperityContactPressure=float(Data['State']['AsperityContactPressure']) # ok
    StateVector[time].HertzianContactPressure=float(Data['State']['HertzianContactPressure']) # ok
    StateVector[time].COF=float(Data['State']['COF']) # ok 
    # StateVector[time].WearDepthRing=float(Data['State']['WearDepthRing']) # NOT ok?
    
    StateVector[time].h= Data['State']['h'] # ok
    StateVector[time].Pressure=Data['State']['Pressure'] # ok
    StateVector[time].Temperature=Data['State']['Temperature'] # ok
    StateVector[time].WallShearStress=Data['State']['WallShearStress'] # ok
    StateVector[time].WearLocationsCylinder=Data['State']['WearLocationsCylinder'] # ok
    # StateVector[time].WearDepthCylinder=Data['State']['WearDepthCylinder']  # NOT ok?
    
    time+=1
    
 
    

"""Post-Processing"""
#################
##### TO DO #####
#################   