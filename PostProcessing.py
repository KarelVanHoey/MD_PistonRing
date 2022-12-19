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
import VisualLib as vis

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
    StateVector[time].WearDepthRing=float(Data['State']['WearDepthRing']) # ok?
    
    StateVector[time].h= Data['State']['h'] # ok
    StateVector[time].Pressure=Data['State']['Pressure'] # ok
    StateVector[time].Temperature=Data['State']['Temperature'] # ok
    StateVector[time].WallShearStress=Data['State']['WallShearStress'] # ok
    StateVector[time].WearLocationsCylinder=Data['State']['WearLocationsCylinder'] # ok
    StateVector[time].WearDepthCylinder=Data['State']['WearDepthCylinder']  # ok?
    
    time+=1
    
 
    

"""Post-Processing"""


# Dimensionless film thickness as function of crank angle

Lambda_values = np.zeros(Time.nt - 1)

for time in range(Time.nt - 1):
    Lambda_values[time] = StateVector[time]

plt.plot(Ops.CranckAngle, Lambda_values, 'bo')
plt.xlabel('Crank angle [rad]')
plt.ylabel('Dimensionless film thickness [-]')
plt.show()


# Stribeck curve

COF_values = np.zeros(Time.nt-1)
Hersey_values = np.zeros(Time.nt - 1)

for time in range(Time.nt - 1):
    COF_values[time] = StateVector[time].COF
    Hersey_values[time] = StateVector[time].Hersey

plt.plot(Hersey_values, COF_values, 'bo')
plt.xlabel('Hersey number [-]')
plt.ylabel('Coefficient of Friction [-]')
plt.show()


# Characteristic pressure & temperature fields at interesting and relevant locations

interesting_timestamps = []
for time in interesting_timestamps:
    vis.Report_PT(Grid, StateVector[time], time=time) # plt.show() has to be uncommented in VisualLib


# 2D vectorplot of flow velocity at relevant locations

    ## To be done:  calculate time derivative of h0 for v2 (v1=0)
    ##              implement formulas from slide 17 in Hydrodyn Theory slides

# interesting_timestamps = []
# for time in interesting_timestamps:
#     x_grid = Grid.x
#     z_n = 100
#     z_grid = np.linspace(0.0, StateVector[time].h[0], z_n)
#     u1 = 0
#     u2 = Ops.PistonVelocity[time] # or SlidingVelocity ? idk
#     u_x = np.zeros((Grid.nx, z_n))
#     u_z = np.zeros((Grid.nx, z_n))
#     for x in range(Grid.nx):
#         for z in range(z_n):
            