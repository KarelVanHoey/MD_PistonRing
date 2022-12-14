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
for time in [120]:# range(1,Time.nt-1):
    FileName='Data_worn/Time_'+str(round(Time.t[time]*1000,4))+'ms.h5' 

    Data=IO.ReadData(FileName)
    StateVector[time].h0=float(Data['State']['h0']) # ok
    # print("Hersey =", Data['State']['Hersey'])
    StateVector[time].Hersey=Data['State']['Hersey'] # ok
    StateVector[time].Lambda=float(Data['State']['Lambda']) # ok
    StateVector[time].HydrodynamicLoad=float(Data['State']['HydrodynamicLoad']) # ok
    StateVector[time].ViscousFriction=float(Data['State']['ViscousFriction'])   # ok
    StateVector[time].AsperityLoad=float(Data['State']['AsperityLoad']) # ok
    StateVector[time].AsperityFriction=float(Data['State']['AsperityFriction']) # ok
    StateVector[time].AsperityContactArea=float(Data['State']['AsperityContactArea']) # ok
    StateVector[time].AsperityContactPressure=float(Data['State']['AsperityContactPressure']) # ok
    StateVector[time].HertzianContactPressure=float(Data['State']['HertzianContactPressure']) # ok
    StateVector[time].COF=float(Data['State']['COF']) # ok 
    StateVector[time].WearDepthRing=float(Data['State']['WearDepthRing']) # ok
    StateVector[time].Viscosity=Data['State']['Viscosity'] # ok
    
    StateVector[time].h= Data['State']['h'] # ok
    StateVector[time].Pressure=Data['State']['Pressure'] # ok
    StateVector[time].Temperature=Data['State']['Temperature'] # ok
    StateVector[time].WallShearStress=Data['State']['WallShearStress'] # ok
    StateVector[time].WearLocationsCylinder=Data['State']['WearLocationsCylinder'] # ok
    StateVector[time].WearDepthCylinder=Data['State']['WearDepthCylinder']  # ok?
    
    time+=1
    
 
    

"""Post-Processing"""


# interesting_timestamps = [1, 100]
plt.plot(Grid.x, StateVector[120].h)
plt.show()

## Dimensionless film thickness as function of crank angle

# for time in range(Time.nt - 1):
#     if time % 100 == 1:
#         # print("Hersey =", StateVector[time].Hersey)
#         plt.plot(Grid.x, StateVector[time].Hersey)
#         plt.show()

# Lambda_values = np.zeros(Time.nt - 1)
# Hersey_values = np.zeros(Time.nt - 1)
# COF_values = np.zeros(Time.nt-1)

# for time in range(Time.nt - 1):
#     Lambda_values[time] = StateVector[time].Lambda
#     Hersey_values[time] = abs(np.mean(StateVector[time].Hersey))
#     COF_values[time] = abs(StateVector[time].COF)

# plt.plot(Ops.CranckAngle[1:], Lambda_values, 'bo')
# plt.xlabel('Crank angle [rad]')
# plt.ylabel('Dimensionless film thickness [-]')
# plt.show()


## Stribeck curve


# Hersey_values = np.zeros(Time.nt - 1)

# plt.plot(Hersey_values, COF_values, 'bo')
# plt.xlabel('Hersey number [-]')
# plt.ylabel('Coefficient of Friction [-]')
# plt.show()


## Characteristic pressure & temperature fields at interesting and relevant locations

# for time in interesting_timestamps:
#     vis.Report_PT(Grid, StateVector[time], time=time) # plt.show() has to be uncommented in VisualLib


## 2D vectorplot of flow velocity at relevant locations

    ## To be done:  calculate time derivative of h0 for v2 (v1=0)
    ##              implement formulas from slide 17 in Hydrodyn Theory slides

# interesting_timestamps = [100]
# for time in interesting_timestamps:
#     visc_x = StateVector[time].Viscosity
#     density = Mixture.Density(StateVector[time])
#     p = StateVector[time].Pressure
#     p_x = Discretization.DDXCentral @ p
#     DDX = Discretization.DDXCentral
#     p_y = 0.0
#     x_grid = Grid.x
#     z_n = Grid.Nx
#     h = StateVector[time].h
    
#     u1 = 0
#     u2 = Ops.SlidingVelocity[time] # or PistonVelocity ? idk
#     Nx = Grid.Nx
#     u_x = np.zeros((z_n, Nx))
#     u_z = np.zeros((z_n, Nx))

#     z_grid = np.linspace(0.0, h[0], z_n)
    
#     for x in range(Grid.Nx):
#         u_x[:, x] = np.array([1/visc_x[x] * p_x[x] * 1/2 * (z**2 - h[x] * z) if z < h[x] else 0 for z in z_grid])  # Poiseuille
#         u_x[:, x] += np.array([(u2 - u1) / h[x] * z + u1 if z < h[x] else 0 for z in z_grid])                      # Couette
#         # Note: list comprehension is used to get zero for points "inside" of the ring in the vector plot

#         u_z[:, x] = np.trapz(1/density[x] * DDX @ (density[x] * u_x[:, x]),z_grid)*10000

#     ## Reduce amount of datapoint used to generate vectors. Note: plt.quiver wants the same amount of elements in X and Y!
#     skip = 10
#     skip1 = (slice(None, None, skip))
#     skip2 = (slice(None, None, skip), slice(None, None, skip))

#     ## Make grid for vectors
#     X, Z = np.meshgrid(x_grid[skip1], z_grid[skip1])

#     X_l, Z_l = np.meshgrid((x_grid[:Nx//2])[skip1], (z_grid)[skip1])
#     X_r, Z_r = np.meshgrid((x_grid[Nx//2+1:])[skip1], (z_grid)[skip1])


#     ## Make vector plot
#     plt.quiver(X,Z,u_x[skip2],u_z[skip2],minlength=0,scale=350)

#     ### (Un)comment following two line to ensure that no vector crosses the ring.
#     # plt.quiver(X_l,Z_l,(u_x[:,:Nx//2])[skip2],(u_z[:,:Nx//2])[skip2],minlength=0,pivot='tip',scale=350)
#     # plt.quiver(X_r,Z_r,(u_x[:,Nx//2+1:])[skip2],(u_z[:,Nx//2+1:])[skip2],minlength=0,pivot='tail',scale=350)

#     plt.plot(x_grid, StateVector[time].h)
#     plt.show()
            


## Wear of Compression Ring and Wear at Cylinder liner after one combustion cycle

# WearDepthRing_values = np.zeros(Time.nt - 1)
# WearDepthCylinder_values = np.zeros(len(interesting_timestamps))

# for time in range(Time.nt - 1):
#     WearDepthRing_values[time] = StateVector[time]

# plt.plot(Time.t, WearDepthRing_values, 'bo')
# plt.xlabel('Crank angle [rad]')
# plt.ylabel('Dimensionless film thickness [-]')
# plt.show()



# for time in interesting_timestamps:
#     plt.plot(StateVector[time].WearLocationsCylinder, StateVector[time].WearDepthCylinder, 'bo')
#     plt.xlabel('Location on cylinder liner [m]')
#     plt.ylabel('Wear depth [m]')
#     plt.show()



