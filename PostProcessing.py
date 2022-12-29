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
from matplotlib.patches import Polygon


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
for time in range(1,Time.nt): # [100]:# 
    FileName='Data_v2/Time_'+str(round(Time.t[time]*1000,4))+'ms.h5' 

    Data=IO.ReadData(FileName)
    StateVector[time].h0=float(Data['State']['h0']) # ok
    # print("Hersey =", Data['State']['Hersey'])
    # StateVector[time].Hersey=Data['State']['Hersey'] # ok
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
    StateVector[time].VapourVolumeFraction=Data['State']['VapourVolumeFraction']
    StateVector[time].Density=Data['State']['Density']

    StateVector[time].Hersey = abs(Ops.SlidingVelocity[time]) * StateVector[time].Viscosity / ( StateVector[time].HydrodynamicLoad + StateVector[time].AsperityLoad)

    
    StateVector[time].h= Data['State']['h'] # ok
    StateVector[time].Pressure=Data['State']['Pressure'] # ok
    StateVector[time].Temperature=Data['State']['Temperature'] # ok
    StateVector[time].WallShearStress=Data['State']['WallShearStress'] # ok
    StateVector[time].WearLocationsCylinder=Data['State']['WearLocationsCylinder'] # ok
    StateVector[time].WearDepthCylinder=Data['State']['WearDepthCylinder']  # ok?
    
    time+=1
    
 
    

"""Post-Processing"""
# ## Hydrodynamic load and asperity load
P_hydro = np.zeros(Time.nt - 1)
P_asp = np.zeros(Time.nt - 1)
P_hydro_c = np.zeros(Time.nt - 1)
P_asp_c = np.zeros(Time.nt - 1)
for time in range(Time.nt - 1):
    P_hydro[time] = StateVector[time].HydrodynamicLoad
    P_asp[time] = StateVector[time].AsperityLoad
    # P_hydro_c[time] = StateVector_c[time].HydrodynamicLoad
    # P_asp_c[time] = StateVector_c[time].AsperityLoad
# plt.plot(P_hydro,linestyle='dashdot')
# plt.plot(P_asp,linestyle='dashdot')
plt.plot(Ops.CranckAngle[1:],P_hydro+P_asp,label='Hydrodynamic load + Asperity load')
plt.plot(Ops.CranckAngle[1:],P_hydro,label='Hydrodynamic load')
plt.plot(Ops.CranckAngle[1:],P_asp,label='Asperity load')
plt.plot(Ops.CranckAngle[1:], Ops.CompressionRingLoad[1:],label='Compression ring load')
# plt.plot(Ops.CranckAngle[1:],P_hydro_c+P_asp_c,label='Hydrodynamic load + Asperity load')
pi = np.pi
psi = np.arange(0, 4 * pi + pi/2, step=(pi/2))
plt.xticks(psi,['0','π/2', 'π', '3π/2', '2π','5π/2', '3π', '7π/2', '4π'])
plt.legend()
plt.savefig('PostProcessing/hydrodynamic_and_asp_load.png',dpi=300)
plt.show()
plt.close()

plt.plot((P_hydro+P_asp)/Ops.CompressionRingLoad[1:])
plt.show()
plt.plot(Ops.CompressionRingLoad[1:]-(P_hydro+P_asp))
plt.show()

interesting_timestamps = np.array([1, 94, 500, 563, 999]) #250, 

## Dimensionless film thickness as function of crank angle


Lambda_values = np.zeros(Time.nt - 1)
Hersey_values = np.zeros(Time.nt - 1)
COF_values = np.zeros(Time.nt-1)

for time in range(Time.nt - 1):
    Lambda_values[time] = StateVector[time].Lambda
    Hersey_values[time] = (abs(np.mean(StateVector[time].Hersey)))
    COF_values[time] = abs(StateVector[time].COF)

# plt.plot(Ops.CranckAngle[1:], Lambda_values, 'bo')
# plt.xlabel('Crank angle $\psi$ [rad]')
# plt.ylabel('$\Lambda$ [-]')
# plt.hlines([1, 2.5],-2,15,'k','--', linewidth=.6)
# plt.xlim([-.5, 13])
# plt.vlines(Ops.CranckAngle[time],0,37,'k','--', linewidth=.6)
# plt.ylim([0,37])
# plt.savefig('PostProcessing/Film_thickness.png',dpi=300)
# plt.show()
# plt.close()

## Multicolor filmthickness
gradient = np.linspace(0,1,len(Ops.CranckAngle[1:]))
for i in range(len(Lambda_values)):
    plt.plot(Ops.CranckAngle[i+1],Lambda_values[i],'o',color=(.47,gradient[i],gradient[i]),markersize=3)
plt.xlabel('Crank angle $\psi$ [rad]')
plt.ylabel('$\Lambda$ [-]')
plt.hlines([ 2.5],-2,15,'k',['dashdot'], linewidth=.8,label='Λ = 2.5')
plt.hlines([1],-2,15,'k',['dotted'], linewidth=.8,label='Λ = 1')
plt.xlim([-.5, 13])
plt.vlines(Ops.CranckAngle[interesting_timestamps],-3,37,'k','--', linewidth=.6)
plt.ylim([-1,37])
pi = np.pi
psi = np.arange(0, 4 * pi + pi/2, step=(pi/2))
plt.xticks(psi,['0','π/2', 'π', '3π/2', '2π','5π/2', '3π', '7π/2', '4π'])
plt.legend()
plt.savefig('PostProcessing/Film_thickness.png',dpi=300)
plt.show()
plt.close()


# # ## Stribeck curve



# plt.plot(Hersey_values, COF_values, 'bo')
# plt.xlabel('Hersey number [-]')
# plt.ylabel('Coefficient of Friction [-]')
# plt.savefig('PostProcessing/Stribeck_curve.png',dpi=300)
# # plt.show()
# plt.close()

# ## Multi color Stribeck
# for i in range(len(COF_values)):
#     plt.plot(Hersey_values[i]*10**4,COF_values[i],'o',color=(.47,gradient[i],gradient[i]),markersize=3)
# plt.xlabel('Hersey number x$ 10^4$ [-]')
# plt.ylabel('Coefficient of Friction [-]')
# # plt.xscale('log')
# plt.savefig('PostProcessing/Stribeck_curve.png',dpi=300)
# plt.show()
# plt.close()



# ## Interesting points
# vis.Report_Ops(Time,Ops,interesting_timestamps)
# plt.savefig('PostProcessing/Interesting_points.png',dpi=300)
# plt.show()
# plt.close()




# ## Characteristic pressure & temperature fields at interesting and relevant locations

# for time in interesting_timestamps:
    
#     vis.Report_PT(Grid, StateVector[time], time=time) # plt.show() has to be uncommented in VisualLib --> kan ook gwn hier
#     figname="PostProcessing/PT@Time_"+"{0:.2f}".format(round(Time.t[time]*1000,5))+"ms.png" 
#     plt.savefig(figname,dpi=300)
#     plt.show()
#     plt.close()


# ### Vapour volume fraction, viscosity, Density at relavant locations
# # interesting_timestamps = [100]
# for time in interesting_timestamps:
#     plt.plot(Grid.x*1000,StateVector[time].VapourVolumeFraction)
#     plt.ylabel(chr(945) + ' [-]'+' at time =' + str(time*5/100) + 'ms')
#     plt.xlabel('x [mm]')
#     figname="PostProcessing/alpha@Time_"+"{0:.2f}".format(round(Time.t[time]*1000,5))+"ms.png" 
#     plt.savefig(figname,dpi=300)
#     # plt.show()
#     plt.close()


#     plt.plot(Grid.x*1000,StateVector[time].Density)
#     plt.ylabel('Density [kg/m³]'+' at time =' + str(time*5/100) + 'ms')
#     plt.xlabel('x [mm]')
#     figname="PostProcessing/rho@Time_"+"{0:.2f}".format(round(Time.t[time]*1000,5))+"ms.png" 
#     plt.savefig(figname,dpi=300)
#     # plt.show()
#     plt.close()


#     plt.plot(Grid.x*1000,StateVector[time].Viscosity)
#     plt.ylabel('Viscosity [Pa s]'+' at time =' + str(time*5/100) + 'ms')
#     plt.xlabel('x [mm]')
#     figname="PostProcessing/µ@Time_"+"{0:.2f}".format(round(Time.t[time]*1000,5))+"ms.png" 
#     plt.savefig(figname,dpi=300)
#     # plt.show()
#     plt.close()

    




## 2D vectorplot of flow velocity at relevant locations

    

# interesting_timestamps = np.array([1, 500, 999]) #250, 

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
#         # print(u_x[:,x])
#         # print(np.array([1/visc_x[x] * p_x[x] * 1/2 * (z**2 - h[x] * z) if z < h[x] else 0 for z in z_grid]))
#         # u_z[:, x] = np.trapz(1/density[x] * DDX @ (density[x] * u_x[:, x]),z_grid)*10000
    
#     # np.average(u_x)
#     ## Reduce amount of datapoint used to generate vectors. Note: plt.quiver wants the same amount of elements in X and Y!
#     skip = 10
#     scale = 1
#     if time in np.array([1, 500, 999]):
#         skip = 12
#         scale = 40
#     skip1 = (slice(None, None, skip))
#     skip2 = (slice(None, None, skip), slice(None, None, skip))

#     ## Make grid for vectors
#     X, Z = np.meshgrid(x_grid[skip1], z_grid[skip1])

#     X_l, Z_l = np.meshgrid((x_grid[:Nx//2])[skip1], (z_grid)[skip1])
#     X_r, Z_r = np.meshgrid((x_grid[Nx//2+1:])[skip1], (z_grid)[skip1])


#     ## Make vector plot
#     pts = [[-0.75,0.02]]
#     for i in range(len(x_grid)):
#         pts.append([x_grid[i]*1000,StateVector[time].h[i]*1000])
#     pts.append([0.75,0.02])
#     p = Polygon(pts,closed=False,ec='black',fc='grey',zorder=.1)
#     plt.gca().add_patch(p)

#     plt.quiver(X*1000,Z*1000,u_x[skip2]*scale,u_z[skip2],pivot='tail',minlength=0,scale=350) #scale=350 for v!=0

#     ### (Un)comment following two line to ensure that no vector crosses the ring.
#     # plt.quiver(X_l,Z_l,(u_x[:,:Nx//2])[skip2],(u_z[:,:Nx//2])[skip2],minlength=0,pivot='tip',scale=350)
#     # plt.quiver(X_r,Z_r,(u_x[:,Nx//2+1:])[skip2],(u_z[:,Nx//2+1:])[skip2],minlength=0,pivot='tail',scale=350)

#     # plt.plot(x_grid*1000, StateVector[time].h*1000)
#     figname="PostProcessing/Vectorplot@Time_"+"{0:.2f}".format(round(Time.t[time]*1000,5))+"ms.png" 
#     plt.xlabel('x [mm]')
#     plt.ylabel('z [mm]')
#     plt.title('Vectorplot for t =' + str(time*5/100) + 'ms')
#     plt.xlim([-.8,.8])
#     plt.ylim([0.0, 0.0175])
#     plt.tight_layout()
#     # plt.savefig(figname,dpi=300)
#     # plt.show()
#     plt.close()

            


## Wear of Compression Ring and Wear at Cylinder liner after one combustion cycle

WearDepthRing_values = np.zeros(Time.nt - 1)
WearDepthCylinder_values = np.zeros(len(interesting_timestamps))

for time in range(Time.nt - 1):
    WearDepthRing_values[time] = StateVector[time].WearDepthRing

plt.plot(Ops.CranckAngle[1:], WearDepthRing_values, 'bo')
plt.xlabel('Crank angle [rad]')
# plt.plot(Time.t[1:], WearDepthRing_values, 'bo')
# plt.xlabel('time [s]')
plt.ylabel('Weardepth ring [mm]')
pi = np.pi
psi = np.arange(0, 4 * pi + pi/2, step=(pi/2))
plt.xticks(psi,['0','π/2', 'π', '3π/2', '2π','5π/2', '3π', '7π/2', '4π'])
plt.savefig('PostProcessing/WearDepth_ring.png',dpi=300)
# plt.show()
plt.close()


# for time in interesting_timestamps:
time = 999 # We are only interested in wear after a full combustion cycle
plt.plot(StateVector[time].WearLocationsCylinder*1000 - 95.5, StateVector[time].WearDepthCylinder, 'bo')
plt.xlabel('Location on cylinder liner [mm]')
plt.ylabel('Wear depth [m]')
plt.savefig('PostProcessing/Wear_cylinder.png',dpi=300)    
# plt.show()
plt.close()

print('Maximim wear depth on cylinder sleeve = ' + str(max(StateVector[time].WearDepthCylinder)))


# # lifetime compression ring

# WearDepth_one_comb_cycle = StateVector[999].WearDepthRing # constant wear rate assumed
# reduction = 0.2 * Engine.CompressionRing.CrownHeight
# nr_comb_cycles = reduction / WearDepth_one_comb_cycle
# rot = nr_comb_cycles * 2 #  1 combustion cycle = 2 rotations of crackshaft
# km = rot / 1200 # 120 km/h @ 2400 rpm --> 1 km/30s --> 1km = 1200 rot
# print('aantal km= ', km)

# lifetime cylinder liner

max_one_comb_cycle = np.max(StateVector[999].WearDepthCylinder)
nr_comb_cycles_2 = 0.000002 / max_one_comb_cycle
rot2 = nr_comb_cycles_2 * 2
km2 = rot2 / 1200
print('aantal km cylinder liner = ', km2)
