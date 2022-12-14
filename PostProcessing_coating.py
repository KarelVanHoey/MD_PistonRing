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
from matplotlib.patches import Polygon, Rectangle


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

StateVector_c=[]
for t in range(Time.nt):
    StateVector_c.append(State(Grid))

""" "Spatial "Discretization by Finite Differences """
Discretization=FiniteDifferences(Grid)



"""Read Data"""
time=0
for time in range(1,Time.nt): # [100]:# 
    FileName = 'Data_v3/Time_'+str(round(Time.t[time]*1000,4))+'ms.h5' 

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

time=0
for time in range(1,Time.nt): # [100]:# 
    FileName = 'Data_coated_v2/Time_'+str(round(Time.t[time]*1000,4))+'ms.h5' 

    Data=IO.ReadData(FileName)
    
    StateVector_c[time].h0=float(Data['State']['h0']) # ok
    # print("Hersey =", Data['State']['Hersey'])
    # StateVector[time].Hersey=Data['State']['Hersey'] # ok
    StateVector_c[time].Lambda=float(Data['State']['Lambda']) # ok
    StateVector_c[time].HydrodynamicLoad=float(Data['State']['HydrodynamicLoad']) # ok
    StateVector_c[time].ViscousFriction=float(Data['State']['ViscousFriction'])   # ok
    StateVector_c[time].AsperityLoad=float(Data['State']['AsperityLoad']) # ok
    StateVector_c[time].AsperityFriction=float(Data['State']['AsperityFriction']) # ok
    StateVector_c[time].AsperityContactArea=float(Data['State']['AsperityContactArea']) # ok
    StateVector_c[time].AsperityContactPressure=float(Data['State']['AsperityContactPressure']) # ok
    StateVector_c[time].HertzianContactPressure=float(Data['State']['HertzianContactPressure']) # ok
    StateVector_c[time].COF=float(Data['State']['COF']) # ok 
    StateVector_c[time].WearDepthRing=float(Data['State']['WearDepthRing']) # ok
    StateVector_c[time].Viscosity=Data['State']['Viscosity'] # ok
    StateVector_c[time].VapourVolumeFraction=Data['State']['VapourVolumeFraction']
    StateVector_c[time].Density=Data['State']['Density']

    StateVector_c[time].Hersey = abs(Ops.SlidingVelocity[time]) * StateVector[time].Viscosity / ( StateVector[time].HydrodynamicLoad + StateVector[time].AsperityLoad)

    
    StateVector_c[time].h= Data['State']['h'] # ok
    StateVector_c[time].Pressure=Data['State']['Pressure'] # ok
    StateVector_c[time].Temperature=Data['State']['Temperature'] # ok
    StateVector_c[time].WallShearStress=Data['State']['WallShearStress'] # ok
    StateVector_c[time].WearLocationsCylinder=Data['State']['WearLocationsCylinder'] # ok
    StateVector_c[time].WearDepthCylinder=Data['State']['WearDepthCylinder']  # ok?

    time+=1


    
 

"""Post-Processing"""

## Hydrodynamic load and asperity load
P_hydro = np.zeros(Time.nt - 1)
P_asp = np.zeros(Time.nt - 1)
P_hydro_c = np.zeros(Time.nt - 1)
P_asp_c = np.zeros(Time.nt - 1)
for time in range(Time.nt - 1):
    P_hydro[time] = StateVector[time].HydrodynamicLoad
    P_asp[time] = StateVector[time].AsperityLoad
    P_hydro_c[time] = StateVector_c[time].HydrodynamicLoad
    P_asp_c[time] = StateVector_c[time].AsperityLoad
# plt.plot(P_hydro,linestyle='dashdot')
# plt.plot(P_asp,linestyle='dashdot')
plt.plot(Ops.CranckAngle[1:],Ops.CompressionRingLoad[1:],linestyle='dashdot',label='Compression ring load')
plt.plot(Ops.CranckAngle[1:],P_hydro_c,label='Hydrodynamic load')
plt.plot(Ops.CranckAngle[1:],P_asp_c,label='Asperity load')
plt.plot(Ops.CranckAngle[1:],P_hydro_c+P_asp_c,label='Hydrodynamic load + Asperity load')
pi = np.pi
psi = np.arange(0, 4 * pi + pi/2, step=(pi/2))
plt.xticks(psi,['0','??/2', '??', '3??/2', '2??','5??/2', '3??', '7??/2', '4??'])
plt.legend()
plt.savefig('PostProcessing_coating/COATING_hydrodynamic_and_asp_load.png',dpi=300)
# plt.show()
plt.close()

interesting_timestamps = np.array([1, 94, 500, 563, 999]) #250, 

## Dimensionless film thickness as function of crank angle


Lambda_values = np.zeros(Time.nt - 1)
Lambda_c = np.zeros(Time.nt - 1)
Hersey_values = np.zeros(Time.nt - 1)
Hersey_c = np.zeros(Time.nt - 1)
COF_values = np.zeros(Time.nt-1)
COF_c = np.zeros(Time.nt-1)
h0 = np.zeros(Time.nt-1)
h0_c = np.zeros(Time.nt-1)

for time in range(Time.nt - 1):
    Lambda_values[time] = StateVector[time].Lambda
    Lambda_c[time] = StateVector_c[time].Lambda
    Hersey_values[time] = (abs(np.mean(StateVector[time].Hersey)))
    COF_values[time] = abs(StateVector[time].COF)
    Hersey_c[time] = (abs(np.mean(StateVector_c[time].Hersey)))
    COF_c[time] = abs(StateVector_c[time].COF)
    h0[time] = StateVector[time].h0
    h0_c[time] = StateVector_c[time].h0


## Single color film thickness
plt.plot(Ops.CranckAngle[1:], Lambda_values, 'o',label='Normal roughness',markersize=3)
plt.plot(Ops.CranckAngle[1:], Lambda_c, 'o',label='Halved roughness',markersize=3)
plt.xlabel('Crank angle $\psi$ [rad]')
plt.ylabel('$\Lambda$ [-]')
plt.hlines([ 2.5],-2,15,'k',['dashdot'], linewidth=.8,label='?? = 2.5')
plt.hlines([1],-2,15,'k',['dotted'], linewidth=.8,label='?? = 1')
plt.xlim([-.5, 13])
# plt.vlines(Ops.CranckAngle[time],-1,60,'k','--', linewidth=.6)
# plt.ylim([-1,60])
plt.legend(loc= (.58,.75))
pi = np.pi
psi = np.arange(0, 4 * pi + pi/2, step=(pi/2))
plt.xticks(psi,['0','??/2', '??', '3??/2', '2??','5??/2', '3??', '7??/2', '4??'])

# plt.vlines(Ops.CranckAngle[interesting_timestamps],-3,60,'k','--', linewidth=.6)

# p = Rectangle((2.156*pi,-3),1.364*pi,70,ec='red',fc='white',zorder=.1,hatch='/') #,fc='red'
# plt.gca().add_patch(p)

plt.savefig('PostProcessing_coating/COATING_Film_thickness_comparison.png',dpi=300)
# plt.show()
plt.close()

## Multicolor filmthickness
gradient = np.linspace(0,1,len(Ops.CranckAngle[1:]))
for i in range(len(Lambda_values)):
    plt.plot(Ops.CranckAngle[i+1],Lambda_c[i],'o',color=(.47,gradient[i],gradient[i]),markersize=3)
plt.xlabel('Crank angle $\psi$ [rad]')
plt.ylabel('$\Lambda$ [-]')
plt.hlines([ 2.5],-2,15,'k',['dashdot'], linewidth=.8,label='?? = 2.5')
plt.hlines([1],-2,15,'k',['dotted'], linewidth=.8,label='?? = 1')
plt.xlim([-.5, 13])
# plt.vlines(Ops.CranckAngle[interesting_timestamps],-3,60,'k','--', linewidth=.6)
# plt.ylim([-1,60])
# p = Rectangle((2.156*pi,-3),1.364*pi,70,ec='red',fc='white',zorder=.1,hatch='/') #,fc='red'
# plt.gca().add_patch(p)
pi = np.pi
psi = np.arange(0, 4 * pi + pi/2, step=(pi/2))
plt.xticks(psi,['0','??/2', '??', '3??/2', '2??','5??/2', '3??', '7??/2', '4??'])
plt.savefig('PostProcessing_coating/COATING_Film_thickness.png',dpi=300)
# plt.show()
plt.close()

## Dimensionfull film thickness!
plt.plot(Ops.CranckAngle[1:], h0, 'o',label='Normal roughness',markersize=3)
plt.plot(Ops.CranckAngle[1:], h0_c, 'o',label='Halved roughness',markersize=3)
plt.xlabel('Crank angle $\psi$ [rad]')
plt.ylabel('h0 [m]')
# plt.hlines([ 2.5],-2,15,'k',['dashdot'], linewidth=.8,label='?? = 2.5')
# plt.hlines([1],-2,15,'k',['dotted'], linewidth=.8,label='?? = 1')
plt.xlim([-.5, 13])
# plt.vlines(Ops.CranckAngle[time],-1,60,'k','--', linewidth=.6)
# plt.ylim([-1,60])
plt.legend(loc= (.58,.75))
pi = np.pi
psi = np.arange(0, 4 * pi + pi/2, step=(pi/2))
plt.xticks(psi,['0','??/2', '??', '3??/2', '2??','5??/2', '3??', '7??/2', '4??'])

# plt.vlines(Ops.CranckAngle[interesting_timestamps],-3,60,'k','--', linewidth=.6)

# p = Rectangle((2.156*pi,-3),1.364*pi,70,ec='red',fc='white',zorder=.1,hatch='/') #,fc='red'
# plt.gca().add_patch(p)

plt.savefig('PostProcessing_coating/COATING_h0_comparison.png',dpi=300)
# plt.show()
plt.close()

## Stribeck curve

plt.plot(Hersey_values*10**4, COF_values, 'o',markersize=3,label='Original roughness')
plt.plot(Hersey_c*10**4, COF_c, 'o',markersize=3,label='Halved roughness')
plt.xlabel('Hersey number x$ 10^4$ [-]')
plt.ylabel('Coefficient of Friction [-]')
plt.legend()
plt.savefig('PostProcessing_coating/COATING_Stribeck_curve_single_color.png',dpi=300)
# plt.show()
plt.close()

## Multi color Stribeck
for i in range(len(COF_values)):
    plt.plot(Hersey_values[i]*10**4,COF_values[i],'o',color=(.47,gradient[i],gradient[i]),markersize=3)
plt.xlabel('Hersey number x$ 10^4$ [-]')
plt.ylabel('Coefficient of Friction [-]')
# plt.xscale('log')
plt.savefig('PostProcessing_coating/COATING_Stribeck_curve.png',dpi=300)
# plt.show()
plt.close()



# ## Interesting points
# vis.Report_Ops(Time,Ops,interesting_timestamps)
# plt.savefig('PostProcessing_worn/WORN_Interesting_points.png',dpi=300)
# plt.show()
# plt.close()




# ## Characteristic pressure & temperature fields at interesting and relevant locations
i=0
for time in interesting_timestamps:
    i+=1
    vis.Report_PT(Grid, StateVector_c[time], time=time) # plt.show() has to be uncommented in VisualLib --> kan ook gwn hier
    figname="PostProcessing_coating/COATING_PT@Time_"+"{0:.2f}".format(round(Time.t[time]*1000,5))+"ms.png" 
    plt.title('Point ' +str(i))
    plt.tight_layout()
    plt.savefig(figname,dpi=300)
    # plt.show()
    plt.close()


## Wear of Compression Ring and Wear at Cylinder liner after one combustion cycle

WearDepthRing_values = np.zeros(Time.nt - 1)
WearDepthRing_c = np.zeros(Time.nt - 1)

for time in range(Time.nt - 1):
    WearDepthRing_values[time] = StateVector[time].WearDepthRing
    WearDepthRing_c[time] = StateVector_c[time].WearDepthRing

plt.plot(Ops.CranckAngle[1:], WearDepthRing_values, 'o',label='Normal roughness',markersize=3)
plt.plot(Ops.CranckAngle[1:], WearDepthRing_c , 'o',label='Halved roughness',markersize=3)
plt.xlabel('Crank angle [rad]')
# plt.plot(Time.t[1:], WearDepthRing_values, 'bo')
# plt.xlabel('time [s]')
pi = np.pi
psi = np.arange(0, 4 * pi + pi/2, step=(pi/2))
plt.xticks(psi,['0','??/2', '??', '3??/2', '2??','5??/2', '3??', '7??/2', '4??'])
plt.ylabel('Weardepth ring [mm]')
# plt.ylim([-.1e-13, 2.5e-13])
plt.legend()
plt.savefig('PostProcessing_coating/COATING_WearDepth_ring_comparison.png',dpi=300)
# plt.show()
plt.close()


# for time in interesting_timestamps:
time = 999 # We are only interested in wear after a full combustion cycle
plt.plot(StateVector[time].WearLocationsCylinder*1000 - 95.5, StateVector[time].WearDepthCylinder, 'o',label='Normal roughness',markersize=3)
plt.plot(StateVector_c[time].WearLocationsCylinder*1000 - 95.5, StateVector_c[time].WearDepthCylinder, 'o',label='Halved roughness',markersize=3)
plt.xlabel('Location on cylinder liner [mm]')
plt.ylabel('Wear depth [m]')
plt.legend()
plt.savefig('PostProcessing_coating/COATING_Wear_cylinder_comparison.png',dpi=300)    
# plt.show()
plt.close()

print('Maximum wear depth on cylinder sleeve = ' + str(max(StateVector[time].WearDepthCylinder)))


## lifetime compression ring

WearDepth_one_comb_cycle = StateVector_c[999].WearDepthRing #- StateVector_c[1].WearDepthRing # constant wear rate assumed
reduction = 0.2 * Engine.CompressionRing.CrownHeight
# reduction = 20e-6 # coating is 20??m thick --> other level is smaller!! 0.2*10??m
nr_comb_cycles = reduction / WearDepth_one_comb_cycle
rot = nr_comb_cycles * 2 #  1 combustion cycle = 2 rotations of crackshaft
km = rot / 1200 # 120 km/h @ 2400 rpm --> 1 km/30s --> 1km = 1200 rot
print('aantal km= ', km)

## lifetime cylinder liner

max_one_comb_cycle = np.max(StateVector_c[999].WearDepthCylinder)
nr_comb_cycles_2 = 0.000002 / max_one_comb_cycle
rot2 = nr_comb_cycles_2 * 2
km2 = rot2 / 1200
print('aantal km cylinder liner = ', km2)


## Quantification of efficiency

F_tan = np.zeros(Time.nt - 1)
F_tan_c = np.zeros(Time.nt - 1)
v = np.zeros(Time.nt - 1)
for time in range(Time.nt - 1):
    v[time] = Ops.SlidingVelocity[time]
    F_tan[time] = StateVector[time].ViscousFriction + StateVector[time].AsperityFriction
    F_tan_c[time] = StateVector_c[time].ViscousFriction + StateVector_c[time].AsperityFriction

W_dot_shear = abs(F_tan * v)
W_dot_shear_c = abs(F_tan_c * v)

### Plotting of shear power

plt.plot(Ops.CranckAngle[1:], W_dot_shear, label='Normal roughness')
plt.plot(Ops.CranckAngle[1:], W_dot_shear_c, label='Halved roughness')
plt.ylabel('Shear power [W/m]')
plt.xlabel('Crank angle $\psi$ [rad]')
pi = np.pi
psi = np.arange(0, 4 * pi + pi/2, step=(pi/2))
plt.xticks(psi,['0','??/2', '??', '3??/2', '2??','5??/2', '3??', '7??/2', '4??'])
plt.legend()
plt.savefig('PostProcessing_coating/COATING_shear_power.png',dpi=300)    
# plt.show()
plt.close()


### Calculation of total shear
W_shear = np.trapz(W_dot_shear,Time.t[1:])
W_shear_c = np.trapz(W_dot_shear_c,Time.t[1:])

print('Total shear work for normal roughness = ' + str(W_shear) + 'Nm/m')
print('Total shear work for Halved roughness = ' + str(W_shear_c) + 'Nm/m')


