
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

# with h5py.File('Data/Time_0.05ms.h5', 'r') as f:
#     for key in f.keys():
#         print(key)

#     data = f['h0']
#     print(min(data))
#     print(max(data))
#     print(data[:15])

# time = [np.round(i,2) for i in np.arange(0.05, 7.50, 0.05 )] #ms

# for t in time:
#     FileName='Data/Time_'+str(t)+'ms.h5'
#     Data=IO.ReadData(FileName)['State']

#     if Data['WearDepthCylinder'].any() != 0:
#         print("Not equal to zero!")
#     # print(Data['State'])

time = 26.90
FileName='Data/Time_'+str(time)+'ms.h5'
Data=IO.ReadData(FileName)['State']
print(Data['WearDepthCylinder'])
# print(round(Time.t[1]*1000,4))
# for i in Data:
#     print(i)

# greek_letterz=[chr(code) for code in range(945,970)]

# print(greek_letterz)

# print(Data['WearLocationsCylinder'])
# print(Data['WearDepthCylinder'])
# print(len(Data['WearLocationsCylinder']))



# print(Grid.x)
# for x in range(Grid.Nx):
    # print(x)
## Voor in postprocessing.py
# p = np.zeros(Time.nt - 1)


# for time in range(Time.nt - 1):
#     p[time] = StateVector[time].HertzianContactPressure
    
# plt.plot(Ops.CranckAngle[1:],p )
# plt.show()
