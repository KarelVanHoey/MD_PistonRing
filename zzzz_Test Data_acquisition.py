import h5py
import numpy as np
from IOHDF5 import IOHDF5
IO=IOHDF5()

# with h5py.File('Data/Time_0.05ms.h5', 'r') as f:
#     for key in f.keys():
#         print(key)

#     data = f['h0']
#     print(min(data))
#     print(max(data))
#     print(data[:15])

time = 0.05 #ms
FileName='Data/Time_'+str(time)+'ms.h5'
Data=IO.ReadData(FileName)

# print(Data['State'])

for key in Data['State']:
    print(key +": " + str(Data['State'][key][0]) )

# print(Data['State']["Viscosity"][0])
# print("ring")
# print(Data['State']["WearDepthRing"])

