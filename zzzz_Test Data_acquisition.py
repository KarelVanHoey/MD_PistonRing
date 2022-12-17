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

# time = [np.round(i,2) for i in np.arange(0.05, 7.50, 0.05 )] #ms

# for t in time:
#     FileName='Data/Time_'+str(t)+'ms.h5'
#     Data=IO.ReadData(FileName)['State']

#     if Data['WearDepthCylinder'].any() != 0:
#         print("Not equal to zero!")
#     # print(Data['State'])

time = 0.25
FileName='Data/Time_'+str(time)+'ms.h5'
Data=IO.ReadData(FileName)
for key in Data['State']:
    print(key +": " + str(Data['State'][key][0]) )


print("ring")
print(Data['State']["WearDepthCylinder"])

