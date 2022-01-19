#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
Last Updated: 8/26/2021, V1
Updates needed to station plots.

This code plots a map of Antarctica based on sliding velocity. 
It then marks the location of several stations where icequakes have been detected. 

'''


# In[ ]:


import sys
import matplotlib.pyplot as plt
import h5py
import numpy as np
import obspy
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


ds = h5py.File('/data/fast0/datasets/ismip6_issm.h5','r')
x = ds['X']
y = ds['Y']
s = ds['SPEED']

plt.subplots()
plt.pcolormesh(x,y,s,vmax=100)


# In[ ]:


'''
Points being used:
Network
Station
Lat
Long
Conversion to EPSG 3031 [x,y]



 XF
DN1N
-82.384995
-142.157425
EPSG 3031: [-508316.13,-654313.34]

 XF
DN2N
-82.671516
-147.520294
EPSG 3031: [-428151.38,-672588.65]

 XF
DN3C
-82.811295
-152.595184
EPSG 3031: [-359959.98,-694290.35]


David's Glacier
-75.333333 
161.256389
EPSG 3031: [30680108.53, -90413806.86]

'''


# In[46]:


'''
Needs adjusting to put both figures in same scale. 

'''

ds = h5py.File('/data/fast0/datasets/ismip6_issm_5km.h5','r')
x = np.array(ds['X'])
y = np.array(ds['Y'])
u = np.array(ds['SPEED'])
m = np.array(ds['MASK'])
u[m < 1] = np.nan

# fig,ax=plt.subplots(1,2)
fig = plt.figure()
ax = fig.add_subplot(111)
plt.pcolormesh(x/1e3,y/1e3,u,vmax=100)
fig2 = plt.figure()
bx = fig2.add_subplot(111)

stations_x = [-1069367.22,-508316.13, -428151.38, -359959.98]
stations_y = [-431167.30, -654313.34, -672588.65, -694290.35]
# ax.set_aspect('equal')
# plt.axes().set_aspect("equal")
# ax.plot(wais_x,wais_y,'^',markerfacecolor="None")
plt.plot(stations_x, stations_y, '^', markerfacecolor = "None")



# In[7]:


'''
Reads and plots traces from Kamb Glacier.

'''

client = Client("IRIS")
t1 = UTCDateTime("1995-11-23T00:00:00")
t2 = UTCDateTime("1996-12-31T00:00:00")
st = client.get_waveforms("XF","DN1*", "--", "EPZ", t1, t2, attach_response = True )
# st = client.get_waveforms("XF","UP1E", "--", "EPZ", t1, t2, attach_response = True )
# st = client.get_waveforms("XF","UP2E", "--", "EPZ", t1, t2, attach_response = True )

pre_filt = (10,20,100,200)
st.remove_response(output = "VEL", pre_filt = pre_filt)


# In[8]:


st[2].plot()


# In[ ]:




