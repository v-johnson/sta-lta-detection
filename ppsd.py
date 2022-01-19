from obspy.io.xseed import Parser
from obspy import read
from obspy import Stream
import matplotlib.pyplot as plt
import re
import os
from scipy.fft import fft, fftfreq
import numpy as np
from time import perf_counter
import pickle

import icequake_processing
import importlib
importlib.reload(icequake_processing)

t_start = perf_counter()
print('Loading Data')

#Creates a list of all data files.
pth = "/data/fast1/time/"
dirs = os.listdir(pth)

for file in dirs[40:42]:
    station_file = file
    st = read(pth + station_file)
    st = st.select(component = "Z")
    station = st[0].stats.station
    tr = st[0].copy()
    t0 = tr.stats.starttime
    parser = Parser('/data/fast1/time/RESP.YE.N303.GPZ.1000SPS.12DB')
    paz = parser.get_paz('GPZ')
    print('Finished loading data in %f seconds'%(perf_counter()-t_start))

    imax = len(st)

    for i in np.arange(0,imax):
        # each st[i] contains one UTC day of data.
        # The first and last days are frequently partial days.
        
        #Finds the date the file is running on regardless of starting on the day or day before. Saves as string for filename.
        length = (st[i].stats.endtime - st[i].stats.starttime)
        new_time = (st[i].stats.starttime + length/2)
        date_string = re.search(r'(.*T)',  str(new_time)).group(0)
        
        
        temp_st = Stream(traces = st[i])
        print('Starting PPSD Calculation...')
        
        # Daily PPSD with 1 hour intervals
        #icequake_processing.better_ppsd(temp_st, new_time, paz, number_of_one_hour_intervals = 24,
                                        #filename = ("test_PPSD_station_%s_%s" %(station, date_string)))
            
        # Hourly PPSD with 5 minute intervals
        icequake_processing.hourly_ppsd(temp_st, new_time, paz, filename = ("ppsd_station_%s_" %station), verbose = 2)
        

        
        print('... Finished day %d of %d in %f seconds.'%(i+1,imax,perf_counter()-t_start))
    print("success!")
    
    
    