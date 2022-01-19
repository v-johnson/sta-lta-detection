import numpy as np
from obspy import Stream
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
import pickle

def hourly_ppsd(st, date, paz = None, start = 0, stop = np.nan, num_of_intervals = 12, filename = None,verbose=0):
    seconds = 12/num_of_intervals * 300
    
    #Determines largest size for arrays, will pad remaining indices with 0.
    sampling_rate = st[0].stats.sampling_rate
    arrsize = seconds * sampling_rate/2
    power_sum = np.zeros(int(arrsize))
 
    if np.isnan(stop):
        stop = st[0].stats.endtime-st[0].stats.starttime
            
    # Creates one hour long pickle files and runs through data in 5 min intervals
    # returns partial trace if trace does not contain data for timeframe.
    
    for i in np.arange(0, 24):
        number_of_good_traces = 0
        fig = plt.figure()
        fig.patch.set_facecolor('white')

        tr_0 = st[0].copy()
        st_chunk = Stream(tr_0)
        starttime = UTCDateTime(date.year,date.month,date.day,i)
        endtime = starttime + (3600)
        st_chunk.trim(starttime,endtime)
        for j in np.arange(0,num_of_intervals):
            tr = st_chunk.copy()
            st_subset = Stream(tr)
            start_1 = starttime + (j * seconds)
            end_1 = start_1 + seconds
            st_subset.trim(start_1,end_1)

            if verbose > 1:
                print("")
                print("     start time is: " + str(start_1))
                print("     end time is: " + str(end_1))
                print(st_subset)

            if paz is not None:
                st_subset.simulate(paz_remove=paz)

            if len(st_subset) > 0:
                new_tr = st_subset[0]
                N = len(new_tr.data)
                if N > 1:
                    number_of_good_traces = number_of_good_traces + 1
                else:
                    continue
                T = tr[0].stats.delta

                yf = fft(new_tr.data)
                xf_orig = fftfreq(N, T)[:N//2]
                xf = xf_orig.copy()
                if verbose > 0:
                    print("Found a time window with data.")
                    print("size of xf: %d"%xf.size)
                    print("size of N: %d"%N)
                ya = 2.0/N * np.abs(yf[0:N//2])
            
#             print("     Success on window %d,%d."%i,%j)

                #Makes ya and power_sum the same size
                if len(power_sum) > len(ya):
                    ya.resize(power_sum.size)
                if len(ya) > len(power_sum):
                    power_sum.resize(ya.size)
                #Makes xf and power_sum the same size, should be redundant after similar check with ya.
                if len(power_sum) > len(xf):
                    xf.resize(power_sum.size)
                if len(xf) > len(power_sum):
                    power_sum.resize(xf.size)    
            
                power_sum = power_sum + ya**2
                plt.plot(xf, ya**2, "-k" , alpha = 0.25)
                plt.yscale("log")
                plt.xscale("log")
            if number_of_good_traces > 0:   
                avg_power = power_sum / number_of_good_traces
                if verbose > 0:
                    print("size of power_sum: %d"%len(power_sum))
                    print("size of avg_power: %d"%len(avg_power))

                if filename == None:
                    plt.show()
                else:
                    filename_ext = filename + "%s"%starttime
                    plt.savefig("hour-ppsd_plots/"+filename_ext+'.png')
                    with open("hour-pickles/"+filename_ext+ '.pkl', "wb") as f:
                        pickle.dump([xf, avg_power], f)
            plt.close(fig)

    
    

def better_ppsd(st, date, paz = None, start = 0, stop = np.nan, number_of_one_hour_intervals = 24, filename = None,verbose=0):
    seconds = 24/number_of_one_hour_intervals * 3600
    
    #Determines largest size for arrays, will pad remaining indices with 0.
    sampling_rate = st[0].stats.sampling_rate
    arrsize = seconds * sampling_rate/2
    power_sum = np.zeros(int(arrsize))
 
    
    fig = plt.figure()
    fig.patch.set_facecolor('white')
    
    if np.isnan(stop):
        stop = st[0].stats.endtime-st[0].stats.starttime
            
    # Works through trace in 1 hour increments from 0:00:00-1:00:00, 
    # returns partial trace if trace does not contain data for timeframe.
    number_of_good_traces = 0
    for i in np.arange(0,number_of_one_hour_intervals):

        tr = st[0].copy()
        st_subset = Stream(tr)
        starttime = UTCDateTime(date.year,date.month,date.day) + (i * seconds)
        endtime = starttime + seconds
        st_subset.trim(starttime,endtime)
        
        if verbose > 1:
            print("")
            print("increment is: " + str(i))
            print("     start time is: " + str(starttime))
            print("     end time is: " + str(endtime))
            print(st_subset)
        

        if paz is not None:
            st_subset.simulate(paz_remove=paz)
        
        if len(st_subset) > 0:

            new_tr = st_subset[0]
            N = len(new_tr.data)
            if N > 1:
                number_of_good_traces = number_of_good_traces + 1
            else:
                continue
            T = tr.stats.delta
            
            yf = fft(new_tr.data)
            xf_orig = fftfreq(N, T)[:N//2]
            xf = xf_orig.copy()
            if verbose > 0:
                print("Found a time window with data.")
                print("size of xf: %d"%xf.size)
                print("size of N: %d"%N)
            ya = 2.0/N * np.abs(yf[0:N//2])
            
#             print("     Success on window %d."%i)

            #Makes ya and power_sum the same size
            if len(power_sum) > len(ya):
                ya.resize(power_sum.size)
            if len(ya) > len(power_sum):
                power_sum.resize(ya.size)
            #Makes xf and power_sum the same size, should be redundant after similar check with ya.
            if len(power_sum) > len(xf):
                xf.resize(power_sum.size)
            if len(xf) > len(power_sum):
                power_sum.resize(xf.size)    
            
            power_sum = power_sum + ya**2
            plt.plot(xf, ya**2, "-k" , alpha = 0.25)
            plt.yscale("log")
            plt.xscale("log")

    avg_power = power_sum/number_of_good_traces
    
    if verbose > 0:
        print("size of power_sum: %d"%len(power_sum))
        print("size of avg_power: %d"%len(avg_power))

#     print("     Average power line successfully plotted")
    plt.plot(xf,avg_power , "-r")
    plt.grid()   
    plt.yscale("log")
    plt.xscale("log")

    if filename == None:
        plt.show()
    else:
        plt.savefig("ppsd_plots/"+filename+'.png')
        with open("pickles/"+filename+ '.pkl', "wb") as f:
            pickle.dump([xf, avg_power], f)
            
    plt.close(fig)
            
    return avg_power,xf
        
    

    
def power_time_series(st, paz = None, start = 0, stop = 1000, step = 5):
    '''
    Returns the time averaged trace RMS.
    '''
    avg_power = []
    t = []
    for i,window_start in enumerate(np.arange(start,stop,step)):
        tr = st[0].copy()
        t1 = tr.stats.starttime
        st_subset = Stream(tr)
        st_subset.trim(t1 + window_start, t1 + 5 + window_start)
        if paz is not None:
            st_subset.simulate(paz_remove=paz)
        
        new_tr = st_subset[0]
        RMS = np.sqrt(np.sum(new_tr.data**2)/len(new_tr.data))
        avg_power.append(RMS)
        t.append(new_tr.stats.starttime)
        
    return avg_power,t
