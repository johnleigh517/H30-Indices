"""

Author: Sean Blake, Trinity College Dublin
Editor: John Malone Leigh, Dublin Institute for Advanced Studies


Date: January 2015
Date Edited: September 2019


Email: blakese@tcd.ie, jmalonel@tcd.ie


The purpose of this module is to provide functions which will facilitate the calculation of k-indices using the FMI method, as well as the time derivate mag fields. See http://swans.meteo.be/sites/default/files/documentation/TN-RMI-2010-01_K-LOGIC.pdf for more details.


Modules Edited: Send_emailz, do_other_plots, do_k_plots


Modules Added: nan_helper, create_folder and archive_maker


Script is upgraded to work for python 3


With the functions in this module, the FMI method can be used as follows:


1) Get minutely time, bx, by arrays using minutely()



2) Get initial K-indices using k_index_func()



3) Use these with initial_smooth()



4) Smooth the initial_smooth



5) Subtract this from minutely data



6) Get second K-indices from this subtracted data



7) Repeat steps 3-6 to get third set of K-indices



These are the final K-Indices needed. Phew!

"""

from numpy import arange

import numpy

import numpy as np

import datetime

import os

import math

#import os.listdir

#import os.path.getsize



from time import strptime

from datetime import timedelta

from scipy import interpolate

from scipy.interpolate import InterpolatedUnivariateSpline



import calendar

import matplotlib.dates as mdates

from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator



from matplotlib import pyplot as plt

import matplotlib.font_manager



##########################################################################

##########################################################################

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]





def time2float(x):

    """converts datetime to float, so that interpolation/smoothing can

       be performed"""

    if (type(x) == numpy.ndarray) or (type(x) == list):

        emptyarray = []

        for i in x:

            z = (i - datetime.datetime(1970, 1, 1, 0)).total_seconds()

            emptyarray.append(z)

        emptyarray = numpy.array([emptyarray])

        return emptyarray[0]

    else:

        return (x - datetime.datetime(1970, 1, 1, 0)).total_seconds()



##########################################################################

##########################################################################



def float2time(x):

    """converts array back to datetime so that it can be plotted with time

       on the axis"""

    if (type(x) == numpy.ndarray) or (type(x) == list):

        emptyarray = []

        for i in x:

            z = datetime.datetime.utcfromtimestamp(i)

            emptyarray.append(z)

        emptyarray = numpy.array([emptyarray])

        return emptyarray[0]

    else:

        return datetime.datetime.utcfromtimestamp(x)



##########################################################################

##########################################################################



def timedatez(date, timez):

    """creating 'timedate' array for date specified, from date + time columns"""

    timedate = [] 

    for i in arange(0, len(date)):

        a = date[i] + timez[i]
        try:
            c = datetime.datetime(*strptime(a, "%Y/%m/%d%H:%M:%S")[0:6])    
        except (ValueError):
            c = datetime.datetime(*strptime(a, "%Y-%m-%d%H:%M:%S.000")[0:6]) 
            
        timedate.append(c)

    return timedate




##########################################################################

##########################################################################
def mag_filter(Bx,By,Bz):
    """
    Filters out noise from cars
    """
    if len(Bx)==len(By) and len(Bx)==len(Bz):
        dF=[]#rate of change of full field
        for i in range(0,len(Bx)-1):

            dF.append(abs(Bx[i+1]-Bx[i])+abs(By[i+1]-By[i])+abs(Bz[i+1]-Bz[i]))
        Bxnew=[]
        Bynew=[]
        print('length',len(dF))
        Bznew=[]
        for j in range(0,len(dF)-60,60):
            if max(dF[j:j+60])>15:
                array=np.full(shape=60,fill_value=99999.99,dtype=float)
                for k in array:
                    Bxnew.append(k)
                    Bynew.append(k)
                    Bznew.append(k)
            else:
                for l in range(0,60):
                    Bxnew.append(Bx[j+l])
                    Bynew.append(By[j+l])
                    Bznew.append(Bz[j+l])

    else:
        print('Error in Mag Filter')
    return Bxnew,Bynew,Bznew


def minute_bin(timedate_float, bx, by, bz, n):



    # Gets the start of the day in seconds

    day_seconds = int(timedate_float[0])-int(timedate_float[0])%(24*3600)



    # Creates array of minutes

    minutes = arange(0, n * 1440)

    minutes = (minutes * 60) + day_seconds



    # master is numpy array with columns for bx, by, bz, count and times

    master = np.zeros((n*1440, 5))

    master[:,-1] = minutes



    # loop over times

    for i, v in enumerate(timedate_float):

        # check which master row it belongs to

        index = int((v - day_seconds)/60) #- 1

        # add to each column

        try:

            master[index][3] += 1

            master[index][0] += bx[i]

            master[index][1] += by[i]

            master[index][2] += bz[i]

        except:

            continue



    # now make empty arrays which will be filled

    minute_bx, minute_by, minute_bz, minute_time = [], [], [], []

    for i, v in enumerate(master):

        if v[3] == 0:   # if count = 0, ignore

            continue

        else:           # otherwise, add average to respective array

            minute_bx.append(v[0]/v[3])

            minute_by.append(v[1]/v[3])

            minute_bz.append(v[2]/v[3])

            minute_time.append(v[4])

    

    return minute_time, minute_bx, minute_by, minute_bz



########################################################################## 

def data_read(saveaddress):

    """ Reads the data in from Magnetometer file. If the file has not finished a row, 

        it copies the complete rows to a new file, and then reads the data from 

        there.



        Returns date, time, bx, by, bz"""



    try:

        datez, timez = np.loadtxt(saveaddress, dtype = str, usecols = (0, 1), 

                        unpack = True, skiprows = 1)    # getting date from file 

        bx, by, bz = np.loadtxt(saveaddress, usecols = (2, 3, 4), unpack = True,

                        skiprows = 1)       # getting everything else



    except (IndexError, ValueError):

        datez, timez = np.loadtxt(saveaddress, dtype = str, usecols = (0, 1), 

                        unpack = True, skiprows = 2)    # getting date from file

        bx, by, bz = np.loadtxt(saveaddress, usecols = (2, 3, 4), unpack = True,

                         skiprows = 3)          # getting everything else



    return (datez, timez, bx, by, bz)

def data_read_intermag(saveaddress):

    """ Reads the data in from Magnetometer file. If the file has not finished a row, 

        it copies the complete rows to a new file, and then reads the data from 

        there.



        Returns date, time, bx, by, bz"""



    try:

        datez, timez = np.loadtxt(saveaddress, dtype = str, usecols = (0, 1), 

                        unpack = True, skiprows = 23)    # getting date from file 

        bx, by, bz = np.loadtxt(saveaddress, usecols = (3, 4, 5), unpack = True,

                        skiprows = 23)       # getting everything else



    except (IndexError, ValueError):

        datez, timez = np.loadtxt(saveaddress, dtype = str, usecols = (0, 1), 

                        unpack = True, skiprows = 24)    # getting date from file

        bx, by, bz = np.loadtxt(saveaddress, usecols = (3, 4, 5), unpack = True,

                         skiprows = 24)          # getting everything else



    return (datez, timez, bx, by, bz)
##########################################################################



def clean2(minute_time, minute_bx, minute_by, minute_bz, sigma_multiple, n):

    #This module is not used, but left for possible future use

    days = int(minute_time[0])/(24*3600)

    day_start = days*24*3600        # start of day in seconds



    # lists to be filled

    bx_chunks, by_chunks, bz_chunks, time_chunks = [], [], [], []



    first_hour = (int(minute_time[0])%(24*3600))/(3600)



    start = 0

    # next code block breaks inputs into hour-long sub-lists

    for i in arange(first_hour, 72, 1):     # loop over 72 hours

        for j in arange(start, len(minute_time), 1):    # loop over data

    

            # if data element is in the next hour, we record its index (j)

            if (minute_time[j] - day_start)/(60*60) > i+1:  



                # get slice of hour for time, bx, by and bz

                time_chunks.append(minute_time[start:j])

                bx_chunks.append(minute_bx[start:j])

                by_chunks.append(minute_by[start:j])

                bz_chunks.append(minute_bz[start:j])



                start = j

                break



    for k in range(n):      # loop over n times

        # loop over bx, by, bz

        for arrayz in (bx_chunks, by_chunks, bz_chunks):

            for index1, value1 in enumerate(arrayz):



                sigma = numpy.std(value1)   # get sigma value for hour chunk

                meanz =np.mean(value1)        # get mean value for hour chunk
                





                for index2, value2 in enumerate(value1[::-1]):

                    # if value > x sigma away from std, delete

                    if abs(value2 - meanz) >= sigma * sigma_multiple:

                        count += 1

                

                        length = len(value1)

                        # delete entries

                        del bx_chunks[index1][length- index2 -1]

                        del by_chunks[index1][length- index2 -1]

                        del bz_chunks[index1][length- index2 -1]

                        del time_chunks[index1][length- index2 -1]



    # finally, collapse 2d lists into 1-d

    time_out = [item for sublist in time_chunks for item1 in time_chunks]

    bx_out = [item for sublist in bx_chunks for item1 in time_chunks]

    by_out = [item for sublist in by_chunks for item1 in time_chunks]

    bz_out = [item for sublist in bz_chunks for item1 in time_chunks]



    return time_out, bx_out, by_out, bz_out



##########################################################################

##########################################################################



def k_index(minute_time, minute_bx, minute_by, k9):



    # lists to be populated

    timestamp, variation, k_index, order = [], [], [], []



    # start of the day in seconds

    day_seconds = int(minute_time[0])-int(minute_time[0])%(24*3600)



    #loop over minute_array and sort them according to 3-hour block

    start = 0

    hour_block1 = int((minute_time[0] - day_seconds)/(3*60*60))

    for index, value in enumerate(minute_time):

        hour_block2 = int((value - day_seconds)/(3*60*60))



        # if hr1 and hr2 and not equal, we have entered a new 3-hr block

        if hour_block2 != hour_block1:

            try:
               
                varx = max(minute_bx[start:index-1]) - min(minute_bx[start:index-1])

                vary = max(minute_by[start:index-1]) - min(minute_by[start:index-1])

    

                # append max variation for that block

                variation.append(max(varx, vary))   

                timestamp.append(day_seconds + (hour_block1*3*60*60))

                order.append(hour_block1+1)

    

                hour_block1 = hour_block2

                start = index

            except:
                
                continue



    # add last entry

    varx = max(minute_bx[start:-1]) - min(minute_bx[start:-1])

    vary = max(minute_by[start:-1]) - min(minute_by[start:-1])
    #print(variation)
    variation.append(max(varx, vary))
    
    #print(variation)
    timestamp.append(day_seconds + (hour_block1*3*60*60))

    order.append(hour_block1+1)



    # now to use these variations to calculate the k-index value

    niemegk = numpy.array([500, 330, 200, 120, 70, 40, 20, 10, 5, 0])   # reference

    thresh = niemegk * k9/500.0 


    #print(variation)
    k_index = []        # k_index list to be populated

    for i in variation:

        for index, j in enumerate(thresh):

            if i >= j:

                z = 9-index

                if z == 0:

                    z = 0.25
                    
                    if i <0.5: #this is to cut out 99999(invalid) data, most are slighly >0 due to fmi_smoothing
                        z=0

                k_index.append(z)

                break



    return k_index, timestamp, order

def h_index(minute_time, minute_bx, minute_by, k9,mode):
    """
    An adapted form of K index, Hp index is used to calculate local indices 
    at 30 or 60 minute scales, with no upper limit 

    """
    
    if mode=='30':
        mult=0.5
    if mode=='60':
        mult=1
    #factor to mutiply by based on whether using 30 min or 60 min
    # lists to be populated

    timestamp, variation, k_index, order = [], [], [], []



    # start of the day in seconds

    day_seconds = int(minute_time[0])-int(minute_time[0])%(24*3600)



    #loop over minute_array and sort them according to 3-hour block

    start = 0

    hour_block1 = int((minute_time[0] - day_seconds)/(mult*60*60))

    for index, value in enumerate(minute_time):

        hour_block2 = int((value - day_seconds)/(mult*60*60))



        # if hr1 and hr2 and not equal, we have entered a new 3-hr block

        if hour_block2 != hour_block1:

            try:
               
                varx = max(minute_bx[start:index-1]) - min(minute_bx[start:index-1])

                vary = max(minute_by[start:index-1]) - min(minute_by[start:index-1])

    

                # append max variation for that block

                variation.append(max(varx, vary))   

                timestamp.append(day_seconds + (hour_block1*mult*60*60))

                order.append(hour_block1+1)

    

                hour_block1 = hour_block2

                start = index

            except:
                
                continue



    # add last entry

    varx = max(minute_bx[start:-1]) - min(minute_bx[start:-1])

    vary = max(minute_by[start:-1]) - min(minute_by[start:-1])
    #print(variation)
    variation.append(max(varx, vary))
    
    #print(variation)
    timestamp.append(day_seconds + (hour_block1*mult*60*60))

    order.append(hour_block1+1)



    # now to use these variations to calculate the k-index value
    if mode=='30':
        niemegk = [267*1.3*1.25,267*1.3,267,190,119,65.7,33.9,17.9,8.89,
                               4.46,2.16,0]
        #niemegk = [267,190,119,65.7,33.9,17.9,8.89,
        #                       4.46,2.16,0]
    if mode=='60':
        niemegk = [337*1.3*1.25,337*1.3,337,218,144,82.7,44.7,24.3,
                               12.4,6.11,2.97,0]
        #niemegk = [337,218,144,82.7,44.7,24.3,
        #                       12.4,6.11,2.97]
    #niemegk=np.array(niemegk)
    thresmax=np.max(niemegk)
    if np.max(variation)>thresmax:
        print('here')
        #need to add more h values as max threshold is exceeded
        
        for i in range(100):
            thresmax=thresmax*1.2
            niemegk.insert(0,thresmax)
            #niemegk=[i for i in niemegk]
            if np.max(variation)>thresmax:
                continue
            else:
                break
                #breaks the loop when a sufficiently large h value is achieved
    print(niemegk)
    niemegk=np.array(niemegk)
    print(thresmax)
    print(niemegk)
    thresh = niemegk * k9/500.0 
    print(thresh)

    #print(variation)
    h_index = []        # k_index list to be populated

    for i in variation:

        for index, j in enumerate(thresh):

            if i >= j:

                z = (len(niemegk)-1)-index

                if z == 0:

                    z = 0.25 #setting a k of 0 to slightly above so can see in a plot
                    
                    if i <0.5: #this is to cut out 99999(invalid) data, most are slighly >0 due to fmi_smoothing
                        z=0

                h_index.append(z)

                break



    return h_index, timestamp, order

#######################################################################

def fmi_smoothed(minute_time, minute_bx, minute_by, minute_time_prev, 

              minute_bx_prev, minute_by_prev, k_index, k_time):



    days = int(minute_time[0])/(24*3600)

    day_seconds = days*24*3600      # start of day in seconds

    data_start = minute_time[0]

    data_end = minute_time[-1]



    # extra minutes, depending on hour  #m value

    extra_time = [120, 120, 120, 60, 60, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 

                          60, 60, 60, 120, 120, 120]



    indices = []

    # loop over minute_array and sort them according to hour block

    start = 0       # index which 

    hour_block1 = int((minute_time[0] - day_seconds)/(60*60))

    for index, value in enumerate(minute_time):

        hour_block2 = int((value - day_seconds)/(60*60))



        # if hr1 and hr2 and not equal, we have entered a new hr block

        if hour_block2 != hour_block1:

            x = (hour_block1)%24    # find what hour of the day it is

            m = extra_time[x]       # get extra minutes from extra_time



            y = k_time.index((hour_block1)/3 + 1)   # find corresponding k-index

            n = k_index[y] ** 3.3       # get extra minutes for k-index

            nm = int(n + m)



            start_time = minute_time[start] - (nm * 60) # index of start of hour

            end_time = minute_time[index] + (nm * 60)   # index of end of hour



            indices.append([start,index, nm, start_time, end_time, hour_block1])

            hour_block1 = hour_block2   # move on to new hour

            start = index               # change index



    # add last index manually and values as above

    x = (hour_block1)%24

    m = extra_time[x]



    y = k_time.index((hour_block1)/3 + 1)

    n = k_index[y] ** 3.3

    nm = int(n + m)



    start_time = minute_time[start] - (nm * 60)

    end_time = minute_time[index] + (nm * 60)



    indices.append([start, len(minute_time), 

                    nm, start_time, end_time, hour_block1])



    # now to add the data

    smoothed_time, smoothed_bx, smoothed_by = [], [], []



    # add point at half-eleven from previous day

    index_start = len(minute_time_prev) - 70



    # get index of prev time where it is > 11:30 of day before

    while minute_time[index_start] < (day_seconds - (30*60)):

        index_start += 1



    # add time to list, mean bx and by values for specified time.

    smoothed_time.append(day_seconds - (30*60)) 

    smoothed_bx.append(numpy.mean(minute_bx_prev[index_start:]))    

    smoothed_by.append(numpy.mean(minute_by_prev[index_start:]))



    for index1, value1 in enumerate(indices): # loop over indices

        start, end, nm, start_time, end_time, hour = value1



        try:

            if start_time < data_start: 

                # if start time is before start of data

                # go on to previous data

    

                # index_start = len(prev_data) - extra time

                index_start = len(minute_time_prev) - nm    

                while minute_time_prev[index_start] < start_time:

                    index_start += 1

                # loop until minute_time_prev[index_start] > start_time

            

                # now do the same for the end data

                index_end = end + nm

                while minute_time[index_end] > end_time:

                    index_end -= 1

                

                # add mean of magnetic data

                smoothed_bx.append(numpy.mean(minute_bx_prev[index_start:] 

                                    + minute_bx[start:index_end]))

                smoothed_by.append(numpy.mean(minute_by_prev[index_start:] 

                                    + minute_by[start:index_end]))



                # add time to list

                smoothed_time.append(((hour+0.5)*(3600))+ day_seconds)  

    

            else:

                index_start = start - nm

                while minute_time[index_start] < start_time:

                    index_start += 1

            

                index_end = end + nm

                if index_end < len(minute_time):

                    while minute_time[index_end] > end_time:

                        index_end -= 1

    

                else:

                    index_end = len(minute_time)

    

                smoothed_bx.append(numpy.mean(minute_bx[index_start:index_end]))

                smoothed_by.append(numpy.mean(minute_by[index_start:index_end]))



                # add time to list

                smoothed_time.append(((hour+0.5)*(3600))+ day_seconds)  

    

                a = ((hour+0.5)*(3600))+ day_seconds



        except:

            continue

    return smoothed_time, smoothed_bx, smoothed_by





##########################################################################

##########################################################################



def fmi_smoothed2(minute_time, minute_bx, minute_by, k_index, k_time):

    extra_time = [120, 120, 120, 60, 60, 60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 

                              60, 60, 60, 120, 120, 120]





    start_time = minute_time[0]

    start_day = start_time - start_time%(24*60*60)

    end_time = minute_time[-1]



    # number of hours to look at:

    hour_number = int((end_time - start_day)/3600) + 1



    hour_indices =[[] for i in range(hour_number)]



    k_index_thing, k_time_thing = [], []

    for i, j in list(zip(k_index, k_time)):

        if i == 0.25:

            k_index_thing.extend((0, 0, 0))

        else:

            k_index_thing.extend((i, i, i))



        l = j*3 

        k_time_thing.extend((l-3, l-2, l-1))



    master = np.zeros((hour_number, 7))

    blank =[]

    for i in range(hour_number):

        try:

            k = k_index_thing[k_time_thing.index(i)]

            n = k**3.3

        except:

            n = 0

        m = extra_time[i%24]



        master[i][0] = i    # hour number

    



        avg_time =start_day + (i*60*60) + (30*60) # half an hour time

        blank.append(avg_time)



        master[i][1] = avg_time



        start_time = avg_time - (30*60) - ((n + m)*60)

        master[i][2] = start_time



        end_time = avg_time + (30*60) + ((n + m)*60)

        master[i][3] = end_time





    for index, value in enumerate(minute_time):



        for i, v in enumerate(master):

            start = master[i][2]

            end = master[i][3]



            if start <= value <= end:

                master[i][4] += minute_bx[index]

                master[i][5] += minute_by[index]

                master[i][6] += 1

            



    nz = (master == 0).sum(1)

    q = master[nz == 0, :]



    smoothed_time = master[:,1]

    smoothed_bx = master[:,4]/master[:,6]

    smoothed_by = master[:,5]/master[:,6]



    clean_time, clean_bx, clean_by = [], [], []

    for t, x, y in zip(smoothed_time,smoothed_bx, smoothed_by):

        if (math.isnan(x) == False) and (math.isnan(y) == False):

            clean_time.append(t)

            clean_bx.append(x)

            clean_by.append(y)





    smt, smx, smy = [], [], []

    for i in clean_time:

         x = np.interp(i, clean_time, clean_bx)

         y = np.interp(i, clean_time, clean_by)

         

         smx.append(x)

         smy.append(y)

         



    return clean_time, smx, smy





##########################################################################

##########################################################################

def smoothed(time_array, y_value, z):   # z = smoothing order <=5

    """ Smooths the rough Sr curve """

    x = time_array

    y = y_value

    xi = np.linspace(x[0], x[-1], 10*len(x))

    ius = InterpolatedUnivariateSpline(x, y, k=z)

    yi = ius(xi)

    return xi, yi



##########################################################################

##########################################################################



def subtracted(minute_time, smooth_time, minute_bx, smooth_bx, minute_by, 

               smooth_by):

    """ Subtracts smooth data from original minute array"""

    subtracted_bx, subtracted_by = [], []

    for index, value in enumerate(minute_time):

        x = numpy.interp(value, smooth_time, smooth_bx)

        y = numpy.interp(value, smooth_time, smooth_by)

        subtracted_bx.append(minute_bx[index]-x)

        subtracted_by.append(minute_by[index]-y)



    return subtracted_bx, subtracted_by



##########################################################################

##########################################################################



def colored(k_index, barlist):

    """Colours yer k-index plots so it looks nice.



    EXAMPLE:

    barlist = ax3.bar(k_timestamps, k_index, width = 0.124)

    colored(k_index)"""

    for i in arange(0, len(k_index), 1):

        if k_index[i] >= 8:

            barlist[i].set_color('deeppink')

            barlist[i].set_edgecolor('k')

            continue

        if k_index[i] >= 6:

            barlist[i].set_color('r')

            barlist[i].set_edgecolor('k')           

            continue

        if k_index[i] >= 5:

            barlist[i].set_color('orange')

            barlist[i].set_edgecolor('k')           

            continue

        if k_index[i] >= 4:

            barlist[i].set_color('g')

            barlist[i].set_edgecolor('k')           

            continue

        if k_index[i] >= 2:

            barlist[i].set_color('c')

            barlist[i].set_edgecolor('k')           

            continue

        if k_index[i] >= 0:

            barlist[i].set_color('b')

            barlist[i].set_edgecolor('k')           

            continue



##########################################################################

##########################################################################



def slope_refined(x, y):    #want it to return y array

    xx, yy = [], []

    new_x = arange(x[0], x[-1], 60)

    for i in new_x:

        y2 = numpy.interp(i, x, y)

        y1 = numpy.interp(i-60, x, y)

        deriv = (y2 - y1) / 1.0

        xx.append(i)

        yy.append(deriv)

    return xx, yy



##########################################################################

##########################################################################



def do_k_plots(k_index3, k_time3, k_timestamp3, minute_time, sitefull, save_address,save_address2):

    # 3-Day K-Index Plot
    

    nowz = datetime.datetime.utcnow() #houdini needs to read utc time
    
    year_str=str(nowz.year)
    month_str="%02d" %(nowz.month)
    day_str="%02d" %(nowz.day)
    date_str= day_str+'/'+month_str+'/'+year_str
    
    #timestamp_string = "Plot updated {}/{}/{} {} UT".format(nowz.day, nowz.month, nowz.year, str(nowz)[10:19]) 
    timestamp_string = "Plot updated "+str(date_str)+ " "+ str(nowz)[10:19] + " UT"


    end = float2time(minute_time[-1]- (minute_time[-1]%(24*3600)) + (24*3600))

    start = end -datetime.timedelta(3)

    middle1 = end - datetime.timedelta(2)

    middle2 = end - datetime.timedelta(1)



    print ("START: ", start)

    print ("MIDDLE1: ", middle1)

    print ("MIDDLE2: ", middle2)

    print ("END: ", end)

 

#   start = minute_time[0] - minute_time[0]%(24*3600)

#   middle1 = float2time(start + (1*24*3600))

#   middle2 = float2time(start + (2*24*3600))

#   end = float2time(start + (3*24*3600))

#   start = float2time(start)

    

    

    #plt.clf()
    #plt.style.use('classic')
    fig1 = plt.figure()

    plt.subplots_adjust(bottom = 0.2, left = 0.07, right = 0.94)

    

    a = plt.subplot(111)

    bartimenew=[]
    for m in range(0,len(k_timestamp3)): #Adjusting figures to centre bar plots, 5400=1h 30m in sec
        bartime2=k_timestamp3[m]+5400
        bartimenew.append(bartime2)
    #print(len(k_timestamp3))
    #print(len(k_index3))
    #print(bartimenew)
    #print(k_index3)
    bartime=float2time(bartimenew)
    

        
        
    barlist = a.bar(bartime, k_index3, width = 0.124)

    colored(k_index3, barlist)

    

    plt.xlabel('{}-{}-{}                         {}-{}-{}                         {}-{}-{}'.format(start.day, start.strftime("%B")[0:3], start.year, middle1.day, middle1.strftime("%B")[0:3], middle1.year, middle2.day, middle2.strftime("%B")[0:3], middle2.year), fontsize = 14)

    plt.ylabel("K-Index", fontsize = 14)

    plt.title(str(sitefull)+" 3-Day Local K-Index", fontsize = 16)

    

    plt.xlim([start, end])

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:00'))

    plt.gca().xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,6)))

    plt.axvline(middle1, color = 'k', lw = 1)

    plt.axvline(middle2, color = 'k', lw = 1)

    

    plt.ylim([0, 9])

    plt.grid(True)

    

    plt.figtext(0.65, 0.02, timestamp_string, fontsize = 11, style='italic')



    

    fig1.set_size_inches(9, 4)
    

    plt.figtext(0.02, 0.02, "SolarTsunamis", fontsize = 11, style = 'italic')
    
    #date_str= "%d%02d%02d" %(nowz.year,nowz.month,nowz.day)    

    #fig1.savefig("C:\\Users\\john_\\Documents\\MagIE\\MagIE_all_files\\K_Indices\\"+str(site_l)+str(middle2str)+"_kindex.png")
    fig1.savefig(str(save_address)+'_kindex.png')
    fig1.savefig(str(save_address2)+'_kindex.png')
    #fig1.savefig('/mnt/data.magie.ie/magnetometer_archive/'+year_str+'/'+month_str+'/'+day_str+'/png/'+str(site_l)+'_kindex'+year_str+month_str+day_str+'_kindex.png')
    


##########################################################################


def do_h_plots(h_index3, h_time3, h_timestamp3, minute_time, sitefull, save_address,save_address2,mode):

    # 3-Day K-Index Plot
    

    nowz = datetime.datetime.utcnow() #houdini needs to read utc time
    
    year_str=str(nowz.year)
    month_str="%02d" %(nowz.month)
    day_str="%02d" %(nowz.day)
    date_str= day_str+'/'+month_str+'/'+year_str
    
    #timestamp_string = "Plot updated {}/{}/{} {} UT".format(nowz.day, nowz.month, nowz.year, str(nowz)[10:19]) 
    timestamp_string = ""#"Plot updated "+str(date_str)+ " "+ str(nowz)[10:19] + " UT"


    end = float2time(minute_time[-1]- (minute_time[-1]%(24*3600)) + (24*3600))

    start = end -datetime.timedelta(3)

    middle1 = end - datetime.timedelta(2)

    middle2 = end - datetime.timedelta(1)

    middle3 = end - datetime.timedelta(hours=30)

    print ("START: ", start)

    print ("MIDDLE1: ", middle1)

    print ("MIDDLE2: ", middle2)

    print ("END: ", end)

 

#   start = minute_time[0] - minute_time[0]%(24*3600)

#   middle1 = float2time(start + (1*24*3600))

#   middle2 = float2time(start + (2*24*3600))

#   end = float2time(start + (3*24*3600))

#   start = float2time(start)

    

    

    #plt.clf()
    #plt.style.use('classic')
    fig1 = plt.figure()

    plt.subplots_adjust(bottom = 0.2, left = 0.07, right = 0.94)

    

    a = plt.subplot(111)

    bartimenew=[]
    for m in range(0,len(h_timestamp3)): #Adjusting figures to centre bar plots, 5400=1h 30m in sec
        if mode=='30':
            bartime2=h_timestamp3[m]+15*60
        if mode=='60':
            bartime2=h_timestamp3[m]+30*60
        bartimenew.append(bartime2)
    #print(len(k_timestamp3))
    #print(len(k_index3))
    #print(bartimenew)
    #print(k_index3)
    bartime=float2time(bartimenew)
    

        
        
    if mode=='30':
        barlist = a.bar(bartime, h_index3, width = 0.124/6)
    if mode=='60':
        barlist = a.bar(bartime, h_index3, width = 0.124/3)

    colored(h_index3, barlist)

    

    #plt.xlabel('{}-{}-{}                                                                                                                     {}-{}-{}                                                                                                                   {}-{}-{}'.format(start.day, start.strftime("%B")[0:3], start.year, middle1.day, middle1.strftime("%B")[0:3], middle1.year, middle2.day, middle2.strftime("%B")[0:3], middle2.year), fontsize = 14)
    plt.xlabel('11/05/2024', fontsize = 14)
    plt.xlabel('10/05/2024                                                                                                                   11/05/2024                                  ', fontsize = 14)
    #plt.ylabel("H30-Index", fontsize = 14)
    plt.ylabel("Hpo30-Index", fontsize = 14)

    #plt.title(str(sitefull)+" 3-Day Local H"+str(mode)+"-Index", fontsize = 16)

    #plt.title(str(sitefull)+" Local H"+str(mode)+"-Index", fontsize = 16)
    plt.title("Hpo"+str(mode)+"-Index", fontsize = 16)
    plt.xlim([middle1, middle3])
    
    plt.xlim([end -datetime.timedelta(days=2,hours=10),end -datetime.timedelta(days=1,hours=10)])

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:00'))

    plt.gca().xaxis.set_major_locator(mdates.HourLocator(byhour=range(0,24,6)))

    plt.axvline(middle1, color = 'k', lw = 1)

    plt.axvline(middle2, color = 'k', lw = 1)

    

    #plt.ylim([0, np.max(h_index3)])
    plt.ylim([0, 15])
    plt.grid(True)

    

    plt.figtext(0.835, 0.02, timestamp_string, fontsize = 11, style='italic')


    plt.axvline(datetime.datetime(2024,5,10,22,45),color='black',linestyle='--',linewidth=3)
    plt.axvline(datetime.datetime(2024,5,11,8,45),color='black',linestyle='--',linewidth=3)
    plt.axvline(datetime.datetime(2024,5,11,12,45),color='black',linestyle='--',linewidth=3)

    plt.axvline(datetime.datetime(2024,10,10,15,15),color='black',linestyle='--',linewidth=3)
    plt.axvline(datetime.datetime(2024,10,10,23,15),color='black',linestyle='--',linewidth=3)
    plt.axvline(datetime.datetime(2024,10,11,10,15),color='black',linestyle='--',linewidth=3)
    

    fig1.set_size_inches(16, 4)
    #fig1.set_size_inches(28, 4)
    

    #plt.figtext(0.07, 0.02, "SolarTsunamis", fontsize = 11, style = 'italic')
    
    #date_str= "%d%02d%02d" %(nowz.year,nowz.month,nowz.day)    

    #fig1.savefig("C:\\Users\\john_\\Documents\\MagIE\\MagIE_all_files\\K_Indices\\"+str(site_l)+str(middle2str)+"_kindex.png")
    fig1.savefig(str(save_address)+'_hindex.png')
    fig1.savefig(str(save_address2)+'_hindex.png')
    #fig1.savefig('/mnt/data.magie.ie/magnetometer_archive/'+year_str+'/'+month_str+'/'+day_str+'/png/'+str(site_l)+'_kindex'+year_str+month_str+day_str+'_kindex.png')
    


##########################################################################
