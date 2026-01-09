"""
Original Authors: Sean Blake, Trinity College Dublin, Johnny Malone Leigh, Dublin Institute for Advanced Studies
Editor: Johnny Malone Leigh, University of Otago

Date First Created: January 2015
Date Last Edited: July 2025


Email: john.m.leigh@otago.ac.nz


The purpose of this script is to calculate h30 and k indices using modules from k_index_pre.py.


The h30 and k-indices (and other plots) are then created and then saved to the server. This script is naturally designed to operate with real time files.
But currently files are older files are manually input, related to data from Malone-Leigh et al., 2026.

"""
import numpy as np
import datetime
import os

#each module needed to run functions in k_index_pre_houdini
from time import strptime
from datetime import timedelta
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import pyplot as plt
import matplotlib.font_manager
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
#Importing files from kindex_pre
from kindex_pre import (nan_helper, time2float,float2time, timedatez,mag_filter, minute_bin, 
                        clean2, k_index, fmi_smoothed, fmi_smoothed2, smoothed,subtracted, colored, 
                        slope_refined, do_k_plots,data_read, h_index,data_read_intermag,do_h_plots)
import seaborn as sns
plt.close('all')
sns.set()
##########################################################################
##########################################################################

sites_obs=['eyr'] 
sitefull_list=['Eyrewell']

k_thres=[750 ,540, 480] 
k_thres=[540]
#the k_threshold is the minimum value to denote a kp 9 storm(in nT)
#needs to match each site


nowz = datetime.datetime.utcnow() #need utc time

site_count=0 #added to read different variables in lists above 
year_str=str(nowz.year)
month_str="%02d" %(nowz.month)
day_str="%02d" %(nowz.day)
date_str= day_str+'/'+month_str+'/'+year_str
intensity_list=[]
k_stamp_list=[]
condition_list=[]




#Creating a loop to read data for each site
for k1, k2 in enumerate(sites_obs):
    
    k2_l = k2.lower() #
    k_max=k_thres[site_count] #ensures next k_thres value ran with each iteration of sites_obs
    sitefull_name=sitefull_list[site_count]
    site_count=site_count+1 
    nowz = datetime.datetime.utcnow() -datetime.timedelta(days=1)#need utc time
    year_str=str(nowz.year)
    month_str="%02d" %(nowz.month)
    day_str="%02d" %(nowz.day)
    
    folder='/home/johnnyl/Documents/Solar-Tsunami/Kptest/'

    save_address=folder +str(k2_l)
    save_address2=folder +str(k2_l)+year_str+month_str+day_str

    
    #save_address='/mnt/data.magie.ie/magnetometer_live/'+str(k2_l)    
    #save_address2="/mnt/data.magie.ie/magnetometer_archive/"+year_str+"/"+month_str+'/'+day_str+'/png/'+str(k2_l)+year_str+month_str+day_str
    print(save_address2)
    #folder="/mnt/data.magie.ie/magnetometer_archive/"
    file_full=[]
    #folder_live="/mnt/data.magie.ie/magnetometer_live/"+str(k2_l)+"_mag_realtime.txt"   

    #file_full.append(folder_live)
    for j in range(0,4,1): #making strings to find files in archive


        
        time=nowz-datetime.timedelta(days=j)
        date_str= "%d%02d%02d" %(time.year,time.month,time.day)
        year_str=str(time.year)
        month_str="%02d" %(time.month)
        day_str="%02d" %(time.day)
        #file_str=year_str+'/'+month_str+'/'+day_str+'/txt/'+str(k2_l)+date_str+'.txt'
        if k2_l!='eyr':
            file_str=str(k2)+date_str+'.txt'
        else:
            file_str=str(k2_l)+date_str+'psec.sec'
        file_str2=folder+file_str
        file_full.append(file_str2)

        # manually setting files
        file_full=['/home/johnnyl/Documents/Solar-Tsunami/Oct2024/awa20241012.txt',
               '/home/johnnyl/Documents/Solar-Tsunami/Oct2024/awa20241011.txt',
               '/home/johnnyl/Documents/Solar-Tsunami/Oct2024/awa20241010.txt',
               '/home/johnnyl/Documents/Solar-Tsunami/Oct2024/awa20241009.txt']
        
        file_full=['/home/johnnyl/Downloads/May2024_Try3/awa20240512.txt',
                   '/home/johnnyl/Downloads/May2024_Try3/awa20240511.txt',
                   '/home/johnnyl/Downloads/May2024_Try3/awa20240510.txt',
                   '/home/johnnyl/Downloads/May2024_Try3/awa20240509.txt']
        
        file_full=['/home/johnnyl/Downloads/eyr20240512psec.sec',
                   '/home/johnnyl/Downloads/eyr20240511psec.sec',
                   '/home/johnnyl/Downloads/eyr20240510psec.sec',
                   '/home/johnnyl/Downloads/eyr20240509psec.sec'
            ]
        
    #try: #try loop added in case site fails due to missing data
    #Add this try loop out when processing in real time
    #And comment out the placeholder 1==1
    if 1==1:
        file_list = []
        for k in file_full:
        	print (k) 
        	#if ".txt" not in k:
        	#	continue
        	l = os.path.getmtime(k)
        	file_list.append((l, k))
        	#file_list = sorted(file_list, reverse = True)
    
        print(len(file_list))
        

        
    ########################################
    #Normal K index section
    
            
                
        current_file = file_list[0]
        current_day_seconds = current_file[0] - current_file[0]%(24*60*60)
        
        datez, timez, bx, by, bz, tempfg = [], [], [], [], [], []
        
        
        for i in range(3, -1, -1): #reads last 4 files in folder, from end to start   
    
            filename = file_list[i][1]
            if k2_l!='eyr':
                date1, time1, bx1, by1, bz1 = data_read(filename)
            else:
                date1, time1, bx1, by1, bz1 = data_read_intermag(filename)
            #if k2_l == 'val' or 'dun':
            #    temp1 = np.loadtxt(filename, usecols = (5,), skiprows = 2)
            #    
            #else:
            #    temp1 = np.loadtxt(filename, usecols = (10,), skiprows = 2)
            #tempfg = np.concatenate((tempfg, temp1), axis = 0)
            #Adding a quick check to ensure that same file isn't repeated two days in a row
            if i<3:
                #print(bx[0])
                #print(bx1[0])
                if bx_record[0] == bx1[0] and by_record[0]==by1[0] and bz_record[0]==bz1[0]:
                    print(filename)
                    print('copy')
                    bx2=np.ones(len(bx1))
                    bx1=np.array([99999.99*i for i in bx2])
                    by1=bx1
                    bz1=bx1
            bx_record=bx1
            by_record=by1
            bz_record=bz1
            bx = np.concatenate((bx,bx1), axis = 0)
            by = np.concatenate((by,by1), axis = 0)
            bz = np.concatenate((bz,bz1), axis = 0)
            datez = np.concatenate((datez,date1), axis = 0)
            timez = np.concatenate((timez,time1), axis = 0)
    
           
            timedate = timedatez(datez, timez)
            timedate_float = time2float(timedate)
            print(bx)
        if k2_l!='awa':
            print('filtering')
            print(len(bx))
            bx,by,bz=mag_filter(bx,by,bz)
            bx=np.array(bx)
            by=np.array(by)
            bz=np.array(bz)
            timedate_float=timedate_float[0:len(bx)]
            print(len(timedate_float),'Timedate len')
            print(len(bx))
            ##########################################################
            #Added to remove errorsome data
            #Note: Don't need for current Valentia setup
            #...Only looking at Variations in Valenti
        print(len(bx))
        bx[bx >= 80000.0] = 'nan'
        bx[bx == 'infs'] = 'nan'
        bx[bx <= 0.0] = 'nan'
        nans, x = nan_helper(bx)
        bx[nans]= np.interp(x(nans), x(~nans), bx[~nans])
       
        
        by[by >= 80000.0] = 'nan'
        by[by == 'infs'] = 'nan'
        by[by <= -80000.0] = 'nan'
        nans, x = nan_helper(by)
        by[nans]= np.interp(x(nans), x(~nans), by[~nans])
        
        bz[bz >= 80000.0] = 'nan'
        bz[bz == 'infs'] = 'nan'
        bz[bz <= -80000.0] = 'nan'
        nans, x = nan_helper(bz)
        bz[nans]= np.interp(x(nans), x(~nans), bz[~nans])
        ###############################################################################
        ###############################################################################
        #Creating plots
        
        print( "Starting analysis\n")
        print ("Getting data in minute bins\n")
        minute_time, minute_bx, minute_by, minute_bz= minute_bin(timedate_float, bx, by, bz, 5)
        
        
        print ("First k_index\n")
        k_index1, k_timestamp1, k_time1 = k_index(minute_time, minute_bx, minute_by, k_max)
        
        
        print ("FMI-Smoothing 1\n")
        smoothed_time, smoothed_bx, smoothed_by = fmi_smoothed2(minute_time, minute_bx,
        minute_by, k_index1, k_time1)
        
        print( "Interpolated Univariate Spline smoothing 1\n")
        smooth_time1, smooth_bx1 = smoothed(smoothed_time, smoothed_bx, 3)
        smooth_time1, smooth_by1 = smoothed(smoothed_time, smoothed_by, 3)
        
        print ("Subtracting data from smoothed\n")
        subtracted_bx1, subtracted_by1 = subtracted(minute_time, smooth_time1,
        minute_bx, smooth_bx1, minute_by, smooth_by1)
        
        print ("Second k-index\n")
        k_index2, k_timestamp2, k_time2 = k_index(minute_time, subtracted_bx1, subtracted_by1, k_max)
        
        print ("FMI-Smoothing 2\n")
        smoothed_time, smoothed_bx, smoothed_by = fmi_smoothed2(minute_time, minute_bx, 
        										minute_by, k_index2, k_time2)
        
        print ("Interpolated Univariate Spline smoothing 2\n")
        smooth_time2, smooth_bx2 = smoothed(smoothed_time, smoothed_bx, 3)
        smooth_time2, smooth_by2 = smoothed(smoothed_time, smoothed_by, 3)
        
        print ("Subtracting data from smoothed\n")
        subtracted_bx2, subtracted_by2 = subtracted(minute_time, smooth_time2, minute_bx, smooth_bx2, minute_by, smooth_by2)
        
        print ("Final k_index\n")
        k_index3, k_timestamp3, k_time3 = k_index(minute_time, subtracted_bx2, subtracted_by2, k_max)
        #Can adapt between H30 and H60 here, '30' or '60'
        mode='30'
        h_indexval, h_timestamp3, h_time3 = h_index(minute_time, subtracted_bx2, subtracted_by2, k_max,mode)
        plt.figure()
        plt.plot(h_indexval[48:])
        
        
        a = minute_time[-1] - minute_time[-1]%(24*60*60) - (2 * 24* 60 * 60)
        for i, v in enumerate(k_timestamp3):
        	if v >= a:
        		break
        k_index3 = k_index3[i:]
        k_timestamp3 = k_timestamp3[i:]
        k_time3 = k_time3[i:]
        
        k_time3 = [x - 8 for x in k_time3]
        
        print ("Plotting...")

        #creating k indices
        do_k_plots(k_index3, k_time3, k_timestamp3, minute_time,sitefull_name, save_address,save_address2)
        
        
        #manually loading Hpo indices
        #Manually editing the rows to allow for cleaner plotting
        hpo=np.loadtxt('/home/johnnyl/Documents/Solar-Tsunami/hp_index2024.txt',usecols=7,skiprows=6270-48)
        
        #hpo=np.loadtxt('/home/johnnyl/Documents/Solar-Tsunami/hp_index2024.txt',usecols=7,skiprows=13566)
        hpo=hpo[0:192]
        
        do_h_plots(hpo, h_time3, h_timestamp3, minute_time,sitefull_name, save_address,save_address2,mode)
        #do_h_plots(h_indexval, h_time3, h_timestamp3, minute_time,sitefull_name, save_address,save_address2,mode)
            
        
    #except:
    #    pass
