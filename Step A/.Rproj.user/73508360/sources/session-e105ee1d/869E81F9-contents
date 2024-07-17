#!/usr/bin/env python
# coding: utf-8

#Provided in 12/14/2022 for ATLANTIS team
#If you have any question, please contact Su Kyong Yun (sukyong.yun@pnnl.gov)
import xarray as xr
import numpy as np
import pandas as pd
import datetime

#Read the solution (one specific filename) or solutions (based on the directory)
#if you are using specific filename
filename="SSM_v43_2018_0001.nc"
ssm_solution=xr.open_dataset(filename, decode_times=False, drop_variables=['siglev','siglay'])

#if you are using directory
#mfdataDIR='../ssm_FVCOMICM_*.nc'
#ssm_solution=xr.open_mfdataset(mfdataDIR) 

node_info=pd.read_excel('ssm_node_info.xlsx',sheet_name='Info')
siglev_center=ssm_solution['siglev_center'].values[:,0]*-1
siglev_diff=np.array([siglev_center[k+1]-siglev_center[k] for k in range(0,10)])

stat_method_list=['mean','sum','min','max','25percentile','75percentile']
variable_list=['temp','salinity']

target_variable='temp' #temperature
target_stat='mean' #mean
target_data=ssm_solution[target_variable].values

#temporal_aggregation(target_data, target_stat)
#depth_average(target_data)
#spatial_aggregation(target_data)

def temporal_aggregation(target_data, target_stat):
    #aggregating based on the time
    #target_data : (time, layer, nodes)
    if target_stat=='mean':
        temp_data=target_data.mean(axis=0) #time, layer, nodes
    elif target_stat=='sum':
        temp_data=target_data.sum(axis=0) #time, layer, nodes
    elif target_stat=='min':
        temp_data=target_data.min(axis=0) #time, layer, nodes
    elif target_stat=='max':
        temp_data=target_data.max(axis=0) #time, layer, nodes
    elif target_stat=='25percentile':
        temp_data=np.percentile(target_data,25,axis=0) #time, layer, nodes
    elif target_stat=='75percentile':
        temp_data=np.percentile(target_data,75,axis=0) #time, layer, nodes
    return temp_data

def depth_average(target_data):
    #depth-average
    if len(target_data.shape)==3:
        #target_data (time, layer, nodes)
        temp_data=(np.transpose(target_data, (0,2,1))*siglev_diff).sum(axis=1)
    elif len(target_data.shape)==2:
        #target_data (layer,nodes)
        temp_data=(np.transpose(target_data, (1, 0))*siglev_diff).sum(axis=1)
    return temp_data

def spatial_aggregation(target_data):
    for unique_node in node_info['region_info'].unique():
        if len(target_data.shape)==3:
            #target_data (time, layer, nodes)
            temp_data=target_data[:,:,node_info['region_info']==unique_node].mean(axis=2)
            #temp_data (time,layer)
        elif len(target_data.shape)==2:
            #target_data (time, nodes) or target_data (layer, nodes) 
            temp_data=target_data[:,node_info['region_info']==unique_node].mean(axis=1)
            #temp_data (time) or temp_data(layer)
        elif len(target_data.shape)==1:
            #target_data (nodes)
            temp_data=target_data[node_info['region_info']==unique_node].mean(axis=1)        
            #temp_data (1)
        
    return temp_data 


# In[ ]:




