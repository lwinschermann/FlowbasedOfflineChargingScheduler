#    Flow-based Offline Charging Scheduler (FOCS)
#    For scheduling of large groups of EVs
#    Part of SmoothEMS met GridShield project - code developed by 
#    Leoni Winschermann, University of Twente, l.winschermann@utwente.nl
#    Leander van der Bijl, University of Twente, l.c.vanderbijl@utwente.nl
#    
#    Copyright (C) 2024 CAES and MOR Groups, University of Twente, Enschede, The Netherlands
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
#    USA

import networkx as nx
import copy
from networkx.algorithms.flow import shortest_augmenting_path
from networkx.algorithms.flow import edmonds_karp
from networkx.algorithms.flow import preflow_push
from networkx.algorithms.flow import dinitz
from networkx.algorithms.flow import boykov_kolmogorov
import math
import pickle
import json
import matplotlib.pyplot as plt
import pandas as pd
import statistics
import csv
from busses import *
import json
import datetime as dt

#for gurobi model
# from gurobipy import *
from itertools import product

#from config_path import Config
import os, sys
import numpy as np
import datetime
import time
import random

# from LP import LP
from FOCS import FOCS, FlowNet, FlowOperations, FOCSinstance
from Bookkeeping import Bookkeeping
from AVR import AVR
from OA import OA
from Price import Price, Greedy



#Some variables that must be defined
maxFlows = [shortest_augmenting_path]#, edmonds_karp, preflow_push, dinitz]
savedata = True
addbaseloadreal = True
addbaseloadrandom = True
co2_source = 'ned.nl' # either 'ned.nl' or 'electricitymaps'
time_intervals = 960
relevant_interval = [97, 96*8 + 1]
createNewInstanceData = True

#max capacity per time interval in kW
capacity = 500

column_vec = ['wf', 'flat_cost', 'cotwocost', 'peak_size']
weights_list = [0.1, 0.3, 0.5, 0.7, 0.9, 1, 2, 4, 6, 8, 10]

cost_df = pd.DataFrame(columns=column_vec)
powerprofile_df = pd.DataFrame(columns = [str(i) for i in range(relevant_interval[0], relevant_interval[1])])



if createNewInstanceData is True:
    # Start creating instance
    instanceData = generateDataWeek()
    instanceData.to_csv('data/output/instance_data.csv', index=False)
else:
    instanceData = pd.read_csv('data/output/GeneratedData.csv')

if addbaseloadreal is True:
    baseload = import_baseload_data()
    instanceData = add_baseload(instanceData, baseload)
    if addbaseloadrandom is True:
        baseloadrandom = generate_baseload_random(time_intervals, 10000, 100000)
        instanceData = add_baseload(instanceData, baseloadrandom)
        if savedata is True:
                df_baseloadrandom = pd.DataFrame(baseloadrandom[96:-191])
                df_baseloadrandom.to_csv('data/output/baseloadrandom.csv', index=False)


#Add the c02 data
if co2_source == 'electricitymaps':
    co2_data = pd.read_csv("data/input/NL_2023_hourly.csv")
    co2_data['Datetime (UTC)'] = pd.to_datetime(co2_data['Datetime (UTC)'], format='%Y-%m-%d %H:%M:%S')
    co2_data = co2_data[co2_data["Datetime (UTC)"] >= dt.datetime(2023,6,1,0,0)]
    #extrapolate to quarters
    co2_data = [int(c) for c in co2_data["Carbon Intensity gCO2eq/kWh (LCA)"].to_list() for i in range(0,4)]
elif co2_source == 'ned.nl':
    co2_data = pd.read_excel("data/input/Dataset CO2 emissiefactor elektriciteit 31 mei 2023- 15 juni 2023.xlsx")
    co2_data['validfrom (UTC)'] = pd.to_datetime(co2_data['validfrom (UTC)'], format='%Y-%m-%d %H:%M:%S')
    co2_data = co2_data[co2_data['validfrom (UTC)'] >= dt.datetime(2023,6,1,0,0)]
    co2_data = [int(1000*c) for c in co2_data['emissionfactor (kg CO2/kWh)'].to_list()]




'''--------------start instance with capacity--------------'''  
instance = FOCSinstance(instanceData, timeStep=900, periodicity= True)   
#instance.add_capacity(cap = capacity)
# make sure same length as number of breakpoints. 
co2 = co2_data[instance.intervals_start[0]:instance.intervals_end[-1]]



'''--------------start Price--------------'''
flowNet = FlowNet()
flowNet.focs_instance_to_network(instance)
flowOp = FlowOperations(flowNet.G, instance)
# flowOp.validate_global_cap()

price = Price(instance, flowNet, flowOp, co2)
price.solve_cost()


f = price.flowOp.partial_flow(price.f_cost, price.instance.I_a[price.instance.I_a.index(relevant_interval[0]):price.instance.I_a.index(relevant_interval[1])], price.instance.I_a, price.instance.jobs, price.instance.J)
cost_df.loc[len(cost_df)] = [0, price.objective(price.f_cost, partial = relevant_interval), price.objective(price.f_cost, co2obj=True, partial = relevant_interval), max([f['i'+str(i)]['t']*(1/price.instance.tau) for i in price.instance.I_a])]


timelist = [0 for i in range(relevant_interval[0], relevant_interval[1])]
indx2 = 0
for indx in range(relevant_interval[0], relevant_interval[1]):
    if 'i'+str(indx) in f:
        timelist[indx2] += (f['i'+str(indx)]['t'])
        indx2 += 1
powerprofile_df.loc[0] = timelist

print('Solved for weight_f = 0')



for weight_f in weights_list:

    '''--------------start Weighted Price--------------'''
    price.solve_flatten_cost(weight_c=1, weight_f= weight_f)
    f = price.flowOp.partial_flow(price.f, price.instance.I_a[price.instance.I_a.index(relevant_interval[0]):price.instance.I_a.index(relevant_interval[1])], price.instance.I_a, price.instance.jobs, price.instance.J)

    print('Solved for weight_f = ', weight_f)

    cost_df.loc[len(cost_df)] = [weight_f, price.objective(partial = relevant_interval), price.objective(co2obj=True, partial = relevant_interval), max([f['i'+str(i)]['t']*(1/price.instance.tau) for i in price.instance.I_a])]

    timelist = [0 for i in range(relevant_interval[0], relevant_interval[1])]
    indx2 = 0
    for indx in range(relevant_interval[0], relevant_interval[1]):
        if 'i'+str(indx) in price.f:
            timelist[indx2] += (price.f['i'+str(indx)]['t'])
            indx2 += 1
    powerprofile_df.loc[len(powerprofile_df)] = timelist


'''--------------start FOCS--------------'''
focs = FOCS(instance,flowNet,flowOp)
f = focs.solve_focs()
f = focs.flowOp.partial_flow(f, focs.instance.I_a[focs.instance.I_a.index(relevant_interval[0]):focs.instance.I_a.index(relevant_interval[1])], focs.instance.I_a, focs.instance.jobs, focs.instance.J)


cost_df.loc[len(cost_df)] = ['inf', price.objective(f = focs.f, partial = relevant_interval), price.objective(f= focs.f, co2obj=True, partial = relevant_interval), max([f['i'+str(i)]['t']*(1/focs.instance.tau) for i in focs.instance.I_a])]

timelist = [0 for i in range(relevant_interval[0], relevant_interval[1])]
indx2 = 0
for indx in range(relevant_interval[0], relevant_interval[1]):
    if 'i'+str(indx) in f:
        timelist[indx2] += (f['i'+str(indx)]['t'])
        indx2 += 1
powerprofile_df.loc[len(powerprofile_df)] = timelist

print('Solved for weight_f = inf')

if savedata is True:
    cost_df.to_csv('data/output/cost_data3.csv', index=False)
    powerprofile_df.to_csv('data/output/power_profiles.csv', index=False)