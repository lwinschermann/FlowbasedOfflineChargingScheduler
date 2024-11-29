#    Flow-based Offline Charging Scheduler (FOCS)
#    For scheduling of large groups of EVs
#    Part of SmoothEMS met GridShield project - code developed by 
#    Leoni Winschermann, University of Twente, l.winschermann@utwente.nl
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
co2_source = 'ned.nl' # either 'ned.nl' or 'electricitymaps'
time_intervals = 960
# relevant_interval = [97, 96*8 + 1]

column_vec = ['wf', 'flat_cost', 'cotwocost', 'peak_size']
weights_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# Real Training data
instanceData = pd.read_excel('data/input/filteredData.xlsx')

instance = FOCSinstance(instanceData[:400], timeStep=900, periodicity=True)
#Add the co2 data
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

co2 = co2_data[instance.intervals_start[0]:instance.intervals_end[-1]]

def determine_power_profile(f, alg):
    # determine powers
    p_i = [((f["i"+str(i)]['t']/(alg.instance.len_i[i]))*alg.instance.timeBase) for i in range(0,alg.instance.m)]
    profile = []
    for i in range(0,alg.instance.m):
        # extrapolate
        profile += [p_i[i] for t in range(0,int(alg.instance.len_i[i]/alg.instance.timeStep))]
    return [0] + profile + [0]

'''--------------start FOCS--------------'''
flowNet = FlowNet()
flowNet.focs_instance_to_network(instance)
flowOp = FlowOperations(flowNet.G, instance)

focs = FOCS(instance,flowNet,flowOp)

'''--------------start Bookkeeping--------------'''
# assuming all results are for the same instance and time horizon
i_timeStep = [focs.instance.intervals_start[i]+freq for i in range(0, focs.instance.m) for freq in range(0,int(focs.instance.len_i[i]/focs.instance.timeStep))]
power_profiles = [[min(i_timeStep)-1]+ i_timeStep + [max(i_timeStep)+1]]
metrics = []
header = ["i"]

'''--------------start Uncontrolled----------------'''
greedy = Price(instance, flowNet, flowOp, [i for i in range(0,len(co2))])
greedy.solve_cost()

power_profiles += [determine_power_profile(greedy.f_cost, greedy)]

metrics += [[greedy.objective(greedy.f_cost), sum([greedy.f_cost['i'+str(i)]['t']*co2[i] for i in focs.instance.I_a]), max([greedy.f_cost['i'+str(i)]['t']*(1/greedy.instance.tau) for i in greedy.instance.I_a])]]
# metrics += [[greedy.objective(greedy.f_cost), sum([co2[i]+power_profiles[-1][i+1]*greedy.instance.tau for i in range(0,len(co2))]), max([greedy.f_cost['i'+str(i)]['t']*(1/greedy.instance.tau) for i in greedy.instance.I_a])]]
header += ["uncontrolled"]

'''--------------start Price--------------'''
price = Price(instance, flowNet, flowOp, co2)
price.solve_cost()

power_profiles += [determine_power_profile(price.f_cost, price)]
metrics += [[price.objective(price.f_cost), price.objective(price.f_cost, co2obj=True), max([price.f_cost['i'+str(i)]['t']*(1/price.instance.tau) for i in price.instance.I_a])]]
header += ["price"]

'''--------------start Weighted Price--------------'''
for w_f in weights_list:
    price.solve_flatten_cost(weight_c=1, weight_f=w_f)

    power_profiles += [determine_power_profile(price.f, price)]
    metrics += [[price.objective(), price.objective(co2obj=True), max([price.f['i'+str(i)]['t']*(1/price.instance.tau) for i in price.instance.I_a])]]
    header += ["weighted_"+str(w_f)]

'''--------------start FOCS--------------'''
flowNet = FlowNet()
flowNet.focs_instance_to_network(instance)
flowOp = FlowOperations(flowNet.G, instance)

focs = FOCS(instance,flowNet,flowOp)
focs.solve_focs()

power_profiles += [determine_power_profile(focs.f, focs)]
metrics += [[focs.objective(), price.objective(f= focs.f, co2obj=True), max([focs.f['i'+str(i)]['t']*(1/focs.instance.tau) for i in focs.instance.I_a])]]
header += ['focs']

'''---------------save results-----------------'''

# make dataframe
power_profiles_df = pd.DataFrame(power_profiles).T
power_profiles_df.columns = header 

metrics_df = pd.DataFrame(metrics)
metrics_df['mode'] = header[1:]
metrics_df['weights'] = [None, None] + weights_list + [None]
metrics_df.columns = ["flatness", "co2", "peak", "mode", "weights"]
metrics_df = metrics_df[["mode", "weights", "flatness", "co2", "peak"]]

emissions_df = pd.DataFrame(co2_data)

# write to csv
path = 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data/output/co2objective/'
print(metrics_df)
power_profiles_df.to_csv(path + "power_profiles_0_900_400.csv", sep = ";", header = True)
metrics_df.to_csv(path + "metrics_0_900_400.csv", sep = ";", header = True)
metrics_df.dropna(subset=['weights'], inplace = True)
metrics_df.to_csv(path + "metrics_weightedonly_0_900_400.csv", sep = ";", header = True)
emissions_df.to_csv(path + "co2.csv", sep = ";", header = True)
