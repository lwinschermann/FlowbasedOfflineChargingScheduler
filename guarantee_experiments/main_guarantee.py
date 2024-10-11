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

#for gurobi model
from gurobipy import *
from itertools import product

from config_path import Config
import os, sys
import numpy as np
import datetime
import time
import random

# from LP import LP
from FOCS import FOCS, FlowNet, FlowOperations, FOCSinstance, Schedule
from Bookkeeping import Bookkeeping
# from AVR import AVR
# from OA import OA
from Price import Price, Greedy

print('import stuff done')

bk = Bookkeeping()
reps = 500
instanceSizes = [n for n in range(1,50)] + [5*n for n in range(10,81)] + [500, 1000, 10000]
timeSteps = [900]
maxFlows = [shortest_augmenting_path]#, edmonds_karp, preflow_push, dinitz]
write = True
randomSample = True
prefix = 'fulpes_'

'''--------- Get power profiles for 400 EVs ----------'''
bk = Bookkeeping()
'''--------------initialize instance--------------'''
# Real Training data
instanceData = pd.read_excel('data/input/filteredData_cumprob.xlsx')
instance = FOCSinstance(instanceData[:400], timeStep=900, periodicity=True, cumprob=False)

instance.add_milestones([[[23,64]] for j in instance.jobs])
instance.milestones_to_segments(toy=False)
flowNet = FlowNet()

'''--------------run FOCS--------------'''
flowNet.focs_instance_to_network(instance.parent_instance)
flowOp = FlowOperations(flowNet.G, instance.parent_instance)
focs = FOCS(instance.parent_instance, flowNet, flowOp)
f = focs.solve_focs()

# assuming all results are for the same instance and time horizon, start by saving 
bk.agg_power_profiles += [[int(focs.instance.intervals_start[i]+freq) for i in range(0, focs.instance.m) for freq in range(0,int(focs.instance.len_i[i]/focs.instance.timeStep))]]
temp = []
# determine powers
p_i = [((focs.f["i"+str(i)]['t']/(focs.instance.len_i[i]))*focs.instance.timeBase) for i in range(0,focs.instance.m)]
for i in range(0,focs.instance.m):
    # extrapolate
    temp += [p_i[i] for t in range(0,int(focs.instance.len_i[i]/focs.instance.timeStep))]
bk.agg_power_profiles += [temp]

'''--------------run FULPES 1 guarantee--------------'''
flowNet.segment_instance_to_network(instance)
flowNet.G = flowNet.G_seg
flowOp = FlowOperations(flowNet.G, instance)
focs = FOCS(instance, flowNet, flowOp)
f_fulpes1 = focs.solve_segmented_focs()
sched_fulpes_1 = Schedule(focs)

temp = []
# determine powers
p_i = [((focs.f["i"+str(i)]['t']/(focs.instance.len_i[i]))*focs.instance.timeBase) for i in range(0,focs.instance.m)]
for i in range(0,focs.instance.m):
    # extrapolate
    temp += [p_i[i] for t in range(0,int(focs.instance.len_i[i]/focs.instance.timeStep))]
bk.agg_power_profiles += [temp]

'''--------------run FULPES 2 guarantees--------------'''
instance = FOCSinstance(instanceData[:400], timeStep=900, periodicity=True, cumprob=False)
instance.add_milestones([[[15,48],[23,64]] for j in instance.jobs])
instance.milestones_to_segments(toy=False)


flowNet = FlowNet()
flowNet.segment_instance_to_network(instance)
flowNet.G = flowNet.G_seg
flowOp = FlowOperations(flowNet.G, instance)
focs = FOCS(instance, flowNet, flowOp)
f = focs.solve_segmented_focs()

temp = []
# determine powers
p_i = [((focs.f["i"+str(i)]['t']/(focs.instance.len_i[i]))*focs.instance.timeBase) for i in range(0,focs.instance.m)]
for i in range(0,focs.instance.m):
    # extrapolate
    temp += [p_i[i] for t in range(0,int(focs.instance.len_i[i]/focs.instance.timeStep))]
bk.agg_power_profiles += [temp]

'''--------------greedy--------------'''
instance = FOCSinstance(instanceData[:400], timeStep=900, periodicity=True, cumprob=False)
flowNet = FlowNet()
flowNet.focs_instance_to_network(instance)
flowOp = FlowOperations(flowNet.G, instance)

greedy = Price(instance, flowNet, flowOp, [i for i in range(0,len(instance.I_a))])
greedy.solve_cost()

temp = []
# determine powers
p_i = [((greedy.f_cost["i"+str(i)]['t']/(greedy.instance.len_i[i]))*greedy.instance.timeBase) for i in range(0,greedy.instance.m)]
for i in range(0,greedy.instance.m):
    # extrapolate
    temp += [p_i[i] for t in range(0,int(greedy.instance.len_i[i]/greedy.instance.timeStep))]
bk.agg_power_profiles += [temp]

power_profiles_df = pd.DataFrame(bk.agg_power_profiles).T
power_profiles_df.columns =  ['i', 'focs', 'fulpes1', 'fulpes2', 'greedy']
# print(power_profiles_df)

if write:
    path = 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data/output/fulpes/'
    power_profiles_df.to_csv(path + 'fulpes_power_profiles_{}.csv'.format('0_900_400'), sep=';')

time.sleep(1)

'''--------- Get running times for up to 400 EVs ----------'''
for instanceSize in instanceSizes:
    bk.empty_temp()
    for rep in range(0,reps):
        '''--------------initialize instance--------------'''
        # Sample real Training data
        sample = sorted(random.sample(range(0,len(instanceData)), instanceSize))
        instance = FOCSinstance(instanceData.iloc[sample], timeStep=900, periodicity=True, cumprob=False)

        instance.add_milestones([[[23,64]] for j in instance.jobs])
        instance.milestones_to_segments(toy=False)
        flowNet = FlowNet()

        '''--------------run FOCS--------------'''
        startFOCS = time.process_time()
        flowNet.focs_instance_to_network(instance.parent_instance)
        flowOp = FlowOperations(flowNet.G, instance.parent_instance)
        focs = FOCS(instance.parent_instance, flowNet, flowOp)
        f = focs.solve_focs()
        endFOCS = time.process_time()

        '''--------------run FULPES 1 guarantee--------------'''
        startFULPES1 = time.process_time()
        flowNet.segment_instance_to_network(instance)
        flowNet.G = flowNet.G_seg
        flowOp = FlowOperations(flowNet.G, instance)
        focs = FOCS(instance, flowNet, flowOp)
        f_fulpes1 = focs.solve_segmented_focs()
        endFULPES1 = time.process_time()

        '''--------------run FULPES 2 guarantees--------------'''
        instance = FOCSinstance(instanceData.iloc[sample], timeStep=900, periodicity=True, cumprob=False)
        instance.add_milestones([[[15,48],[23,64]] for j in instance.jobs])
        instance.milestones_to_segments(toy=False)

        flowNet = FlowNet()
        startFULPES2 = time.process_time()
        flowNet.segment_instance_to_network(instance)
        flowNet.G = flowNet.G_seg
        flowOp = FlowOperations(flowNet.G, instance)
        focs = FOCS(instance, flowNet, flowOp)
        f = focs.solve_segmented_focs()
        endFULPES2 = time.process_time()

        '''--------------save runtimes--------------'''
        bk.rtFOCStemp += [endFOCS - startFOCS]
        bk.rtFULPES1temp += [endFULPES1 - startFULPES1]
        bk.rtFULPES2temp += [endFULPES2 - startFULPES2]
        exit()
        print("End rep ", rep)
    bk.write_runtimes_full_problem()
    if write:
        bk.write_fulpes_instance_to_csv('0_900_'+str(instanceSize))
    print("Done n = ", instanceSize, "\n")
if write:
    bk.write_fulpes_to_csv(pf = '0_900', fc = instanceSizes)

#optional operations
#focs.write(type = 'pickle', suffix = '')
#instance.validate_solution(f)
#flowOp.validate_global_cap()
#focsMPC.validate_MPC()
#lp.print_results()
#focs.write()

'''plot histogram'''
# plt.figure("Histogram runtime FOCS")
# plt.hist(bk.rtFOCS[-1], bins = 100)
# plt.xlabel('Runtime [s]')
# plt.ylabel('Number of runs')
# plt.show()
# plt.close()
# '''plot histogram'''
# plt.figure("Histogram runtime LP")
# plt.hist(bk.rtLP[-1], bins = 100)
# plt.xlabel('Runtime [s]')
# plt.ylabel('Number of runs')
# plt.show()
# plt.close()
# '''plot histogram'''
# plt.figure("Histogram runtime FOCS")
# plt.hist(bk.rtFOCSmpc[-1], bins = 100)
# plt.xlabel('Runtime [s]')
# plt.ylabel('Number of runs')
# plt.show()
# plt.close()


'''
TODO
General
    synthetic instances
    change weight factors of objective function by length of intervals. Makes the values smaller too ;)
Gurobi
FOCS
    warm-start function
'''