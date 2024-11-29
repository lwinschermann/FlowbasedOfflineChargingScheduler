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

from LP import LP
from FOCS import FOCS, FlowNet, FlowOperations, FOCSinstance, Schedule
from Bookkeeping import Bookkeeping
from AVR import AVR
from OA import OA
from Price import Price, Greedy

bk = Bookkeeping()
reps = 500
instanceSizes = [400]
timeSteps = [900]
maxFlows = [shortest_augmenting_path]#, edmonds_karp, preflow_push, dinitz]
write = True
randomSample = True
prefix = '_mean_feas'

'''--------------initialize instance--------------'''
# Real Training data
instanceData = pd.read_excel('data/input/filteredData_cumprob.xlsx')
# filter instanceData to EVs that occur at least 50 times
instanceData = instanceData.groupby('card_id').filter(lambda x: len(x) > 50)
# print(len(instanceData['card_id'].unique()))

bk = Bookkeeping()
days = [x for x in range(0,30)]
for x in days:
    x = len(days) - 1 - x
    # choose the k^th occurance of each card_id and initialize instance
    instance = FOCSinstance(instanceData[instanceData.groupby('card_id').cumcount()==x], timeStep=900, periodicity=True, cumprob=True, prob_t1=prefix)

    '''--------------start FOCS--------------'''
    flowNet = FlowNet()
    flowNet.focs_instance_to_network(instance)
    flowOp = FlowOperations(flowNet.G, instance)

    focs = FOCS(instance,flowNet,flowOp)
    focs.solve_focs()
    schedule = Schedule(focs)

    '''naive FOCS schedule'''
    # flow to dataframe - format schedule
    s = schedule.flow_to_dataframe(focs.f)
    if write:
        s.to_csv("data/output/schedules/s{}_{}.csv".format(prefix,x), sep=';')
    schedule.derive_early_departure_schedule(f=focs.f)
    schedule.objective(focs.f)
    bk.write_objectives_schedule(schedule, weighting_method = 's', day = x, idx = instanceData[instanceData.groupby('card_id').cumcount()==x].index, n_idx = instanceData[instanceData.groupby('card_id').cumcount()==x]['count so far'])

    '''L-LYNCS'''
    schedule.solve_schedule(how = 'linear')
    # flow to dataframe - format schedule
    s_lin = schedule.flow_to_dataframe(schedule.f_sched)
    if write:
        s_lin.to_csv("data/output/schedules/s_lin{}_{}.csv".format(prefix,x), sep=';')
    schedule.derive_early_departure_schedule()
    schedule.objective(schedule.f_sched)
    bk.write_objectives_schedule(schedule, weighting_method = 's_lin', day = x, idx = instanceData[instanceData.groupby('card_id').cumcount()==x].index, n_idx = instanceData[instanceData.groupby('card_id').cumcount()==x]['count so far'])

    '''Q-LYNCS'''
    schedule.solve_schedule(how = 'quadratic')
    # flow to dataframe - format schedule
    s_q = schedule.flow_to_dataframe(schedule.f_sched)
    if write:
        s_q.to_csv("data/output/schedules/s_q{}_{}.csv".format(prefix,x), sep=';')
    schedule.derive_early_departure_schedule()
    schedule.objective(schedule.f_sched)
    bk.write_objectives_schedule(schedule, weighting_method = 's_q', day = x, idx = instanceData[instanceData.groupby('card_id').cumcount()==x].index, n_idx = instanceData[instanceData.groupby('card_id').cumcount()==x]['count so far'])

    '''L-LYNCS with weights scaled between 0 and 1'''
    schedule.solve_schedule(how = 'linear_proportional')
    # flow to dataframe - format schedule
    s_linp = schedule.flow_to_dataframe(schedule.f_sched)
    if write:
        s_linp.to_csv("data/output/schedules/s_linp{}_{}.csv".format(prefix,x), sep=';')
    schedule.derive_early_departure_schedule()
    schedule.objective(schedule.f_sched)
    bk.write_objectives_schedule(schedule, weighting_method = 's_linp', day = x, idx = instanceData[instanceData.groupby('card_id').cumcount()==x].index, n_idx = instanceData[instanceData.groupby('card_id').cumcount()==x]['count so far'])

    '''Q-LYNCS with weights scaled between 0 and 1'''
    schedule.solve_schedule(how = 'quadratic_proportional')
    # flow to dataframe - format schedule
    s_qp = schedule.flow_to_dataframe(schedule.f_sched)
    if write:
        s_qp.to_csv("data/output/schedules/s_qp{}_{}.csv".format(prefix,x), sep=';')
    schedule.derive_early_departure_schedule()
    schedule.objective(schedule.f_sched)
    bk.write_objectives_schedule(schedule, weighting_method = 's_qp', day = x, idx = instanceData[instanceData.groupby('card_id').cumcount()==x].index, n_idx = instanceData[instanceData.groupby('card_id').cumcount()==x]['count so far'])

    '''H-LYNCS'''
    schedule.solve_schedule(how = 'probabilistic')
    # flow to dataframe - format schedule
    s_prob = schedule.flow_to_dataframe(schedule.f_sched)
    if write:
        s_prob.to_csv("data/output/schedules/s_prob{}_{}.csv".format(prefix,x), sep=';')
    schedule.derive_early_departure_schedule()
    schedule.objective(schedule.f_sched)
    bk.write_objectives_schedule(schedule, weighting_method = 's_prob', day = x, idx = instanceData[instanceData.groupby('card_id').cumcount()==x].index, n_idx = instanceData[instanceData.groupby('card_id').cumcount()==x]['count so far'])

    '''H-LYNCS without delay till after z_j'''
    schedule.solve_schedule(how = 'probabilistic_full')
    # flow to dataframe - format schedule
    s_probfull = schedule.flow_to_dataframe(schedule.f_sched)
    if write:
        s_probfull.to_csv("data/output/schedules/s_probfull{}_{}.csv".format(prefix,x), sep=';')
    schedule.derive_early_departure_schedule()
    schedule.objective(schedule.f_sched)
    bk.write_objectives_schedule(schedule, weighting_method = 's_probfull', day = x, idx = instanceData[instanceData.groupby('card_id').cumcount()==x].index, n_idx = instanceData[instanceData.groupby('card_id').cumcount()==x]['count so far'])
# bk save all results
bk.write_objectives_schedule_to_csv(data = instanceData, weighting_methods=['s', 's_lin', 's_linp', 's_q', 's_qp', 's_prob', 's_probfull'], days = days, jobs = instanceData[instanceData.groupby('card_id').cumcount()==0]['card_id'])





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