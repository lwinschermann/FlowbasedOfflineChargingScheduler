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
import pandas as pd
import statistics
import csv

#for gurobi model
from gurobipy import *
from itertools import product

import os, sys
import datetime
import time
import random

from LP import LP
from FOCS import FOCS, FlowNet, FlowOperations, FOCSinstance, Bookkeeping

bk = Bookkeeping()
idle = 0.01
reps = 500
instanceSizes = [n for n in range(1,20)] + [n for n in range(20,501,10)] + [n for n in range(600, 1001, 100)] 
timeSteps = [60]
maxFlows = [preflow_push]
write = True
randomSample = True

# Real Training data
instanceData = pd.read_excel('data/input/filteredData.xlsx')

timeTotal = time.process_time()
for maxFlowIndex, maxFlowAlg in enumerate(maxFlows): # algorithm used to determine max flow 
    maxFlowIndex = 2
    for timeStep in timeSteps: # set interval length. 900 = 15 min intervals, 1 = second-based
        for instanceSize in instanceSizes: # set number of jobs. Sample the last x sessions
            bk.prefix += [str(maxFlowIndex) + "_" + str(timeStep) + "_" + str(instanceSize) + "_"]
            if not randomSample:
                instance = FOCSinstance(instanceData[:instanceSize], timeStep)
            bk.empty_temp()
            for rep in range(0,reps):
                if randomSample:
                    sample = sorted(random.sample(range(0,len(instanceData)), instanceSize))
                    instance = FOCSinstance(instanceData.iloc[sample], timeStep)  

                '''--------------start LP formulation--------------'''
                time.sleep(idle)
                startLP = time.process_time()
                lp = LP(instance)
                lp.build_model()
                bk.mbLPtemp += [time.process_time() - startLP]
                lp.solve_model()
                bk.rtLPtemp += [time.process_time() - startLP]
                bk.solLPtemp += [bk.rtLPtemp[-1] - bk.mbLPtemp[-1]]

                '''--------------start FOCS--------------'''
                time.sleep(idle)
                startFOCS = time.process_time()
                flowNet = FlowNet()
                flowNet.focs_instance_to_network(instance)
                flowOp = FlowOperations(flowNet.G, instance)
                focs = FOCS(instance, flowNet, flowOp)
                focs.flow_func = maxFlowAlg
                bk.mbFOCStemp += [time.process_time() - startFOCS]

                f = focs.solve_focs(MPCstopper=False, MPCcondition=0)
                bk.rtFOCStemp += [time.process_time() - startFOCS]
                bk.solFOCStemp += [bk.rtFOCStemp[-1] - bk.mbFOCStemp[-1]]

                print('Done FOCS full day version')

                '''--------------start FOCS with pre-mature stop --------------'''
                time.sleep(idle)
                startFOCSmpc = time.process_time()
                flowNet = FlowNet()
                flowNet.focs_instance_to_network(instance)
                flowOp = FlowOperations(flowNet.G, instance)
                focsMPC = FOCS(instance, flowNet, flowOp)
                focsMPC.flow_func = maxFlowAlg
                bk.mbFOCSmpctemp += [time.process_time() - startFOCSmpc]

                fMPC = focsMPC.solve_focs(MPCstopper=True, MPCcondition=0)
                bk.rtFOCSmpctemp += [time.process_time() - startFOCSmpc]
                bk.solFOCSmpctemp += [bk.rtFOCSmpctemp[-1] - bk.mbFOCSmpctemp[-1]]
                print('Done FOCS pre-mature stop full day version')

                ''' ----------------same procdure with problem starting at noon -----------------'''
                # add this check to avoid empty input
                if 12/instance.tau < instance.breakpoints[-1]:
                    instance.reduced_problem(f, start_h = 12)

                    '''--------------start LP formulation--------------'''
                    time.sleep(idle)
                    startLP = time.process_time()
                    lp = LP(instance)
                    lp.build_model()
                    bk.mbLPtempred += [time.process_time() - startLP]

                    lp.solve_model()
                    bk.rtLPtempred += [time.process_time() - startLP]
                    bk.solLPtempred += [bk.rtLPtempred[-1] - bk.mbLPtempred[-1]]

                    '''--------------start FOCS--------------'''
                    time.sleep(idle)
                    startFOCS = time.process_time()
                    flowNet = FlowNet()
                    flowNet.focs_instance_to_network(instance)
                    flowOp = FlowOperations(flowNet.G, instance)
                    focs = FOCS(instance, flowNet, flowOp)
                    focs.flow_func = maxFlowAlg
                    bk.mbFOCStempred += [time.process_time() - startFOCS]

                    f = focs.solve_focs(MPCstopper=False, MPCcondition=0)
                    bk.rtFOCStempred += [time.process_time() - startFOCS]
                    bk.solFOCStempred += [bk.rtFOCStempred[-1] - bk.mbFOCStempred[-1]]

                    '''--------------start FOCS with pre-mature stop --------------'''
                    time.sleep(idle)
                    startFOCSmpc = time.process_time()
                    flowNet = FlowNet()
                    flowNet.focs_instance_to_network(instance)
                    flowOp = FlowOperations(flowNet.G, instance)
                    focsMPC = FOCS(instance, flowNet, flowOp)
                    focsMPC.flow_func = maxFlowAlg
                    bk.mbFOCSmpctempred += [time.process_time() - startFOCSmpc]

                    fMPC = focsMPC.solve_focs(MPCstopper=True, MPCcondition=0)
                    bk.rtFOCSmpctempred += [time.process_time() - startFOCSmpc]
                    bk.solFOCSmpctempred += [bk.rtFOCSmpctempred[-1] - bk.mbFOCSmpctempred[-1]]      
                else:
                    pass      

                print('finish rep ', rep)    

            if write:
                bk.write_runtimes_full_problem(focs, lp, focsMPC)
                bk.write_instance_to_csv(pf = 'full')
                bk.write_runtimes_partial_problem(focs, lp, focsMPC)                
                bk.write_instance_to_csv(pf = 'partial')

            print('Done ' + bk.prefix[-1])
        if write:
            bk.write_timestep_to_csv(fc = instanceSizes)
    if write:
        bk.write_flowmethod_to_csv(fc = instanceSizes, sc = timeSteps)    
print(time.process_time() - timeTotal)

#optional operations
#focs.write(type = 'pickle', suffix = '')
#instance.validate_solution(f)
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
Gurobi
FOCS
    warm-start function
'''