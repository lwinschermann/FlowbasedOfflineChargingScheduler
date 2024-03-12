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
import math
import matplotlib.pyplot as plt
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
from FOCS import FOCS, FlowNet, FlowOperations, FOCSinstance
from Bookkeeping import Bookkeeping

instanceSize = 10 #number of EVs/jobs in instance
timeStep = 900 #quarterly granularity
maxFlowAlg = shortest_augmenting_path #alternatively use e.g., edmonds_karp, preflow_push, or dinitz
randomSample = True

# Real Training data
instanceData = pd.read_excel('data/input/filteredData.xlsx')

if not randomSample:
    instance = FOCSinstance(instanceData[:instanceSize], timeStep)
if randomSample:
    sample = sorted(random.sample(range(0,len(instanceData)), instanceSize))
    instance = FOCSinstance(instanceData.iloc[sample], timeStep)  

'''--------------start FOCS--------------'''
flowNet = FlowNet()
flowNet.focs_instance_to_network(instance)
flowOp = FlowOperations(flowNet.G, instance)
focs = FOCS(instance, flowNet, flowOp)
focs.flow_func = shortest_augmenting_path
f = focs.solve_focs(MPCstopper=False, MPCcondition=0)

obj_val = focs.objective()
print('FOCS objective value = ', obj_val)
print('FOCS flow (schedule in middle edge layer): \n', focs.f)

#optional operations
#instance.toy_instance_2() #overwrites (potentially empty) FOCS instance with toy instance
#instance.reduced_problem(f, start_h = 12) #creates instance remaining at noon
#instance.validate_solution(f)
#flowOp.validate_global_cap()
#focsMPC.validate_MPC()
#lp.print_results()
#focs.write()


'''
TODO
General
    synthetic instances
Gurobi
FOCS
    warm-start function
'''