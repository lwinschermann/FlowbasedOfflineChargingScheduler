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

from LP import LP
from FOCS import FOCS, FlowNet, FlowOperations, FOCSinstance, Bookkeeping

bk = Bookkeeping()
reps = 5
instanceSizes = [n for n in range(1,20)] + [n for n in range(20,501,10)] + [n for n in range(600, 1001, 100)] 
timeSteps = [3600]
maxFlows = [shortest_augmenting_path]
write = True

for maxFlowIndex, maxFlowAlg in enumerate(maxFlows): # algorithm used to determine max flow 
    for timeStep in timeSteps: # set interval length. 900 = 15 min intervals, 1 = second-based
        for instanceSize in instanceSizes: # set number of jobs. Sample the last x sessions
            bk.prefix += [str(maxFlowIndex) + "_" + str(timeStep) + "_" + str(instanceSize) + "_"]

            bk.empty_temp()
            bk.read_temp('data/output/instances/instance_{}_{}_{}_full.csv'.format(maxFlowIndex, timeStep, instanceSize))
            bk.write_runtimes_full_problem()

            bk.empty_temp()
            bk.read_temp('data/output/instances/instance_{}_{}_{}_partial.csv'.format(maxFlowIndex, timeStep, instanceSize))
            bk.write_runtimes_partial_problem()
        bk.write_timestep_to_csv(fc=instanceSizes)
    # bk.write_flowmethod_to_csv(fc = instanceSizes, sc = timeSteps)    
