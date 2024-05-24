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
import math
import pickle
import json
import statistics
import pandas as pd
import csv
import numpy as np   

from FOCS import FlowNet, FlowOperations, FOCS

class OA:
    def __init__(self, instance, flowNet, flowOp):
        self.instance = instance
        self.flowNet = flowNet
        self.flowOp = flowOp
        self.f = flowOp.empty_flow(flowOp.G)
        self.J_a = [] #active jobs. Those available at any current time

    #function takes self.instance and job indices, and adds them to the object instance.partial_instance which in turn is a FOCS instance object
    def add_jobs(self, new_job_indices, start = 0):
        #initialize jobs 
        #BUG there is no check to see whether jobs are already in instance!!!
        self.instance.partial_instance.jobs += [self.instance.jobs[j] for j in new_job_indices] 
        self.instance.partial_instance.jobs_cap += [self.instance.jobs_cap[j] for j in new_job_indices] # in W
        self.instance.partial_instance.jobs_demand += [self.instance.jobs_demand[j] for j in new_job_indices]  #in kWh
        self.instance.partial_instance.n = len(self.instance.partial_instance.jobs) #number of jobs

        #initialize intervals
        #FIXME !!!!!! check how done in reduced problem!!
        self.instance.partial_instance.breakpoints = self.instance.partial_instance.breakpoints
        self.instance.partial_instance.intervals_start = self.instance.partial_instance.breakpoints[:-1]
        self.instance.partial_instance.intervals_end = self.instance.partial_instance.breakpoints[1:]
        self.instance.partial_instance.m = len(self.instance.partial_instance.intervals_start) #number of atomic intervals
        self.instance.partial_instance.I_a = self.active_intervals(jobs = self.instance.partial_instance.jobs, J_inverse = self.instance.J_inverse, start = start)
        self.instance.partial_instance.len_i = self.instance.len_i#[self.instance.len_i[i] for i in self.instance.partial_instance.I_a]#[(self.instance.partial_instance.intervals_end[i] - self.instance.partial_instance.intervals_start[i])*self.instance.partial_instance.timeStep for i in self.instance.partial_instance.I_a]

        self.instance.find_partial_J_inverse()
        self.instance.find_partial_J()
        return
    
    def active_intervals(self, jobs, J_inverse, start = 0):
        temp = []
        for j in jobs:
            temp += J_inverse['j'+str(j)]
        temp = [i for i in sorted(list(dict.fromkeys(temp))) if i >= start]
        return temp
    
    def solve_oa(self):
        #initialize empty partial instance
        self.instance.empty_instance()
        self.instance.backup = copy.deepcopy(self.instance)
        self.instance.partial_instance.breakpoints = self.instance.breakpoints
        jobs_arrived = []

        #start MPC loop. Triggered at breakpoints
        for tickIndex, tick in enumerate(self.instance.intervals_start):
            jobs_new = [j for j in self.instance.J['i'+str(tickIndex)] if j not in jobs_arrived]
            jobs_arrived += jobs_new

            #update instance with new arrivals
            self.add_jobs(jobs_new, start = tickIndex)

            #make network
            flowNet = FlowNet()
            flowNet.focs_instance_to_network(self.instance.partial_instance)
            flowOp = FlowOperations(flowNet.G, self.instance.partial_instance)
            G_tick = copy.deepcopy(flowNet.G) #focs.solve will alter network. We need it to add next time step to overall solution

            #initialize FOCS
            focs = FOCS(self.instance.partial_instance, flowNet, flowOp)

            #determine flow. If instance is empty, use empty flow
            if len(self.instance.partial_instance.jobs) == 0:
                #save 0-flow to f
                f_tick = focs.f
            else:
                #solve FOCS
                f_tick = focs.solve_focs(MPCstopper=False, MPCcondition=0)

            #isolate first active interval flow flow_dict_tick
            flow_dict_tick = self.flowOp.partial_flow(f_tick, self.instance.partial_instance.I_a[:1], self.instance.partial_instance.I_a, self.instance.partial_instance.jobs, self.instance.partial_instance.J)

            #save corresponding flow flow_dict_tick
            self.f = self.flowOp.add_flows(self.f,flow_dict_tick,G_tick)

            #reduce problem instance by corresponding flow
            self.instance.partial_instance.reduced_problem_mpc(f_tick, start_idx = tickIndex)

        return self.f
    
    #copied from FOCS class
    def objective(self):
        #normalized objective for validation
        #instead of the sum of squares, we first weight the square of the power (e/len_i[i]) by len_i[i]/timeStep      
        #determine power per interval
        p_i = [((self.f["i"+str(i)]['t']/(self.instance.len_i[i]))*self.instance.timeBase) for i in range(0,self.instance.m)]
        #determine weighted squared powers
        powerSquare = [(p_i[i]**2)*(self.instance.len_i[i]/self.instance.timeStep) for i in range(0,self.instance.m)]
        self.objNormalized = sum(powerSquare)
        return self.objNormalized