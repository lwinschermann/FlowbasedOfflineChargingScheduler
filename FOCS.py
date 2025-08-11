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
import random
import matplotlib.pyplot as plt

class FOCSinstance:
    def __init__(self, instanceData, timeStep, timeBase = 3600, periodicity = False, cumprob=False, prob_t0='', prob_t1=''):
        self.data = instanceData
        self.timeStep = timeStep
        self.timeBase = timeBase
        self.tau = timeStep/timeBase # conversion factor power to energy for one time step
        self.global_cap_active = False
        self.periodicity = periodicity
        self.cumprob = cumprob
        self.prob_t0 = prob_t0
        self.prob_t1 = prob_t1
        self.milestones_active = False
        self.milestones = None

        #initialize jobs
        self.n = len(instanceData) #number of jobs
        self.jobs = [j for j in range(0,self.n)] 
        self.jobs_cap = [22 if instanceData["average_power_W"].iloc[j]>instanceData["maxPower"].iloc[j] else instanceData["maxPower"].iloc[j]/1000 for j in self.jobs] # in W
        self.jobs_demand = [instanceData["total_energy_Wh"].iloc[j]/1000 for j in self.jobs]  #in kWh
        self.jobs_departure = [instanceData["t1_" + str(timeStep)+str(prob_t1)].iloc[j] for j in self.jobs]
        self.jobs_arrival = [instanceData["t0_" + str(timeStep)+str(prob_t0)].iloc[j] for j in self.jobs]

        #initialize intervals
        self.breakpoints = sorted(instanceData["t0_"+str(timeStep)]._append(instanceData["t1_" + str(timeStep)]).drop_duplicates().tolist())
        if self.periodicity:
            self.augment_breakpoints()
            if self.cumprob:
                # ! only in data for 900s timeSteps. In this setup, only works for simulation of 1 day.
                self.jobs_cumprob = [[instanceData['cumprob_i'+str(i+1)].iloc[j] for i in range(self.breakpoints[0],self.breakpoints[-1])] for j in self.jobs]
        self.intervals_start = self.breakpoints[:-1]
        self.intervals_end = self.breakpoints[1:]
        self.m = len(self.intervals_start) #number of atomic intervals
        self.I_a = [i for i in range(0,self.m)] 
        self.len_i = [(self.intervals_end[i] - self.intervals_start[i])*timeStep for i in self.I_a]

        self.find_J_inverse()
        self.find_J()

    def add_milestones(self, milestones = None, activate = True):
        self.milestones_active = activate
        if milestones is None:
            # add asr milestones: 23kWh by 4pm (if quarters then by quarter 64)
            self.milestones = [[[23,64]] for j in self.jobs] # if multiple guarantees, energy is additive!
        else:
            # add milestones
            self.milestones = milestones
            if len(self.milestones) != self.n:
                print('[WARNING]: milestones should have format list with #jobs entries and each jobs entry a list of lists of the form [energy(aggregated), deadline].')

    def milestones_to_segments(self, toy=False):
        if self.milestones is None:
            print('[WARNING]: You are trying to use segmentation, but have not defined milestones...')
        self.segments = []
        for j in self.jobs:
            segments_j = self.milestones[j]

            # if an EV departs way before the milestone, we adapt the segment to charge as much as HAS to be charged before departure, so that it would have been possible to charge the entire guarantee by the milestone.
            # NOTE Design choice. May as well have chosen to charge to as close as possible to the guarantee (full power till departure)
            feas_j = [max(0,stone[0] - max(0,stone[1]-self.jobs_departure[j])*self.jobs_cap[j]*self.tau ) for stone in segments_j]

            # ensure milestone feasibility
            # NOTE we also capped milestones at the demand of the vehicle, and at what is feasibly possible before the milestone deadline
            segments_j = [[min(stone[0], feas_j[id], self.jobs_demand[j], max(0,(stone[1]-self.jobs_arrival[j])*self.jobs_cap[j]*self.tau)), stone[1]] for id,stone in enumerate(segments_j)]

            # check if we need to add a segment with the actual departure
            if max([0]+[segment[0] for segment in segments_j if segment[1] <= self.jobs_departure[j]]) < self.jobs_demand[j]:
                segments_j += [[self.jobs_demand[j],self.jobs_departure[j]]]
            
            # sort based on deadline
            segments_j = sorted(segments_j, key=lambda x: x[1])

            # make non-additive
            temp_segments_j = copy.deepcopy(segments_j) # or else overwrite the previous entries already before fixing the next
            for id in range(1,len(segments_j)):
                segments_j[id][0] = max(0,temp_segments_j[id][0] - temp_segments_j[id-1][0])
            
            # drop zero-demand jobs. Else get errors in heuristic 1.
            segments_j = [segment for segment in segments_j if segment[0]>0]

            self.segments += [segments_j]
        if not toy:
            self.find_seg_J_inverse()
            self.find_seg_J()
        
            # transform milestones into segments
            self.seg_n = [len(self.segments[j]) for j in self.jobs] #number of segments
            self.seg_jobs = [[id for id,k in enumerate(self.segments[j])] for j in self.jobs]
            self.seg_jobs_demand = [[milestone[0] for milestone in self.segments[j]] for j in self.jobs]
        
        # do segment adaptations to instance
        self.parent_instance = copy.deepcopy(self)
        self.parent_jobs = [j for j in self.jobs for k in self.seg_jobs[j]]
        self.jobs = [j for j in range(0,len(self.parent_jobs))]
        self.jobs_demand = [dem for j in self.parent_instance.jobs for dem in self.seg_jobs_demand[j]]
        self.jobs_cap = [self.jobs_cap[j] for j in self.parent_instance.jobs for k in self.seg_jobs_demand[j]]
        
        self.J_inverse = {}
        for j in self.jobs:
            self.J_inverse["j"+str(j)] = self.parent_instance.seg_J_inverse["j"+str(self.parent_jobs[j])+str("-"+str(self.parent_jobs[:j].count(self.parent_jobs[j])))] 
        self.find_J()
        
        return 

    def find_seg_J_inverse(self):
        # example: self.seg_J_inverse = {"j0-1": [0,1], "j0-1": [0,1,2], "j1-0": [0], "j1-1": [0,1,2]}
        self.seg_J_inverse = {}
        for j in self.jobs:
            for idx, segment in enumerate(self.segments[j]):
                self.seg_J_inverse["j"+str(j)+"-"+str(idx)] = [i for i in self.J_inverse["j"+str(j)] if i < self.breakpoints.index(segment[1])]
        return
    
    def find_seg_J(self):
        self.seg_J = {}
        for i in self.I_a:
            self.seg_J["i"+str(i)] = {}
            for j in self.J["i"+str(i)]:
                self.seg_J["i"+str(i)][j] = [id for id,k in enumerate(self.segments[j]) if i in self.seg_J_inverse["j"+str(j)+"-"+str(id)]]
        return


    def augment_breakpoints(self):
        #number of units of timeStep from first to last breakpoint
        if len(self.breakpoints) <= 0:
            print('[WARNING]: will not augment empty breakpoints. Assumed to be empty instance.')
            return
        nr = int(max(self.breakpoints) - min(self.breakpoints))
        if abs(nr - max(self.breakpoints) + min(self.breakpoints)) > 0.0001:
            print('[WARNING]: Breakpoints non-integer! Check augment_breakpoints() function in FOCSinstance class.')
        self.breakpoints = [int(min(self.breakpoints) + bp) for bp in range(0,nr + 1)]

    def empty_instance(self):
        self.partial_instance = FOCSinstance(self.data, self.timeStep, self.timeBase)
        self.partial_instance.jobs = []
        self.partial_instance.n = 0
        self.partial_instance.jobs_cap = []
        self.partial_instance.jobs_demand = []

        self.partial_instance.breakpoints = []
        self.partial_instance.intervals_start = []
        self.partial_instance.intervals_end = []
        self.partial_instance.m = 0
        self.partial_instance.I_a = [] 
        self.partial_instance.len_i = []

        self.partial_instance.J_inverse = {}
        self.partial_instance.J = {}

    #based on a flow, determine the feasible 
    def feasible_subinstance(self, f, I):
        # update jobs
        self.partial_instance.jobs = self.jobs
        self.partial_instance.n = self.n
        self.partial_instance.jobs_cap = self.jobs_cap
        self.partial_instance.jobs_demand = [sum([f['j'+str(j)]['i'+str(i)] for i in self.J_inverse['j'+str(j)] if i in I]) for j in self.jobs]

        # update intervals
        self.partial_instance.breakpoints = self.breakpoints
        self.partial_instance.intervals_start = self.breakpoints[:-1]
        self.partial_instance.intervals_end = self.breakpoints[1:]
        self.partial_instance.m = len(self.intervals_start) 
        self.partial_instance.I_a = I
        self.partial_instance.len_i = self.len_i

        # update J and J_inverse
        self.find_partial_J_inverse()
        self.find_partial_J()

        return self.partial_instance

    def reduced_problem(self, flow, start_h, err = 0.0000001):
        #choose breakpoint closest to start
        start = start_h/self.tau
        iStartBp = min(self.breakpoints, key= lambda x:abs(x-start))
        iStart = self.breakpoints.index(iStartBp)

        #demand for all later than iStart
        supplied = [sum([flow['j'+str(j)]['i'+str(i)] for i in self.J_inverse['j'+str(j)] if i< iStart]) for j in self.jobs]

        #correct demand
        #if bigger --> start new starting point, add some demand
        if iStartBp > start:
            fraction = (iStartBp - start)/(self.len_i[iStart - 1]/self.timeStep)
            supplied = [supplied[j] - (flow['j'+str(j)]['i'+str(iStart-1)])*fraction if iStart-1 in self.J_inverse['j'+str(j)] else supplied[j] for j in self.jobs]
            self.breakpoints = sorted(self.breakpoints + [start])

        #if smaller --> start later. Remove some additional demand. Start new starting point.
        elif iStartBp < start:
            #FIXME add check if there is any later breakpoints.
            fraction = (start - iStartBp)/(self.len_i[iStart]/self.timeStep)
            supplied = [supplied[j] + (flow['j'+str(j)]['i'+str(iStart)]*fraction) if iStart in self.J_inverse['j'+str(j)] else supplied[j] for j in self.jobs]
            self.breakpoints = sorted(self.breakpoints + [start])

        self.jobs_demand = [self.jobs_demand[j] - supplied[j] for j in range(0,len(self.jobs))]

        # drop jobs that are completed.
        # jobMask = [(self.jobs_demand[j] > err) for j in range(0,len(self.jobs))]

        # # update jobs
        # self.jobs_cap = [self.jobs_cap[j] for j in range(0,len(self.jobs)) if jobMask[j]]
        # self.jobs_demand = [self.jobs_demand[j] for j in range(0,len(self.jobs)) if jobMask[j]]
        # self.jobs = [self.jobs[j] for j in range(0,len(self.jobs)) if jobMask[j]]
        # self.n = len(self.jobs)

        # update intervals
        self.breakpoints = [self.breakpoints[i] for i in range(0,len(self.breakpoints)) if self.breakpoints[i] >= start]
        self.intervals_start = self.breakpoints[:-1]
        self.intervals_end = self.breakpoints[1:]
        self.m = len(self.intervals_start) #number of atomic intervals
        self.I_a = [i for i in range(0,self.m)] 
        self.len_i = [(self.intervals_end[i] - self.intervals_start[i])*self.timeStep for i in self.I_a]
        self.find_J_inverse(startCon = True)
        self.find_J()
        return self

    def reduced_problem_mpc(self, flow, start_idx, err = 0.0000001):
        # demand for all after this timestep
        supplied = [sum([flow['j'+str(j)]['i'+str(i)] for i in self.J_inverse['j'+str(j)] if i<=start_idx]) for j in self.jobs]
        self.jobs_demand = [self.jobs_demand[j] - supplied[j] for j in range(0,len(self.jobs))]

        # drop jobs that are completed.
        jobMask = [(self.jobs_demand[j] > err) for j in range(0,len(self.jobs))]

        # update jobs
        self.jobs_cap = [self.jobs_cap[j] for j in range(0,len(self.jobs)) if jobMask[j]]
        self.jobs_demand = [self.jobs_demand[j] for j in range(0,len(self.jobs)) if jobMask[j]]
        self.jobs = [self.jobs[j] for j in range(0,len(self.jobs)) if jobMask[j]]
        self.n = len(self.jobs)

        # update intervals
        self.breakpoints = [self.breakpoints[i] for i in range(0,len(self.breakpoints)) if self.breakpoints[i] >= start_idx]
        self.intervals_start = self.breakpoints[:-1]
        self.intervals_end = self.breakpoints[1:]
        self.m = len(self.intervals_start) #number of atomic intervals
        self.I_a = [i+start_idx for i in range(0,self.m)] 

        return self

    def find_J_inverse(self, startCon = False):
        self.J_inverse = {}
        if not startCon:
            for j in self.jobs:
                self.J_inverse["j"+str(j)] = [i for i in range(self.intervals_start.index(self.data["t0_"+str(self.timeStep)].iloc[j]), self.intervals_end.index(self.data["t1_"+str(self.timeStep)].iloc[j])+1)] 
        else:
            for j in self.jobs:
                try:
                    temp = self.intervals_start.index(self.data["t0_"+str(self.timeStep)+self.prob_t0].iloc[j])
                except:
                    if self.data["t0_"+str(self.timeStep)+self.prob_t0].iloc[j] < min(self.intervals_start):
                        temp = 0
                    else: 
                        print('[WARNING]: Unconsidered case. Check find_J_inverse() function in FOCSinstance class.')
                self.J_inverse["j"+str(j)] = [i for i in range(temp, self.intervals_end.index(self.data["t1_"+str(self.timeStep)+self.prob_t1].iloc[j])+1)] 
        return

    def find_J(self):
        self.J = {}
        for i in self.I_a:
            self.J["i"+str(i)] = [j for j in self.jobs if i in self.J_inverse["j"+str(j)]]
        return

    def find_partial_J_inverse(self):
        self.partial_instance.J_inverse = {}
        #take only those jobs that are in self.instance.partial_instance.jobs
        for j in self.partial_instance.jobs:
            #take self.instance.J_inverse and remove intervals after start moment
            self.partial_instance.J_inverse['j'+str(j)] = [i for i in self.J_inverse['j'+str(j)] if i in self.partial_instance.I_a]
        return

    def find_partial_J(self):
        self.partial_instance.J = {}
        #take only those intervals that are in self.instance.partial_instance.I_a
        for i in self.partial_instance.I_a:
            #take self.instance.J and remove jobs that are not in partial instance
            self.partial_instance.J['i'+str(i)] = [j for j in self.J['i'+str(i)] if j in self.partial_instance.jobs]
        return
    
    def add_capacity(self,cap, constant = True, activate = True):
        #define list with global maximum powers per interval, in kW/timestep
        
        print('[WARNING]: global cap feature broken! Dont use. We set self.global_cap_active to False now.')
        activate = False

        if constant == True:
            self.global_cap = [cap for i in self.I_a]
        else:
            self.global_cap = cap
        if len(self.global_cap) != len(self.I_a):
            print('[WARNING]: global capacity defined for a number of intervals different from the number of active atomic intervals in the problem instance. Check add_capacity() documentation.')

        #turn on constraint
        self.global_cap_active = activate
        return

    def validate_solution(self,flow_dict, err = 0.0000001):
        print("################################################")
        print("Begin validation of solution")
        invalid = False
        demandNotMet = []
        for j in self.jobs:
            charge = 0
            for i in self.J_inverse["j" + str(j)]:
                charge += flow_dict["j" + str(j)]["i" + str(i)]

                #check feasibility - charge at no more than max power
                if flow_dict["j" + str(j)]["i" + str(i)] > self.jobs_cap[j]*self.len_i[i]:
                    print("[WARNING] Schedule exceeds maximum power of job {} during interval {}. Charge = {}".format(j,i,flow_dict["j" + str(j)]["i" + str(i)]))
                    invalid = True

            demandNotMet = demandNotMet + [self.jobs_demand[j] - charge]
            if demandNotMet[-1] > err:
                print("[WARNING] Energy not served to EV {} exceeds defined error. ENS_j = {}".format(j,demandNotMet[-1]))
                invalid = True

            #check feasibility - charge only if available
            for i in [i for i in self.I_a if i not in self.J_inverse["j"+ str(j)]]:
                try: 
                    flow_dict["j"+str(j)]["i"+str(i)]
                    print("[WARNING] Job {} is charging in interval {} without being available. Charge = {}".format(j,i,flow_dict["j"+str(j)]["i"+str(i)]))
                    invalid = True
                except:
                    pass

        #check feasibility - receive approximately all charge
        if sum(demandNotMet) > err:
            print("[WARNING] Total energy not served exceeds defined error. ENS = ", sum(demandNotMet), "(in validation check)", sum(self.jobs_demand) - sum(flow_dict["s"].values()), "(total demand - outflow of s)")
            invalid = True
        if not invalid:
            print("Validation check did not detect any infeasibilities.")
        print("Conclude validation of solution")
        print("################################################")
        return invalid
    
    def toy_instance_1(self):
        # overwrite instance with toy example to validate
        self.timeStep = 3600
        self.tau = self.timeStep/3600

        #initialize jobs
        self.n = 2 #number of jobs
        self.jobs = [j for j in range(0,self.n)] 
        self.jobs_cap = [2,2]
        self.jobs_demand = [2,2]

        #initialize intervals
        self.breakpoints = [0,1,2,3]
        self.intervals_start = [0,1,2]
        self.intervals_end = [1,2,3]
        self.m = len(self.intervals_start) #number of atomic intervals
        self.I_a = [i for i in range(0,self.m)]
        self.len_i = [(self.intervals_end[i] - self.intervals_start[i])*self.timeStep for i in self.I_a]

        self.J_inverse = {"j0": [0,1,2], "j1": [1]}
        self.J = {"i0": [0], "i1": [0,1], "i2": [0]}
        
        return
    
    def toy_instance_2(self):
        #toy instance with non-unit size intervals and realistic timeStep/timeBase. Result should be entirely flat.

        # overwrite instance with toy example to validate
        self.timeStep = 900
        self.timeBase = 3600
        self.tau = self.timeStep/self.timeBase

        #initialize jobs
        self.n = 2 #number of jobs
        self.jobs = [j for j in range(0,self.n)] 
        self.jobs_cap = [2,2] #in kWh
        self.jobs_demand = [2,2] #in kWh

        #initialize intervals
        self.breakpoints = [0,4,12,16]
        self.intervals_start = [0,4,12]
        self.intervals_end = [4,12,16]
        self.m = len(self.intervals_start) #number of atomic intervals
        self.I_a = [i for i in range(0,self.m)]
        self.len_i = [(self.intervals_end[i] - self.intervals_start[i])*self.timeStep for i in self.I_a]
        self.J_inverse = {"j0": [0,1,2], "j1": [1]}
        self.J = {"i0": [0], "i1": [0,1], "i2": [0]}
        return

    def toy_instance_3(self):
        #toy instance 2 with periodicity. Result should be entirely flat.

        # overwrite instance with toy example to validate
        self.timeStep = 900
        self.timeBase = 3600
        self.tau = self.timeStep/self.timeBase

        #initialize jobs
        self.n = 2 #number of jobs
        self.jobs = [j for j in range(0,self.n)] 
        self.jobs_cap = [2,2] #in kW
        self.jobs_demand = [2,2] #in kWh

        #initialize intervals
        self.breakpoints = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
        self.intervals_start = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
        self.intervals_end = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
        self.m = len(self.intervals_start) #number of atomic intervals
        self.I_a = [i for i in range(0,self.m)]
        self.len_i = [(self.intervals_end[i] - self.intervals_start[i])*self.timeStep for i in self.I_a]
        self.J_inverse = {"j0": [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15], "j1": [4,5,6,7,8,9,10,11]}
        self.J = {"i0": [0], "i1": [0], "i2": [0], "i3": [0], "i4": [0,1], "i5": [0,1], "i6": [0,1], "i7": [0,1], "i8": [0,1], "i9": [0,1], "i10": [0,1], "i11": [0,1], "i12": [0], "i13": [0], "i14": [0], "i15": [0]}

    def toy_instance_4(self):
        # overwrite instance with toy example to validate
        self.timeStep = 3600
        self.tau = self.timeStep/3600

        #initialize jobs
        self.n = 2 #number of jobs
        self.jobs = [j for j in range(0,self.n)] 
        self.jobs_cap = [4,4]
        self.jobs_demand = [6,6]

        #segments
        self.seg_n = [2,2] #number of segments
        self.seg_jobs = [[0,1],[0,1]]
        self.seg_jobs_demand = [[4.5,1.5],[3,3]]

        #initialize intervals
        self.breakpoints = [0,1,2,3]
        self.intervals_start = [0,1,2]
        self.intervals_end = [1,2,3]
        self.m = len(self.intervals_start) #number of atomic intervals
        self.I_a = [i for i in range(0,self.m)]
        self.len_i = [(self.intervals_end[i] - self.intervals_start[i])*self.timeStep for i in self.I_a]

        self.J_inverse = {"j0": [0,1,2], "j1": [0,1,2]}
        self.J = {"i0": [0,1], "i1": [0,1], "i2": [0,1]}

        self.seg_J_inverse = {"j0-0": [0,1], "j0-1": [0,1,2], "j1-0": [0], "j1-1": [0,1,2]}
        self.seg_J = {"i0": {0: [0,1], 1: [0,1]}, "i1": {0: [0,1], 1: [1]}, "i2": {0:[1], 1:[1]}}
        return
    
    def toy_instance_5(self):
        # overwrite instance with toy example to validate
        self.timeStep = 3600
        self.tau = self.timeStep/3600

        #initialize jobs
        self.n = 1 #number of jobs
        self.jobs = [j for j in range(0,self.n)] 
        self.jobs_cap = [4]
        self.jobs_demand = [4]

        #segments
        self.seg_n = [2] #number of segments
        self.seg_jobs = [[0,1]]
        self.seg_jobs_demand = [[3,1]]

        #initialize intervals
        self.breakpoints = [0,1,2]
        self.intervals_start = [0,1]
        self.intervals_end = [1,2]
        self.m = len(self.intervals_start) #number of atomic intervals
        self.I_a = [i for i in range(0,self.m)]
        self.len_i = [(self.intervals_end[i] - self.intervals_start[i])*self.timeStep for i in self.I_a]

        self.J_inverse = {"j0": [0,1]}
        self.J = {"i0": [0], "i1": [0]}

        self.seg_J_inverse = {"j0-0": [0], "j0-1": [0,1]}
        self.seg_J = {"i0": {0: [0,1]}, "i1": {0: [1]}}
        return

class FlowNet:
    def __init__(self):
        self.G = nx.DiGraph()

    def focs_instance_to_network(self, instance):
        #add nodes
        self.G.add_nodes_from(["s"] + ["j"+str(j) for j in instance.jobs] + ["i"+str(i) for i in instance.I_a] + ["t"])

        #add edges
        #later, initialize without edges to t and add those based on known capacity
        self.G.add_edges_from(
                            [("s", "j"+str(j), {"capacity": instance.jobs_demand[instance.jobs.index(j)]}) for j in instance.jobs]   #D_s
                        +  sum([[("j"+str(j), "i"+str(i), {"capacity": instance.jobs_cap[instance.jobs.index(j)]*instance.len_i[i]/instance.timeBase}) for i in instance.J_inverse["j"+str(j)]] for j in instance.jobs],[])   #D_0
                        +  [("i"+str(i), "t") for i in instance.I_a] #D_t
                        )
        
    def focs_solution_to_network(self, instance, f, err = 0.000001, how = 'linear', factor = 1000, structure = "focs", custom = None):
        #make empty network
        if structure == "focs":
            self.focs_instance_to_network(instance)
            self.G_sched = copy.deepcopy(self.G)
        elif structure == "segment":
            self.segment_instance_to_network(instance)
            self.G_sched = copy.deepcopy(self.G_seg)
        else: 
            print("[WARNING]: structure undefined. Try focs or segment.")
        # self.G_sched = copy.deepcopy(self.G)
        #add capacities corresponding to focs.f solution
        for i in instance.I_a:
            self.G_sched['i'+str(i)]['t']['capacity'] = f['i'+str(i)]['t']

        #per job, determine costs 
        if not instance.periodicity:
            print('[WARNING]: the method does not normalize for interval lengths. Periodicity assumed and detected as False.')   
        if structure == "focs":
            for j in instance.jobs:
                n_min = int(instance.jobs_demand[j]/(instance.tau*instance.jobs_cap[j])) #rounded minimal number of intervals needed for charging
                n_slack = len(instance.J_inverse['j'+str(j)][n_min+1:])
                if custom is not None:
                    how = custom[instance.jobs.index(j)]
                for i_id, i in enumerate(instance.J_inverse['j'+str(j)][n_min+1:]):
                    if how == 'linear':
                        self.G_sched['j'+str(j)]['i'+str(i)]['weight'] = i_id + 1
                    elif how == 'quadratic':
                        self.G_sched['j'+str(j)]['i'+str(i)]['weight'] = int((i_id + 1)**2)
                    elif how == 'linear_proportional':
                        self.G_sched['j'+str(j)]['i'+str(i)]['weight'] = int(factor*(i_id +1) /(n_slack))
                    elif how == 'quadratic_proportional':
                        self.G_sched['j'+str(j)]['i'+str(i)]['weight'] = int(factor*(i_id +1)**2 /(n_slack))
                    elif how == 'probabilistic':
                        self.G_sched['j'+str(j)]['i'+str(i)]['weight'] = int(factor/(instance.jobs_cumprob[j][i]+err))
                    elif how == 'probabilistic_full':
                        for i in instance.J_inverse['j'+str(j)]:
                            self.G_sched['j'+str(j)]['i'+str(i)]['weight'] = int(factor/(instance.jobs_cumprob[j][i]+err))
                    elif how == 'reverse':
                        self.G_sched['j'+str(j)]['i'+str(instance.J_inverse['j'+str(j)][-i_id - n_min - 1 -1])]['weight'] = i_id + 1
                    else:
                        print('[WARNING]: how undefined. Try \'linear\' or \'quadratic\'. No weight assigned.')
                        print('how = ', how)
        #add demands
        demand = max(0, sum([self.G_sched['s']['j'+str(j)]['capacity'] for j in instance.jobs]) - err)
        nx.set_node_attributes(self.G_sched, {'s': - demand, 't': demand}, name = 'demand')

        return self.G_sched

    def segment_instance_to_network(self,instance):
        self.G_seg = nx.DiGraph()
        # add nodes
        self.G_seg.add_nodes_from(["s"] + ["j"+str(j) for j in instance.jobs] + ["p"+str(p)+"i"+str(i) for i in instance.I_a for p in instance.parent_instance.J['i'+str(i)]] + ["i"+str(i) for i in instance.I_a] + ["t"])

        # add edges
        self.G_seg.add_edges_from(
                        [("s", "j"+str(j), {"capacity": instance.jobs_demand[instance.jobs.index(j)]}) for j in instance.jobs]   #D_0
                        +  sum([[("j"+str(j), "p"+str(instance.parent_jobs[j])+"i"+str(i), {"capacity": instance.jobs_cap[instance.jobs.index(j)]*instance.len_i[i]/instance.timeBase}) for i in instance.J_inverse["j"+str(j)]] for j in instance.jobs],[])   #D_1
                        +  sum([[("p"+str(p)+"i"+str(i), "i"+str(i), {"capacity": instance.parent_instance.jobs_cap[instance.parent_instance.jobs.index(p)]*instance.len_i[i]/instance.timeBase}) for i in instance.parent_instance.J_inverse["j"+str(p)]] for p in instance.parent_instance.jobs],[])   #D_2
                        +  [("i"+str(i), "t") for i in instance.I_a] #D_t
                        )
                        
class FlowOperations:
    def __init__(self, G, instance):
        self.G = G
        self.instance = instance

    #initialize empty flow
    def empty_flow(self, G = None):
        if G is None:
            G = self.G
        f = {}
        for u in G.nodes:
            f[u] = {}
        for u,v in G.edges:
            f[u][v] = 0
        return f
    
    def solve_flow(self):
        self.flow_val, self.flow_dict = nx.maximum_flow(self.G, "s", "t", flow_func=shortest_augmenting_path)

    #adding one flow to another
    def add_flows(self,flow1, flow2, G2):
        for (u,v) in G2.edges:
            flow1[u][v] += flow2[u][v]
        return flow1
    
    #adding price to focs flows
    def add_price_to_focs_flow(self,flow1,flow2,G1,G2):
        for (u,v) in G2.edges:
            try:
                flow1[u][v] += flow2[u][v]
            except:
                # print('sanity check only edges through c-nodes are skipped: ', u,'\t',v)
                pass
        # make feasible - incoming flow in t:
        for u in G1.nodes:
            if (u,'t') in G1.edges:
                temp = 0
                for v in G1.nodes:
                    if (v,u) in G1.edges:
                        temp += flow1[v][u]
                flow1[u]['t'] = temp

    #Reduce flow network by critical flow
    def reduce_network(self,Gr,fr_dict,I_critr):
        #FIXME make graph object, self.G should be changed.
        for i in I_critr:
            Gr.remove_node("i{}".format(i))
        for (u,v) in Gr.edges:
            #FIXME how to loop over flow elements, not use try except construction?
            try:
                Gr[u][v]["capacity"] -= fr_dict[u][v]
            except:
                pass
        return Gr #=G_r+1

    #L(I) function in paper. Returns 0 if I_list is empty.
    def length_sum_intervals(self,I_list, L_list):
        return sum([L_list[i] for i in I_list])

    def partial_flow(self,full_flow, node_set, intervals, jobs, J, structure = "focs", parent_J = {}, parent_jobs = []):
        partial_flow = copy.deepcopy(full_flow)
        if structure == "focs":
            for i in intervals:
                if i not in node_set:
                    partial_flow["i{}".format(i)]["t"] = 0
                    for j in J["i" + str(i)]:
                        #FIXME double check definition J[i]. Need to update function if we remove nodes... Active jobs only.
                        partial_flow["j{}".format(j)]["i{}".format(i)] = 0
            for j in jobs:
                partial_flow["s"]["j{}".format(j)] = sum(partial_flow["j{}".format(j)].values())
        elif structure == "segment":
            for i in intervals:
                if i not in node_set: #node set = critical intervals
                    partial_flow["i"+str(i)]["t"] = 0
                    for p in parent_J["i" + str(i)]:
                        partial_flow["p"+str(p)+"i"+str(i)]["i" + str(i)] = 0
                    for j in J["i"+str(i)]:
                        partial_flow["j"+str(j)]["p"+str(parent_jobs[j])+"i"+str(i)] = 0
            for j in jobs:
                partial_flow["s"]["j"+str(j)] = sum(partial_flow["j"+str(j)].values())
        else:
            print("[WARNING]: structure for partial flow unknown. Try focs or segment.")
        return partial_flow
    
    def validate_global_cap(self, G=None, global_cap=None, I_a = None, demand = None, err = 0.00000001):
        #given a power profile, check that there exists a feasible solution
        print("################################################")
        print("Begin validation of feasibility of global power profile")
        if G is None:
            G = self.G
        if global_cap is None:
            global_cap = self.instance.global_cap
        if I_a is None:
            I_a = self.instance.I_a
        if demand == None:
            demand = sum([G["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs])

        #define power capacities
        for i in I_a:
            G["i{}".format(i)]["t"]["capacity"] = self.instance.global_cap[i]*self.instance.len_i[i]*self.instance.tau/self.instance.timeStep   

        #run max flow
        flow_val, flow_dict = nx.maximum_flow(G, "s", "t", flow_func=shortest_augmenting_path)

        #feasibility check
        if demand - flow_val > err:
            print('[WARNING]: The demand cannot be met given the defined global power profile and an error margin of ', err)
        else: 
            print('[MESSAGE]: There exists a feasible solution given the global power profile')
        print("################################################")
        return


class FOCS:
    def __init__(self, instance, flowNet, flowOp):
        self.instance = instance
        self.flowNet = flowNet
        self.flowOp = flowOp
        self.I_p = []
        self.I_a = instance.I_a
        self.global_cap_active = instance.global_cap_active
        self.no_demand = False # sanity check

        self.I_crit = [] #end of each round, save the list of active intervals here
        self.J = instance.J #per interval, specify a list of jobs available
        self.J_inverse = instance.J_inverse #per job, specify which intervals it is available
        self.rd = 0 #round counter
        self.it = 0 #iteration counter
        self.terminate = False
        self.G = flowNet.G
        self.G_r = self.G
        self.total_demand_r = sum([self.G_r["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs])
        self.flow_func = shortest_augmenting_path #(or e.g., edmonds_karp)
        self.flow_val = 0
        self.f = flowOp.empty_flow(self.G)

    #update network capacities (g function in paper)
    def update_network_capacities_g(self,G_r,flow_val): #doesn't have to be same G_r as in self. May be G_rk
        if self.it > 0:
            #print("update capacities of round ", self.rd, "iteration ", self.it)
            demand = self.maxDiff
            demand_normalized = demand/self.flowOp.length_sum_intervals(self.I_a,self.instance.len_i)
            for i in self.I_a:
                G_r["i{}".format(i)]["t"]["capacity"] += demand_normalized * self.instance.len_i[i]
        else:
            #print("initialize capacities of round ", self.rd)
            demand = sum([G_r["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs])
            if demand == 0:
                print('[WARNING]: instance has no energy demand. Return empty flow as solution.')
                self.no_demand = True
                return G_r
            demand_normalized = demand/self.flowOp.length_sum_intervals(self.I_a, self.instance.len_i)         
            for i in self.I_a:
                G_r["i{}".format(i)]["t"]["capacity"] = demand_normalized * self.instance.len_i[i]

        #if global maximum power constraint is active, trigger it here
        if self.instance.global_cap_active:
            for i in self.I_a:
                G_r["i{}".format(i)]["t"]["capacity"] = min(G_r["i{}".format(i)]["t"]["capacity"], self.instance.global_cap[i]*self.instance.len_i[i]*self.instance.tau/self.instance.timeStep)        
        return G_r
    
    def solve_focs(self, err = 0.0000001, MPCstopper = False, MPCcondition = 0):
        while not self.terminate:
            #begin round
            #initiate capacities
            if self.it == 0:
                G_rk = copy.deepcopy(self.G_r)
            G_rk = self.update_network_capacities_g(G_rk, self.flow_val)
            if self.no_demand:
                return self.f

            #determine max flow
            self.flow_val, flow_dict = nx.maximum_flow(G_rk, "s", "t", flow_func=self.flow_func)
            self.maxDiff = sum([G_rk["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs]) - self.flow_val
            if self.total_demand_r - self.flow_val < err:
                #end round

                #Check for subcrit. For some instances, there are still active intervals that only now don't reach the max anymore #FIXME (for now)
                subCrit_mask = [G_rk["i{}".format(i)]["t"]["capacity"]-flow_dict["i{}".format(i)]["t"] > err for i in self.I_a]
                subCrit = [self.I_a[i] for i in range(0,len(self.I_a)) if subCrit_mask[i]]
                
                #Update I_p
                self.I_p += subCrit
                #Update I_crit
                self.I_crit = self.I_crit + [sorted([self.I_a[i] for i in range(0,len(self.I_a)) if not subCrit_mask[i]])]

                #check terminate
                if len(self.I_p) == 0:
                    self.terminate = True
                    #Add critical flow to current flow
                    self.f = self.flowOp.add_flows(self.f, flow_dict, G_rk)
                else:
                    #determine critical flow
                    I_crit_r = self.I_crit[-1] 
                    flow_dict_crit = self.flowOp.partial_flow(flow_dict, I_crit_r, self.I_a, self.instance.jobs, self.instance.J) #FIXME check if J needs to be reduced?
                    
                    #Add critical flow to current flow
                    self.f = self.flowOp.add_flows(self.f,flow_dict_crit,G_rk)

                    #Reduce flow network by critical flow
                    self.G_r = self.flowOp.reduce_network(self.G_r,flow_dict_crit,I_crit_r)

                    #Update active and parked sets
                    self.I_a = sorted(self.I_p)
                    self.I_p = []

                    if MPCstopper:
                        self.MPCcondition = MPCcondition
                        if MPCcondition < min(self.I_a):
                            # print('[MESSAGE]: Solution optimal for first {} timesteps only.'.format(MPCcondition + 1))
                            return self.f

                    #amount of work yet to be scheduled
                    self.total_demand_r = sum([self.G_r["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs])

                    self.rd += 1
                    self.it = 0
            else:
                #initiate next iteration
                #Determine subcritical sets
                if self.instance.global_cap_active:
                    subCrit_mask = [(G_rk["i{}".format(i)]["t"]["capacity"] - flow_dict["i{}".format(i)]["t"] > err) or (self.instance.global_cap[i]*self.instance.len_i[i]*self.instance.tau/self.instance.timeStep == G_rk["i{}".format(i)]["t"]["capacity"]) for i in self.I_a]
                else:
                    subCrit_mask = [G_rk["i{}".format(i)]["t"]["capacity"] - flow_dict["i{}".format(i)]["t"] > err for i in self.I_a]
                subCrit = [self.I_a[i] for i in range(0,len(self.I_a)) if subCrit_mask[i]]

                #Reduce network G_rk
                flow_dict_sub = self.flowOp.partial_flow(flow_dict, subCrit, self.I_a, self.instance.jobs, self.instance.J)
                G_rk = self.flowOp.reduce_network(G_rk,flow_dict_sub,subCrit)
                self.total_demand_r = sum([G_rk["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs])
                
                #Update I_p
                self.I_p += subCrit

                #Update I_a
                self.I_a = sorted([self.I_a[i] for i in range(0,len(self.I_a)) if not subCrit_mask[i]])
                self.it += 1

        return self.f
    
    def solve_segmented_focs(self, err = 0.0000001, MPCstopper = False, MPCcondition = 0):
        while not self.terminate:
            #begin round
            #initiate capacities
            if self.it == 0:
                G_rk = copy.deepcopy(self.G_r)
            G_rk = self.update_network_capacities_g(G_rk, self.flow_val)

            #determine max flow
            self.flow_val, flow_dict = nx.maximum_flow(G_rk, "s", "t", flow_func=self.flow_func)
            self.maxDiff = sum([G_rk["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs]) - self.flow_val
            if self.total_demand_r - self.flow_val < err:
                #end round

                #Check for subcrit. For some instances, there are still active intervals that only now don't reach the max anymore #FIXME (for now)
                subCrit_mask = [G_rk["i{}".format(i)]["t"]["capacity"]-flow_dict["i{}".format(i)]["t"] > err for i in self.I_a]
                subCrit = [self.I_a[i] for i in range(0,len(self.I_a)) if subCrit_mask[i]]
                
                #Update I_p
                self.I_p += subCrit
                #Update I_crit
                self.I_crit = self.I_crit + [sorted([self.I_a[i] for i in range(0,len(self.I_a)) if not subCrit_mask[i]])]

                #check terminate
                if len(self.I_p) == 0:
                    self.terminate = True
                    #Add critical flow to current flow
                    self.f = self.flowOp.add_flows(self.f, flow_dict, G_rk)
                else:
                    #determine critical flow
                    I_crit_r = self.I_crit[-1] 
                    flow_dict_crit = self.flowOp.partial_flow(flow_dict, I_crit_r, self.I_a, self.instance.jobs, self.instance.J, structure="segment", parent_J=self.instance.parent_instance.J, parent_jobs=self.instance.parent_jobs) 
                    
                    #Add critical flow to current flow
                    self.f = self.flowOp.add_flows(self.f,flow_dict_crit,G_rk)

                    #Reduce flow network by critical flow
                    self.G_r = self.flowOp.reduce_network(self.G_r,flow_dict_crit,I_crit_r)

                    #Update active and parked sets
                    self.I_a = sorted(self.I_p)
                    self.I_p = []

                    if MPCstopper:
                        self.MPCcondition = MPCcondition
                        if MPCcondition < min(self.I_a):
                            # print('[MESSAGE]: Solution optimal for first {} timesteps only.'.format(MPCcondition + 1))
                            return self.f

                    #amount of work yet to be scheduled
                    self.total_demand_r = sum([self.G_r["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs])

                    self.rd += 1
                    self.it = 0
            else:
                #initiate next iteration
                #Determine subcritical sets
                subCrit_mask = [G_rk["i{}".format(i)]["t"]["capacity"] - flow_dict["i{}".format(i)]["t"] > err for i in self.I_a]
                subCrit = [self.I_a[i] for i in range(0,len(self.I_a)) if subCrit_mask[i]]

                #Reduce network G_rk
                flow_dict_sub = self.flowOp.partial_flow(flow_dict, subCrit, self.I_a, self.instance.jobs, self.instance.J, structure = "segment", parent_J=self.instance.parent_instance.J, parent_jobs = self.instance.parent_jobs)
                G_rk = self.flowOp.reduce_network(G_rk,flow_dict_sub,subCrit)
                self.total_demand_r = sum([G_rk["s"]["j{}".format(j)]["capacity"] for j in self.instance.jobs])
                
                #Update I_p
                self.I_p += subCrit

                #Update I_a
                self.I_a = sorted([self.I_a[i] for i in range(0,len(self.I_a)) if not subCrit_mask[i]])
                self.it += 1

        return self.f

    def objective(self):
        #normalized objective for validation
        #instead of the sum of squares, we first weight the square of the power (e/len_i[i]) by len_i[i]/timeStep      
        #determine power per interval
        p_i = [((self.f["i"+str(i)]['t']/(self.instance.len_i[i]))*self.instance.timeBase) for i in range(0,self.instance.m)]
        #determine weighted squared powers
        powerSquare = [(p_i[i]**2)*(self.instance.len_i[i]/self.instance.timeStep) for i in range(0,self.instance.m)]
        self.objNormalized = sum(powerSquare)
        return self.objNormalized

    def partial_objective(self, start_h):
        #normalized objective for validation of remaining problem setup
        #choose breakpoint closest to start
        start = start_h/self.instance.tau
        iStartBp = min(self.instance.breakpoints, key= lambda x:abs(x-start))
        iStart = self.instance.breakpoints.index(iStartBp)
        
        #determine power per interval
        
        p_i = [(self.f["i"+str(i)]['t']/(self.instance.len_i[i]))*self.instance.timeBase for i in range(0,self.instance.m)]
        if iStartBp > start: #too short already.
            p_i = p_i[iStart-1:]
            factor = [(iStartBp - start)] + [(self.instance.len_i[iStart+i]/self.instance.timeStep) for i in range(0,len(p_i)-1)] 
        elif iStartBp < start: #case checked against full version and MPC partial
            #determine weighted squared powers
            p_i = p_i[iStart:]
            factor = [(self.instance.breakpoints[iStart +1] - start)] + [self.instance.len_i[iStart+1+i]/self.instance.timeStep for i in range(0,len(p_i)-1)]
        else: #case checked: MPC against itself
            p_i = p_i[iStart:]
            factor = [self.instance.len_i[iStart+i]/self.instance.timeStep for i in range(0,len(p_i))]
        powerSquare = [(p_i[i]**2)*factor[i] for i in range(0,len(p_i))]
        self.objNormalized = sum(powerSquare)
        return self.objNormalized
    
    def write(self, type = 'txt', suffix = ''):
        if type == 'json':
            print('[WARNING] Saving to json not supported yet')
            return
            with open("my.json","w") as file:
                json.dump(self.f,file)
            file.close()

        if type == 'pickle':
            file = open('C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data/output/f_{}.pickle'.format(suffix), 'wb')
            pickle.dump(self.f, file)
            file.close()

            #to retrieve:
            #file = open('store.pckl', 'rb')
            #obj = pickle.load(file)
            #file.close()
            #print(obj)

        if type == 'txt':
            file = open("C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data/output/f_{}.txt".format(suffix), "w+")
            file.write(json.dumps(self.f, indent=4, sort_keys=True))
            file.close()
        return
    
    def validate_MPC(self):
        print("[MESSAGE]: Begin validation of solution")
        invalid = False
        Loc = []

        #MPC solution
        flow = copy.deepcopy(self.f)

        #full solution
        self.flowNet.focs_instance_to_network(self.instance)
        focs = FOCS(self.instance, self.flowNet, self.flowOp)
        f = focs.solve_focs(MPCstopper=False)

        #compare solutions
        for i in range(0,self.MPCcondition + 1):
            if flow["i"+str(i)]['t'] != f["i"+str(i)]['t']:
                invalid = True
                Loc = Loc + [i]
        #jaay it's consistent
        if invalid:
            print('[WARNING]: Inconsistency detected for intervals i + {}'.format(Loc))
        else:
            print("Validation check did not detect any inconsistencies.")
        print("[MESSAGE]: Conclude validation of solution")
        return
    
class Schedule: #class introduced to a) extract specific schedules from optimal solutions b) choose a schedule with LYNCS c) evaluate schedules
    def __init__(self, focs):
        self.instance = focs.instance
        self.flowNet = focs.flowNet
        self.flowOp = focs.flowOp
        self.focs = focs
    
    def solve_schedule(self, how = 'linear', custom = None):
        #make min cost flow network
        self.G_sched = self.flowNet.focs_solution_to_network(self.instance, self.focs.f, how=how, custom=custom)
        #solve min cost flow
        self.f_sched = nx.min_cost_flow(self.G_sched)
        return self.f_sched
    
    def flow_to_dataframe(self, f = None):
        if f is None:
            f = self.f_sched
        # flow to dataframe - format schedule
        s = [self.focs.instance.I_a] + [self.focs.instance.intervals_start] + [[f['i'+str(i)]['t'] for i in self.focs.instance.I_a]] + [[0]*len(self.focs.instance.I_a) for j in self.focs.instance.jobs]
        for jid, j in enumerate(self.focs.instance.jobs):
            for i in self.focs.instance.J_inverse['j'+str(j)]:
                s[3+jid][i] = f['j'+str(j)]['i'+str(i)] / self.focs.instance.tau
        s = pd.DataFrame(s).T
        s.columns = ['time'] + ['breakpoints'] + ['agg'] + ['j'+str(j) for j in self.focs.instance.jobs]
        return s

    def derive_early_departure_schedule(self, key='t1_900', f=None):
        if f is None:
            f = self.f_sched
        # determine G and empty flow
        self.focs.flowNet.focs_instance_to_network(self.focs.instance)
        flowOp = FlowOperations(self.focs.flowNet.G, self.focs.instance)
        self.jobs_dep_real = [self.focs.instance.data[key].iloc[j] for j in self.focs.instance.jobs] 
        self.f_ed = flowOp.empty_flow()
        # fill empty flow
        for j in self.focs.instance.jobs:
            if self.jobs_dep_real[j] in self.focs.instance.breakpoints:
                temp_idx = self.focs.instance.breakpoints.index(self.jobs_dep_real[j])
            else:
                temp_idx = max(self.focs.instance.breakpoints)
            for i in range(min(self.focs.instance.J_inverse['j'+str(j)]), min(self.focs.instance.breakpoints.index(self.focs.instance.jobs_departure[j]), temp_idx)):
                self.f_ed['j'+str(j)]['i'+str(i)] = f['j'+str(j)]['i'+str(i)]
            self.f_ed['s']['j'+str(j)] = sum(self.f_ed['j'+str(j)].values())
        for i in self.focs.instance.I_a:
            self.f_ed['i'+str(i)]['t'] = sum([self.f_ed['j'+str(j)]['i'+str(i)] for j in self.focs.instance.J['i'+str(i)]])
        return self.f_ed
    
    def objective(self, f_plan = None, f_ed = None):
        if f_plan is None:
            f_plan = self.f_sched
        if f_ed is None:
            f_ed = self.derive_early_departure_schedule(f = f_plan)

        ## Metrics per job

        # ENS exact
        self.jobs_ens_abs_exact = [self.focs.f['s']['j'+str(j)] - f_ed['s']['j'+str(j)] for j in self.focs.instance.jobs]
        self.jobs_ens_rel_exact = [self.jobs_ens_abs_exact[j]/self.focs.flowNet.G['s']['j'+str(j)]['capacity'] for j in self.focs.instance.jobs]
        # ENS rounded
        self.jobs_ens_abs = [round(self.focs.f['s']['j'+str(j)] - f_ed['s']['j'+str(j)],6) for j in self.focs.instance.jobs]
        self.jobs_ens_rel = [round(self.jobs_ens_abs_exact[j]/self.focs.flowNet.G['s']['j'+str(j)]['capacity'],6) for j in self.focs.instance.jobs]
        
        # QOS1 : relative energy served
        self.jobs_es_rel = [1 - self.jobs_ens_rel[j] for j in self.focs.instance.jobs]

        # QOS2 : relative wait till charging starts
        self.jobs_qos2_waiting_real, self.jobs_qos2_waiting_plan = self.qos_2(f_plan=f_plan)

        # QOS3 : FIXME variation of charging power over time
        self.jobs_qos3_powervar_real, self.jobs_qos3_powervar_plan = self.qos_3(f_plan = f_plan, f_ed = f_ed)

        # QOE : 
        self.jobs_qoe_real_ed, self.jobs_qoe_real_pot, self.jobs_qoe_plan_ed, self.jobs_qoe_plan_pot = self.qoe(f_plan=f_plan,f_ed=f_ed)

        ## Global metrics

        # #QOE total
        self.qoe_real_ed_total_exact = sum(self.jobs_qoe_real_ed)
        self.qoe_real_pot_total_exact = sum(self.jobs_qoe_real_pot)
        self.qoe_plan_ed_total_exact = sum(self.jobs_qoe_plan_ed)        
        self.qoe_plan_pot_total_exact = sum(self.jobs_qoe_plan_pot)        
        self.qoe_real_ed_total_rel = self.qoe_real_ed_total_exact/self.focs.instance.n
        self.qoe_real_pot_total_rel = self.qoe_real_pot_total_exact/self.focs.instance.n
        self.qoe_plan_ed_total_rel = self.qoe_plan_ed_total_exact/self.focs.instance.n
        self.qoe_plan_pot_total_rel = self.qoe_plan_pot_total_exact/self.focs.instance.n

        # ENS max
        self.jobs_ens_abs_max = max(self.jobs_ens_abs)
        self.jobs_ens_rel_max = max(self.jobs_ens_rel)

        # QOS1 min
        self.qos_1_min = min(self.jobs_es_rel)
        # QOS2 min
        self.qos_2_plan_min = min(self.jobs_qos2_waiting_plan)
        self.qos_2_real_min = min(self.jobs_qos2_waiting_real)
        # QOS3 min
        self.qos_3_plan_min = min(self.jobs_qos3_powervar_plan)
        self.qos_3_real_min = min(self.jobs_qos3_powervar_real)

        # Jain's fairness index based on relative ENS
        self.jain_ens_rel = self.Jain(self.jobs_ens_rel)
        self.jain_ens_rel_exact = self.Jain(self.jobs_ens_rel_exact)

        # Jain's fairness index based on qos and qoe
        self.jain_qos_1 = self.Jain(self.jobs_es_rel)
        self.jain_qos_2_plan = self.Jain(self.jobs_qos2_waiting_plan)
        self.jain_qos_2_real = self.Jain(self.jobs_qos2_waiting_real)
        self.jain_qos_3_plan = self.Jain(self.jobs_qos3_powervar_plan)
        self.jain_qos_3_real = self.Jain(self.jobs_qos3_powervar_real)

        # Hoßfeld's fairness index based on relative ENS
        self.hossfeld_ens_rel = self.Hossfeld(self.jobs_ens_rel)
        self.hossfeld_ens_rel_exact = self.Hossfeld(self.jobs_ens_rel_exact)

        # Hoßfeld's fairness index based on qos and qoe
        self.hossfeld_qos_1 = self.Hossfeld(self.jobs_es_rel)
        self.hossfeld_qos_2_plan = self.Hossfeld(self.jobs_qos2_waiting_plan)
        self.hossfeld_qos_2_real = self.Hossfeld(self.jobs_qos2_waiting_real)
        self.hossfeld_qos_3_plan = self.Hossfeld(self.jobs_qos3_powervar_plan)
        self.hossfeld_qos_3_real = self.Hossfeld(self.jobs_qos3_powervar_real)

        # cycle switches TODO

        # energy served
        self.es_exact = sum(self.focs.instance.jobs_demand) - sum(self.jobs_ens_abs_exact)
        self.es = round(sum(self.focs.instance.jobs_demand) - sum(self.jobs_ens_abs),6)
        # energy not served
        self.ens_abs_exact = sum(self.jobs_ens_abs_exact)
        self.ens_abs = sum(self.jobs_ens_abs)
        self.ens_rel_exact_avg = sum(self.jobs_ens_rel_exact)/len(self.jobs_ens_rel_exact)
        self.ens_rel_avg = sum(self.jobs_ens_rel)/len(self.jobs_ens_rel)

        return
    
    def Jain(self, x):
        return sum(x)**2/(len(x) * sum([x_i**2 for x_i in x]))
    
    def Hossfeld(self, x, H = 1, h = 0):
        return 1 - 2* statistics.pstdev(x)/(H-h)
    
    def qos_2(self, f_plan = None): # qos_2 in Danner and de Meer (2021)
        if f_plan is None:
            f_plan = self.f_sched

        # determine first interval with positive charge
        first_pos_charging_power = [list(f_plan['j'+str(j)].values()).index([i for i in list(f_plan['j'+str(j)].values()) if i != 0][0]) for j in self.focs.instance.jobs]
        idx_start_charging = [self.focs.instance.J_inverse['j'+str(j)][first_pos_charging_power[self.focs.instance.jobs.index(j)]] for j in self.focs.instance.jobs]
        t_start_charging = [self.focs.instance.breakpoints[idx_start_charging[self.focs.instance.jobs.index(j)]] for j in self.focs.instance.jobs]

        # determine 
        self.jobs_qos2_waiting_plan = [max(0,1 - (t_start_charging[j] - self.focs.instance.jobs_arrival[j])/(self.focs.instance.jobs_departure[j] - self.focs.instance.jobs_arrival[j])) for j in self.focs.instance.jobs]

        # determine 
        self.jobs_qos2_waiting_real = [max(0,1 - (t_start_charging[j] - self.focs.instance.jobs_arrival[j])/(self.jobs_dep_real[j] - self.focs.instance.jobs_arrival[j])) for j in self.focs.instance.jobs]

        return self.jobs_qos2_waiting_real, self.jobs_qos2_waiting_plan
    
    def qos_3(self,f_plan = None, f_ed = None):
    # NOTE Checking correctness of the statements by Danner and de Meer with Maria rn.
    # NOTE Checked. They use the biased sample standard deviation (divide by N). See https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9639.1980.tb00398.x for proof of bound.
        if f_plan is None:
            f_plan = self.f_sched 
        if f_ed is None:
            f_ed = self.f_ed

        # create lists to store results. Index = job index.
        jobs_qos3_powervar_real = []
        jobs_qos3_powervar_plan = []

        # make schedule
        s = self.flow_to_dataframe(f_plan)
        s_ed = self.flow_to_dataframe(f_ed)

        for j in self.focs.instance.jobs:
            # isolate list of powers
            s_j = list(s['j'+str(j)])
            s_ed_j = list(s_ed['j'+str(j)])

            for case in [[s_j, jobs_qos3_powervar_plan, 'plan'], [s_ed_j, jobs_qos3_powervar_real, 'real']]:

                # retrieve first and last non-zero index
                temp = [i for i,e in enumerate(case[0]) if e!=0]

                # raise exception if no charging at all (statistics error std with empty sample)
                if len(temp) == 0:
                    print('[WARNING]: EV j = {} does not charge in {} case. QoS_3 set to 1.'.format(j, case[2]))
                    case[1] += [1]
                    continue

                t_start = min(temp)
                t_end = max(temp)+1

                # qos for j add to list
                case [1] += [1 - (2*statistics.pstdev(case[0][t_start:t_end])/self.focs.instance.jobs_cap[j]) ]

        return jobs_qos3_powervar_real, jobs_qos3_powervar_plan
    
    def supplied_by_milestone_ed(self,f,j,milestone_t):
            # determine 
            # FIXME this assumes periodicity
            if not self.focs.instance.periodicity:
                print('[WARNING]: Function supplied_by_milestone_ed() assumed periodicity == True. Here it is False.')
            # FIXME this assumes milestones in timestep units
            
            supplied = sum([f['j'+str(j)]['i'+str(i)] for i in self.focs.instance.J_inverse['j'+str(j)] if i < min(milestone_t, self.jobs_dep_real[j])])

            if milestone_t < min(self.focs.instance.J_inverse['j'+str(j)])+1:
                print('[WARNING]: milestone deadline before arrival of the EV.') 
            # early departure case: EV leaves before milestone becomes applicable. 
            elif milestone_t > min(self.jobs_dep_real[j], max(self.focs.instance.J_inverse['j'+str(j)])+1):
                supplied_pot = supplied + (milestone_t - min(self.jobs_dep_real[j], max(self.focs.instance.J_inverse['j'+str(j)])+1))*self.focs.instance.jobs_cap[j]*self.focs.instance.tau # if had charged at full power after time of departure, would have made milestone?
            else: 
                supplied_pot = supplied

            return supplied, supplied_pot #supplied is what has been supplied. supplied_pot is what could have been supplied if driver had stayed until milestone. 

    def qoe_j_ed(self,f,j,milestones_j, condition = 'supplied', err=0.0000001):
        for milestone in milestones_j:
                supplied, supplied_pot = self.supplied_by_milestone_ed(f,j,milestone[1])
                # check if made milestone or fully charged for this milestone. 
                if condition == 'supplied':
                    if supplied < min(milestone[0],self.focs.instance.jobs_demand[j])-err:
                        return 0
                elif condition == 'supplied_pot':
                    if supplied_pot < min(milestone[0],self.focs.instance.jobs_demand[j])-err:
                        return 0
                else:
                    print('[WARNING]: Condition value has to be supplied or supplied_pot. Other values not known.')
        return 1
    
    def supplied_by_milestone_plan(self,f,j,milestone_t):
            # determine 
            # FIXME this assumes periodicity
            if not self.focs.instance.periodicity:
                print('[WARNING]: Function supplied_by_milestone_plan() assumed periodicity == True. Here it is False.')
            # FIXME this assumes milestones in timestep units
            
            supplied = sum([f['j'+str(j)]['i'+str(i)] for i in self.focs.instance.J_inverse['j'+str(j)] if i < milestone_t])

            if milestone_t < min(self.focs.instance.J_inverse['j'+str(j)])+1:
                print('[WARNING]: milestone deadline before arrival of the EV.') 
            # early departure case: EV leaves before milestone becomes applicable. 
            elif milestone_t > max(self.focs.instance.J_inverse['j'+str(j)])+1:
                supplied_pot = supplied + (milestone_t - (max(self.focs.instance.J_inverse['j'+str(j)])+1))*self.focs.instance.jobs_cap[j]*self.focs.instance.tau # if had charged at full power after time of (planned/expected) departure, would have made milestone?
            else: 
                supplied_pot = supplied

            return supplied, supplied_pot #supplied is what has been supplied. supplied_pot is what could have been supplied if driver had stayed until milestone. 

    def qoe_j_plan(self,f,j,milestones_j, condition = 'supplied', err=0.0000001):
        for milestone in milestones_j:
                supplied, supplied_pot = self.supplied_by_milestone_plan(f,j,milestone[1])
                # check if made milestone or fully charged for this milestone. 
                if condition == 'supplied':
                    if supplied < min(milestone[0],self.focs.instance.jobs_demand[j])-err:
                        return 0
                elif condition == 'supplied_pot':
                    if supplied_pot < min(milestone[0],self.focs.instance.jobs_demand[j])-err:
                        return 0
                else:
                    print('[WARNING]: Condition value has to be supplied or supplied_pot. Other values not known.')
        return 1

    def qoe(self, f_plan = None, f_ed = None, err = 0.0000001): # check milestones and guarantees.
        # Binary metric per EV. 1 if requirements met. 0 if not.
        if f_plan is None:
            f_plan = self.f_sched
        if f_ed is None:
            f_ed = self.f_ed
        # if button is off, make guarantee static asr
        if not self.focs.instance.milestones_active:
            # add asr milestones: 23kWh by 4pm (if quarters then by quarter 64)
            milestones = [[[23,64]] for j in self.focs.instance.jobs] # if multiple guarantees, energy is additive!
        else:
            milestones = self.focs.instance.milestones

        jobs_qoe_real_ed = [self.qoe_j_ed(f_ed, j, milestones[j],condition='supplied') for j in self.focs.instance.jobs]
        jobs_qoe_real_pot = [self.qoe_j_ed(f_ed, j, milestones[j],condition='supplied_pot')for j in self.focs.instance.jobs]
        jobs_qoe_plan_ed = [self.qoe_j_plan(f_plan, j, milestones[j],condition='supplied') for j in self.focs.instance.jobs]
        jobs_qoe_plan_pot = [self.qoe_j_plan(f_plan, j, milestones[j],condition='supplied_pot')for j in self.focs.instance.jobs]

        return jobs_qoe_real_ed, jobs_qoe_real_pot, jobs_qoe_plan_ed, jobs_qoe_plan_pot