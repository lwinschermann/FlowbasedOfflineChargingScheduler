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

class FOCSinstance:
    def __init__(self, instanceData, timeStep, timeBase = 3600):
        self.data = instanceData
        self.timeStep = timeStep
        self.timeBase = timeBase
        self.tau = timeStep/timeBase # conversion factor power to energy for one time step
        self.global_cap_active = False

        #initialize jobs
        self.n = len(instanceData) #number of jobs
        self.jobs = [j for j in range(0,self.n)] 
        self.jobs_cap = [22 if instanceData["average_power_W"].iloc[j]>instanceData["maxPower"].iloc[j] else instanceData["maxPower"].iloc[j]/1000 for j in self.jobs] # in W
        self.jobs_demand = [instanceData["total_energy_Wh"].iloc[j]/1000 for j in self.jobs]  #in kWh

        #initialize intervals
        self.breakpoints = sorted(instanceData["t0_"+str(timeStep)].append(instanceData["t1_" + str(timeStep)]).drop_duplicates().tolist())
        self.intervals_start = self.breakpoints[:-1]
        self.intervals_end = self.breakpoints[1:]
        self.m = len(self.intervals_start) #number of atomic intervals
        self.I_a = [i for i in range(0,self.m)] 
        self.len_i = [(self.intervals_end[i] - self.intervals_start[i])*timeStep for i in self.I_a]

        self.find_J_inverse()
        self.find_J()

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
        jobMask = [(self.jobs_demand[j] > err) for j in range(0,len(self.jobs))]

        # update jobs
        self.jobs_cap = [self.jobs_cap[j] for j in range(0,len(self.jobs)) if jobMask[j]]
        self.jobs_demand = [self.jobs_demand[j] for j in range(0,len(self.jobs)) if jobMask[j]]
        self.jobs = [self.jobs[j] for j in range(0,len(self.jobs)) if jobMask[j]]
        self.n = len(self.jobs)

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
                    temp = self.intervals_start.index(self.data["t0_"+str(self.timeStep)].iloc[j])
                except:
                    if self.data["t0_"+str(self.timeStep)].iloc[j] < min(self.intervals_start):
                        temp = 0
                    else: 
                        print('[WARNING]: Unconsidered case. Check find_J_inverse() function in FOCSinstance class.')
                self.J_inverse["j"+str(j)] = [i for i in range(temp, self.intervals_end.index(self.data["t1_"+str(self.timeStep)].iloc[j])+1)] 
        return

    def find_J(self):
        self.J = {}
        for i in self.I_a:
            self.J["i"+str(i)] = [j for j in self.jobs if i in self.J_inverse["j"+str(j)]]
        return
    
    def add_capacity(self,cap, constant = True, activate = True):
        #define list with global maximum powers per interval
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

class FlowNet:
    def __init__(self):
        self.G = nx.DiGraph()

    def focs_instance_to_network(self, instance):
        #add nodes
        self.G.add_nodes_from(["s"] + ["j"+str(j) for j in instance.jobs] + ["i"+str(i) for i in instance.I_a] + ["t"])

        #add edges
        #later, initialize without edges to t and add those based on known capacity
        #can we use capacity function?
        self.G.add_edges_from(
                            [("s", "j"+str(j), {"capacity": instance.jobs_demand[instance.jobs.index(j)]}) for j in instance.jobs]   #D_s
                        +  sum([[("j"+str(j), "i"+str(i), {"capacity": instance.jobs_cap[instance.jobs.index(j)]*instance.len_i[i]/instance.timeBase}) for i in instance.J_inverse["j"+str(j)]] for j in instance.jobs],[])   #D_0
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

    def partial_flow(self,full_flow, node_set, intervals, jobs, J):
        partial_flow = copy.deepcopy(full_flow)
        for i in intervals:
            if i not in node_set:
                partial_flow["i{}".format(i)]["t"] = 0
                for j in J["i" + str(i)]:
                    #FIXME double check definition J[i]. Need to update function if we remove nodes... Active jobs only.
                    partial_flow["j{}".format(j)]["i{}".format(i)] = 0
        for j in jobs:
            partial_flow["s"]["j{}".format(j)] = sum(partial_flow["j{}".format(j)].values())
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
    

