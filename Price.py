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

class Price:
    def __init__(self, instance, flowNet, flowOp, co2):
        self.instance = instance
        self.flowNet = flowNet
        self.flowOp = flowOp
        self.co2 = co2
        self.f = flowOp.empty_flow(flowOp.G)
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

        self.co2_unique = list(set(sorted(self.co2)))
        self.co2_unique_index = [i for i in range(0,len(self.co2_unique))]
        self.co2_inverse_index = self.find_co2_inverse()

    def find_co2_inverse(self):
        co2_inverse_index = {}
        for c in self.co2_unique_index:
            co2_inverse_index['c'+str(c)] = [i for i in self.I_a if self.co2[i] == self.co2_unique[c]]
        return co2_inverse_index

    def solve_cost(self, err=0.0000001):
        self.G_cost = copy.deepcopy(self.flowOp.G)

        print('[WARNING]: Weighting in min cost flow does not normalize for interval length! Use unit-sized intervals.')

        # make network with capacities based on cost
        for i in self.instance.I_a:
            self.G_cost['i'+str(i)]['t']['weight'] = self.co2[i]
            pass

        # if applicable, add global capacity constraints
        if self.instance.global_cap_active:
            for i in self.instance.I_a:
                self.G_cost['i'+str(i)]['t']['capacity'] = self.instance.global_cap[i]*self.instance.len_i[i]*self.instance.tau/self.instance.timeStep

        # add demands
        demand = sum([self.G_cost['s']['j'+str(j)]['capacity'] for j in self.instance.jobs]) - err
        nx.set_node_attributes(self.G_cost, {'s': - demand, 't': demand}, name = 'demand')

        # solve min cost flow
        self.f_cost = nx.min_cost_flow(self.G_cost)

        # extract list with energy per co2 value
        co2_cap = {}
        co2_idx = [i for i in range(0,len(self.co2))]
        used_intervals = []
        for uniqueVal in range(0,len(set(self.co2))):
            # which co2idx apply
            cheap_intervals = [co2_idx[i] for i in co2_idx if self.co2[i] == min([self.co2[k] for k in co2_idx if co2_idx[k] not in used_intervals])]
            used_intervals += cheap_intervals
            # total energy provided with that co2 cost
            co2_cap['c'+str(uniqueVal)] = sum([self.f_cost['i'+str(i)]['t'] for i in cheap_intervals])
        self.co2_cap = co2_cap

        return co2_cap
    
    def augment_instance(self, dummy_profile, err=0.000001):
        augmented_instance = copy.deepcopy(self.instance)

        augmented_instance.dummy_indices = [i for i in range(0,len(dummy_profile)) if dummy_profile[i] > 0]

        #add jobs 
        augmented_instance.n_dummy = len(augmented_instance.dummy_indices)
        augmented_instance.jobs += [augmented_instance.n + i for i in range(0,augmented_instance.n_dummy)] 
        augmented_instance.jobs_cap += [dummy_profile[augmented_instance.dummy_indices[i]]*(1+1/self.instance.tau) for i in range(0,augmented_instance.n_dummy)] # in W. Deliberately made greater than it is for feasibility
        augmented_instance.jobs_demand += [dummy_profile[augmented_instance.dummy_indices[i]] for i in range(0,augmented_instance.n_dummy)]  #in kWh

        # not updating intervals since we assume periodicity == True
        #update intervals
        
        for i in range(0,len(augmented_instance.dummy_indices)):
            augmented_instance.J_inverse['j'+str(augmented_instance.n + i)] = [augmented_instance.dummy_indices[i]]
            augmented_instance.J['i'+str(augmented_instance.dummy_indices[i])] += [augmented_instance.n + i]

        #update global capacity
        if augmented_instance.global_cap_active:
            for i in augmented_instance.dummy_indices:
                augmented_instance.global_cap[i] += dummy_profile[i]*(err+1/self.instance.tau) 
            # augmented_instance.global_cap = [augmented_instance.global_cap[i] + dummy_profile[i]*(err+1/self.instance.tau) for i in augmented_instance.I_a]

        augmented_instance.n = len(augmented_instance.jobs) #number of jobs

        return augmented_instance
    
    def solve_flatten_cost(self, weight_c = 1, weight_f = 0.0001):
        # make profile based on co2 and weights
        weight_term = weight_c/(2*weight_f)
        co2_dummy_profile = [weight_term*self.co2[i] for i in self.instance.I_a]
      
        # augment instance with price weights
        augmented_instance = self.augment_instance(co2_dummy_profile)

        # solve for augmented instance
        flowNet = FlowNet()
        flowNet.focs_instance_to_network(augmented_instance)
        flowOp = FlowOperations(flowNet.G, augmented_instance)
        focs = FOCS(augmented_instance, flowNet, flowOp)
        f = focs.solve_focs()

        # transform into un-augmented instance
        for j in self.instance.jobs:
            self.f['s']['j'+str(j)] = f['s']['j'+str(j)]
            for i in self.instance.J_inverse['j'+str(j)]:
                self.f['j'+str(j)]['i'+str(i)] = f['j'+str(j)]['i'+str(i)]
        for i in self.instance.I_a:
            self.f['i'+str(i)]['t'] = sum([self.f['j'+str(j)]['i'+str(i)] for j in self.instance.J['i'+str(i)]])

        return self.f
     
    #copied from FOCS class
    def objective(self, f=None, structure='focs', co2obj = False, priceobj = None, partial=[]):
        if f is None:
            f = self.f
        if len(partial)==2:
            # only take interval between partial[0] (including) and partial[1] (excluding)
            node_set = self.instance.I_a[self.instance.I_a.index(partial[0]):self.instance.I_a.index(partial[1])]
            # reduce flow
            f = self.flowOp.partial_flow(f, node_set, self.instance.I_a, self.instance.jobs, self.instance.J)
        else:
            node_set = self.instance.I_a
        #normalized objective for validation
        #instead of the sum of squares, we first weight the square of the power (e/len_i[i]) by len_i[i]/timeStep      
        #determine power per interval
        if structure == 'focs':
            p_i = [((f["i"+str(i)]['t']/(self.instance.len_i[i]))*self.instance.timeBase) for i in range(0,self.instance.m)] 
            if co2obj:
                #assuming co2 in g CO2eq/kWh
                energyPriced = [self.co2_unique[c]*sum([f['i'+str(i)]['t'] for i in self.co2_inverse_index['c'+str(c)]]) for c in self.co2_unique_index]
                self.objCO2 = sum(energyPriced)
                return self.objCO2
            if priceobj is not None:
                energyPriced = [priceobj[i]*f['i'+str(i)]['t'] for i in node_set]
                self.objPrice = sum(energyPriced)
                return self.objPrice
        elif structure == 'price':
            p_i = [((f["i"+str(i)]["c"+str(self.co2_unique.index(self.co2[i]))]/(self.instance.len_i[i]))*self.instance.timeBase) for i in range(0,self.instance.m)]
            if co2obj:
                #assuming co2 in g CO2eq/kWh
                energyPriced = [self.co2[c]*f['c'+str(c)]['t'] for c in self.co2_unique_index]
                self.objCO2 = sum(energyPriced)
                return self.objCO2        
        else:
            print('[WARNING]: Structure for flow undefined. Try focs or price.')
        #determine weighted squared powers
        powerSquare = [(p_i[i]**2)*(self.instance.len_i[i]/self.instance.timeStep) for i in range(0,self.instance.m)]
        self.objNormalized = sum(powerSquare)

        return self.objNormalized
    

# Note: This class is a heuristic. 
# You greedily choose the cheapest intervals, flatten them out, SAVE A PARTIAL SCHEDULE, then repeat. By saving the partial schedule (without later rescheduling) you loose optimality.
class Greedy:
    def __init__(self, instance, flowNet, flowOp, co2):
        self.instance = instance
        self.flowNet = flowNet
        self.flowOp = flowOp
        self.co2 = co2
        self.f = flowOp.empty_flow(flowOp.G)

    def solve_greedy(self, err = 0.000000001):
        co2_idx = [i for i in range(0,len(self.co2))]
        used_intervals = []

        #initialize full copy instance
        self.instance.remaining_instance = copy.deepcopy(self.instance)

        #initialize partial instance
        self.instance.empty_instance()

        # loop through unique values in co2:
        for uniqueVal in range(0,len(set(self.co2))):
            
            # choose cheapest intervals
            cheap_intervals = [co2_idx[i] for i in co2_idx if self.co2[i] == min([self.co2[k] for k in co2_idx if co2_idx[k] not in used_intervals])]
            print('cheap intervals = ', cheap_intervals)
            used_intervals += cheap_intervals

            # make network
            flowNet = FlowNet()
            flowNet.focs_instance_to_network(self.instance.remaining_instance)

            # reduce such that only consider cheapest intervals
            for i in co2_idx:
                if i in cheap_intervals:
                    #if global maximum power constraint is active, trigger it here
                    if self.instance.global_cap_active:
                        flowNet.G["i{}".format(i)]["t"]["capacity"] = self.instance.global_cap[i]*self.instance.len_i[i]*self.instance.tau/self.instance.timeStep  
                else:
                    flowNet.G.remove_node("i{}".format(i))

            flowOp = FlowOperations(flowNet.G, self.instance.remaining_instance)
            G_cheap = copy.deepcopy(flowNet.G) #focs.solve will alter network. We need it to add next time step to overall solution

            # solve max flow
            flow_val, f = nx.maximum_flow(flowNet.G, "s", "t", flow_func=shortest_augmenting_path)

            if len(cheap_intervals) > 1:
                # reduce instance to that max flow. 
                partial_instance = self.instance.feasible_subinstance(f, cheap_intervals)
                
                # initiate and solve FOCS for reduced instance
                flowNet = FlowNet()
                flowNet.focs_instance_to_network(self.instance.partial_instance)
                flowOp = FlowOperations(flowNet.G, self.instance.partial_instance)
                G_cheap = copy.deepcopy(flowNet.G) #focs.solve will alter network. We need it to add next time step to overall solution
                focs = FOCS(self.instance.partial_instance, flowNet, flowOp)
                focs.solve_focs()
                
                f = focs.f

            flowOp.add_flows(self.f, f, G_cheap)

            # reduce remaining instance
            self.instance.remaining_instance.jobs_demand = [self.instance.jobs_demand[j] - sum([f['j'+str(j)]['i'+str(i)] for i in self.instance.remaining_instance.J_inverse['j'+str(j)] if i in cheap_intervals]) for j in self.instance.remaining_instance.jobs]
    
            # check if done. 
            if sum(self.instance.jobs_demand) - sum([self.f["s"]["j{}".format(j)] for j in self.instance.jobs]) < err:
                print('doneeee')
                print('[WARNING]: Heuristic.')
                print('\n',self.f)
                return self.f
        print('[WARNING]: Looped through all co2 levels, but algorithm did not terminate.')
        #DONT FORGET PRICE OBJECTIVE FUNCTION
        return self.f

    #copied from FOCS class
    def objective(self, co2obj = False):
        if co2obj:
            #assuming co2 in g CO2eq/kWh
            energyPriced = [self.co2[c]*f['c'+str(c)]['t'] for c in self.co2_unique_index]
            self.objCO2 = sum(energyPriced)
            return self.objCO2
        
        #normalized objective for validation
        #instead of the sum of squares, we first weight the square of the power (e/len_i[i]) by len_i[i]/timeStep      
        #determine power per interval
        p_i = [((self.f["i"+str(i)]['t']/(self.instance.len_i[i]))*self.instance.timeBase) for i in range(0,self.instance.m)]
        #determine weighted squared powers
        powerSquare = [(p_i[i]**2)*(self.instance.len_i[i]/self.instance.timeStep) for i in range(0,self.instance.m)]
        self.objNormalized = sum(powerSquare)
        
        return self.objNormalized