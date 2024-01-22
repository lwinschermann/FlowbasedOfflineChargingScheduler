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

from gurobipy import *

class LP():
    def __init__(self, instance):
        self.instance = instance
    
    def build_model(self):
        model = Model("Create EV charging schedule minimizing aggregated power 2-norm")
        model.modelSense = GRB.MINIMIZE
        model.Params.OutputFlag = 0

        #What variables are given
        n = self.instance.n             # number of jobs
        T = self.instance.m             # time horizon
        jobs = self.instance.jobs
        len_i = [self.instance.len_i[i]/self.instance.timeStep for i in range(0,T)]    # interval lengths in units of timeStep. Default quarters.
        aj = [self.instance.J_inverse["j" + str(j)][0] for j in self.instance.jobs]    # arrival times per job
        dj = [self.instance.J_inverse["j" + str(j)][-1] + 1 for j in self.instance.jobs]   # departure times per job
        ej = self.instance.jobs_demand    # energy demand per job in kWh (h modelled by timeBase for conversion)
        pmaxj = self.instance.jobs_cap # maximum power per job in kW
        pj = []                        # empty list for decision variables by job. In kWh

        for jId, j in enumerate(jobs):
            # Decision variables 
            pj = pj + [model.addVars(tuplelist([(j,t) for t in self.instance.J_inverse["j"+str(j)]]), name = "decision variables pj", vtype = GRB.CONTINUOUS, lb = 0)] #in power units
            model.update()
            # add constraints
            # supply all energy. kWh == kW*units of timeStep*timeStep/timeBase ok.
            model.addConstr(ej[jId] == quicksum(pj[jId][j,t]*len_i[t]*self.instance.tau for t in self.instance.J_inverse["j" +str(j)]), name = 'energy job')

            # stay within power limit
            for t in self.instance.J_inverse["j"+ str(j)]:
                # kW <= kW ok.
                model.addConstr(pj[jId][j,t] <= pmaxj[jId], name = 'power limit j = ' +str(j))

        # objective function
        # sum(((sum kW)^2)*units of timeStep) --> kW^2
        model.setObjective(quicksum((quicksum(pj[jobs.index(j)][j,t] for j in self.instance.J['i' + str(t)])**2)*len_i[t] for t in range(0,T)), sense = GRB.MINIMIZE) # normalized for intervals length. 
        model.update()
        self.pj = pj
        self.model = model
        return self.model

    def solve_model(self):
        # optimize
        self.model.optimize()

    def print_results(self):
        # print resulting schedules
        for v in self.model.getVars():
            print("{}: {}".format(v.varName, v.X))
        # print objective value
        print('objective LP = ', self.model.ObjVal)


