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

from networkx.algorithms.flow import shortest_augmenting_path
import pandas as pd

import random

from FOCS import FOCS, FlowNet, FlowOperations, FOCSinstance
from Bookkeeping import Bookkeeping
from AVR import AVR
from OA import OA
from Price import Price

bk = Bookkeeping()
reps = 500
instanceSizes = [400]
timeSteps = [900]
maxFlows = [shortest_augmenting_path]
write = False
randomSample = True

# Real Training data
instanceData = pd.read_excel('data/input/filteredData.xlsx')

for timeStep in timeSteps:
    for instanceSize in instanceSizes:
        for rep in range(0,reps):
            print('rep = ', rep)
            bk.prefix += ['rev_' + str(0) + '_' + str(timeStep) + '_' + str(instanceSize) + '_' + str(rep)]
            '''--------------start instance with capacity--------------'''
            if randomSample is True:
                sample = sorted(random.sample(range(0,len(instanceData)), instanceSize))
                instance = FOCSinstance(instanceData.iloc[sample], timeStep)
                instance_periodic = FOCSinstance(instanceData.iloc[sample], timeStep, periodicity=True)  
            else:
                instance = FOCSinstance(instanceData[:instanceSize], timeStep=timeStep)      
                instance_periodic = FOCSinstance(instanceData[:instanceSize], timeStep=timeStep, periodicity=True)      
            '''--------------start Greedy----------'''
            flowNet = FlowNet()
            flowNet.focs_instance_to_network(instance_periodic)
            flowOp = FlowOperations(flowNet.G, instance_periodic)

            greedy = Price(instance_periodic, flowNet, flowOp, [i for i in range(0,instance_periodic.breakpoints[-1]-instance_periodic.breakpoints[0])])
            greedy.solve_cost()

            '''--------------start OA--------------'''
            flowNet = FlowNet()
            flowNet.focs_instance_to_network(instance)
            flowOp = FlowOperations(flowNet.G, instance)
            
            oa = OA(instance, flowNet, flowOp)

            f = oa.solve_oa()

            '''--------------start AVR--------------'''
            avr = AVR(instance, flowNet, flowOp)

            f = avr.solve_avr()

            '''--------------start FOCS--------------'''
            focs = FOCS(instance, flowNet, flowOp)

            f = focs.solve_focs()

            '''--------------write results of rep----------'''
            bk.write_cf(focs, avr, oa, greedy)
            bk.write_power(focs, avr, oa, greedy)

        if write is True:
            bk.write_online_to_csv()

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