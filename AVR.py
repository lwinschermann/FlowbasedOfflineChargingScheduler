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

class AVR:
    def __init__(self, instance, flowNet, flowOp):
        self.instance = instance
        self.flowNet = flowNet
        self.flowOp = flowOp
        self.f = flowOp.empty_flow(flowOp.G)
    
    def solve_avr(self):
        avr_power = []
        for job in self.instance.jobs:
            #length of stay job j in [units of timeBase]. Usually hours
            dwell_j = self.flowOp.length_sum_intervals(self.instance.J_inverse['j'+str(job)], self.instance.len_i)/self.instance.timeBase

            #calculate average power
            avr_power_j = self.instance.jobs_demand[job]/dwell_j

            #loop over i 
            for i in self.instance.J_inverse['j'+str(job)]:
                #calculate energy per interval connection
                energy_i_j = avr_power_j * self.instance.len_i[i]/self.instance.timeBase

                #add j-i
                self.f['j'+str(job)]['i'+str(i)] = energy_i_j
            
                #ADD i-t
                self.f['i'+str(i)]['t'] += energy_i_j
            
            #add s-j
            self.f['s']['j'+str(job)] = sum([self.f['j'+str(job)]['i'+str(i)] for i in self.instance.J_inverse['j'+str(job)]])
        
            avr_power = avr_power + [avr_power_j]
        self.avr_power = avr_power
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