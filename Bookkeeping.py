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



class Bookkeeping():
    def __init__(self):
        # empty lists to save results in
        self.prefix = []

        # model building runtimes. Will be a list of lists
        self.mbLP = []
        self.mbFOCS = []
        self.mbFOCSmpc = []

        # solving runtimes. Will be a list of lists
        self.solLP = []
        self.solFOCS = []
        self.solFOCSmpc = []

        # total runtimes. Will be a list of lists
        self.rtLP = []
        self.rtFOCS = []
        self.rtFOCSmpc = []

        # model building runtimes. Will be a list of lists
        self.mbLP_median = []
        self.mbFOCS_median = []
        self.mbFOCSmpc_median = []

        # solving runtimes. Will be a list of lists
        self.solLP_median = []
        self.solFOCS_median = []
        self.solFOCSmpc_median = []

        # total runtimes. Will be a list of lists
        self.rtLP_median = []
        self.rtFOCS_median = []
        self.rtFOCSmpc_median = []

        # objective values
        self.ovLP = []
        self.ovFOCS = []
        self.ovAVR = []
        self.ovOA = []
        self.ovFOCSmpc = []
        self.ovRatio = []
        self.ovRatioAVR = []
        self.ovRatioOA = []

        # power profiles
        self.power_profiles = pd.DataFrame({'i': []})

        ''' for remaining problem case '''
        # model building runtimes. Will be a list of lists
        self.mbLPred = []
        self.mbFOCSred = []
        self.mbFOCSmpcred = []

        # solving runtimes. Will be a list of lists
        self.solLPred = []
        self.solFOCSred = []
        self.solFOCSmpcred = []

        # total runtimes. Will be a list of lists
        self.rtLPred = []
        self.rtFOCSred = []
        self.rtFOCSmpcred = []

        # model building runtimes. Will be a list of lists
        self.mbLPred_median = []
        self.mbFOCSred_median = []
        self.mbFOCSmpcred_median = []

        # solving runtimes. Will be a list of lists
        self.solLPred_median = []
        self.solFOCSred_median = []
        self.solFOCSmpcred_median = []

        # total runtimes. Will be a list of lists
        self.rtLPred_median = []
        self.rtFOCSred_median = []
        self.rtFOCSmpcred_median = []

        # objective values
        self.ovLPred = []
        self.ovFOCSred = []
        self.ovFOCSmpcred = []
        self.ovRatiored = []
    
    def empty_temp(self):
        # model building runtimes. Will be a list of lists
        self.mbLPtemp = []
        self.mbFOCStemp = []
        self.mbFOCSmpctemp = []

        # solving runtimes. Will be a list of lists
        self.solLPtemp = []
        self.solFOCStemp = []
        self.solFOCSmpctemp = []

        # total runtimes. Will be a list of lists
        self.rtLPtemp = []
        self.rtFOCStemp = []
        self.rtFOCSmpctemp = []

        ## for partial problem
        # model building runtimes. Will be a list of lists
        self.mbLPtempred = []
        self.mbFOCStempred = []
        self.mbFOCSmpctempred = []

        # solving runtimes. Will be a list of lists
        self.solLPtempred = []
        self.solFOCStempred = []
        self.solFOCSmpctempred = []

        # total runtimes. Will be a list of lists
        self.rtLPtempred = []
        self.rtFOCStempred = []
        self.rtFOCSmpctempred = []
    
    def read_temp(self, path):
        data = pd.read_csv(path, delimiter=';')

        # model building runtimes. Will be a list of lists
        self.mbLPtemp += data['mbLP'].tolist()
        self.mbFOCStemp += data['mbFOCS'].tolist()
        self.mbFOCSmpctemp += data['mbFOCSmpc'].tolist()

        # solving runtimes. Will be a list of lists
        self.solLPtemp += data['solLP'].tolist()
        self.solFOCStemp += data['solFOCS'].tolist()
        self.solFOCSmpctemp += data['solFOCSmpc'].tolist()

        # total runtimes. Will be a list of lists
        self.rtLPtemp += data['rtLP'].tolist()
        self.rtFOCStemp += data['rtFOCS'].tolist()
        self.rtFOCSmpctemp += data['rtFOCSmpc'].tolist()

    def read_tempred(self, path):
        data = pd.read_csv(path, delimiter=';')

        # model building runtimes. Will be a list of lists
        self.mbLPtempred += data['mbLP'].tolist()
        self.mbFOCStempred += data['mbFOCS'].tolist()
        self.mbFOCSmpctempred += data['mbFOCSmpc'].tolist()

        # solving runtimes. Will be a list of lists
        self.solLPtempred += data['solLP'].tolist()
        self.solFOCStempred += data['solFOCS'].tolist()
        self.solFOCSmpctempred += data['solFOCSmpc'].tolist()

        # total runtimes. Will be a list of lists
        self.rtLPtempred += data['rtLP'].tolist()
        self.rtFOCStempred += data['rtFOCS'].tolist()
        self.rtFOCSmpctempred += data['rtFOCSmpc'].tolist()

    def write_runtimes_full_problem(self, focs=None, lp=None, focsmpc=None):

        # model building runtimes. Will be a list of lists
        self.mbLP += [self.mbLPtemp]
        self.mbFOCS += [self.mbFOCStemp]
        self.mbFOCSmpc += [self.mbFOCSmpctemp]

        # solving runtimes. Will be a list of lists
        self.solLP += [self.solLPtemp]
        self.solFOCS += [self.solFOCStemp]
        self.solFOCSmpc += [self.solFOCSmpctemp]

        # total runtimes. Will be a list of lists
        self.rtLP += [self.rtLPtemp]
        self.rtFOCS += [self.rtFOCStemp]
        self.rtFOCSmpc += [self.rtFOCSmpctemp]

        # model building runtimes. Will be a list of lists
        self.mbLP_median += [statistics.median(self.mbLPtemp)]
        self.mbFOCS_median += [statistics.median(self.mbFOCStemp)]
        self.mbFOCSmpc_median += [statistics.median(self.mbFOCSmpctemp)]

        # solving runtimes. Will be a list of lists
        self.solLP_median += [statistics.median(self.solLPtemp)]
        self.solFOCS_median += [statistics.median(self.solFOCStemp)]
        self.solFOCSmpc_median += [statistics.median(self.solFOCSmpctemp)]

        # total runtimes. Will be a list of lists
        self.rtLP_median += [statistics.median(self.rtLPtemp)]
        self.rtFOCS_median += [statistics.median(self.rtFOCStemp)]
        self.rtFOCSmpc_median += [statistics.median(self.rtFOCSmpctemp)]

        # objective values
        if focs is not None:
            self.ovFOCS += [focs.objective()]
        if lp is not None:
            self.ovLP += [lp.model.ObjVal]
            if focs is not None:
                self.ovRatio += [lp.model.ObjVal/focs.objective()]
        if focsmpc is not None:
            self.ovFOCSmpc +=[focsmpc.objective()]
        return

    def write_runtimes_partial_problem(self, focs=None, lp=None, focsmpc=None):

        # model building runtimes. Will be a list of lists
        self.mbLPred += [self.mbLPtempred]
        self.mbFOCSred += [self.mbFOCStempred]
        self.mbFOCSmpcred += [self.mbFOCSmpctempred]

        # solving runtimes. Will be a list of lists
        self.solLPred += [self.solLPtempred]
        self.solFOCSred += [self.solFOCStempred]
        self.solFOCSmpcred += [self.solFOCSmpctempred]

        # total runtimes. Will be a list of lists
        self.rtLPred += [self.rtLPtempred]
        self.rtFOCSred += [self.rtFOCStempred]
        self.rtFOCSmpcred += [self.rtFOCSmpctempred]

        # model building runtimes. Will be a list of lists
        self.mbLPred_median += [statistics.median(self.mbLPtempred)]
        self.mbFOCSred_median += [statistics.median(self.mbFOCStempred)]
        self.mbFOCSmpcred_median += [statistics.median(self.mbFOCSmpctempred)]

        # solving runtimes. Will be a list of lists
        self.solLPred_median += [statistics.median(self.solLPtempred)]
        self.solFOCSred_median += [statistics.median(self.solFOCStempred)]
        self.solFOCSmpcred_median += [statistics.median(self.solFOCSmpctempred)]

        # total runtimes. Will be a list of lists
        self.rtLPred_median += [statistics.median(self.rtLPtempred)]
        self.rtFOCSred_median += [statistics.median(self.rtFOCStempred)]
        self.rtFOCSmpcred_median += [statistics.median(self.rtFOCSmpctempred)]

        # objective values
        if focs is not None:
            self.ovFOCSred += [focs.objective()]
        if lp is not None:
            self.ovLPred += [lp.model.ObjVal]
            if focs is not None:
                self.ovRatiored += [lp.model.ObjVal/focs.objective()]
        if focsmpc is not None:
            self.ovFOCSmpcred +=[focsmpc.objective()]
        return
    
    def write_cf(self, focs, avr, oa):
        if avr is not None:
            self.ovAVR += [avr.objective()]

        if oa is not None:
            self.ovOA += [oa.objective()]

        if focs is not None:
            self.ovFOCS += [focs.objective()]
        
        if (focs is not None) and (avr is not None):
            self.ovRatioAVR += [avr.objective()/focs.objective()]

        if (focs is not None) and (oa is not None):
            self.ovRatioOA += [oa.objective()/focs.objective()]

        #load duration is later only a sorted list of this. Since x1 and x2 are 'independent' in the sense that we don't want them in the same plot.
        return
    
    def write_power(self,focs,avr,oa):
        if len(self.prefix) > 0:
            header = ['i', 'focs_'+str(self.prefix[-1]), 'avr_'+str(self.prefix[-1]), 'oa_'+str(self.prefix[-1])]
        else:
            header = ['i', 'focs_', 'avr_', 'oa_']
        
        # assuming all results are for the same instance and time horizon
        i_timeStep = [focs.instance.intervals_start[i]+freq for i in range(0, focs.instance.m) for freq in range(0,int(focs.instance.len_i[i]/focs.instance.timeStep))]
        power_profiles = [i_timeStep]
        for alg in [focs, avr, oa]:
            temp = []
            # determine powers
            p_i = [((alg.f["i"+str(i)]['t']/(alg.instance.len_i[i]))*alg.instance.timeBase) for i in range(0,alg.instance.m)]
            for i in range(0,alg.instance.m):
                # extrapolate
                temp += [p_i[i] for t in range(0,int(alg.instance.len_i[i]/alg.instance.timeStep))]
            power_profiles += [temp]
        
        # make dataframe
        df = pd.DataFrame(power_profiles).T
        df.columns = header 

        # merge with other dataframe
        self.power_profiles = pd.merge(self.power_profiles, df, on = 'i', how = 'outer')

        return
    
    def write_online_to_csv(self):
        path = 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data/output/online/'
        #write power
        self.power_profiles.to_csv(path + 'power_profiles_{}.csv'.format(self.prefix[-1].rsplit('_', 1)[0]), sep=';')

        index = [self.prefix[x].split("_")[-1] for x in range(0,len(self.prefix))]

        #write competitive ratios
        with open(path + 'cf_{}.csv'.format(self.prefix[-1].rsplit('_', 1)[0]), 'w', newline='' ) as f:
            writer = csv.writer(f, delimiter=';')
            #header
            writer.writerow(['index', 'ovFOCS', 'ovAVR', 'ovOA', 'ovRatioAVR', 'ovRatioOA'])
            #content
            writer.writerows(np.array([index, self.ovFOCS, self.ovAVR, self.ovOA, self.ovRatioAVR, self.ovRatioOA]).T.tolist())
        
        #write competitive ratios sorted
        with open(path + 'cf_{}_sorted.csv'.format(self.prefix[-1].rsplit('_', 1)[0]), 'w', newline='' ) as f:
            writer = csv.writer(f, delimiter=';')
            #header
            writer.writerow(['index', 'ovFOCS', 'ovAVR', 'ovOA', 'ovRatioAVR', 'ovRatioOA'])
            #content
            writer.writerows(np.array([index, sorted(self.ovFOCS), sorted(self.ovAVR), sorted(self.ovOA), sorted(self.ovRatioAVR), sorted(self.ovRatioOA)]).T.tolist())
        
        return

    def write_instance_to_csv(self, pf = ''):
        path = 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data/output/instances/'
        with open(path + 'instance_{}{}.csv'.format(self.prefix[-1], pf), 'w', newline='' ) as f:
            writer = csv.writer(f, delimiter=';')
            #header
            writer.writerow(['mbLP', 'mbFOCS', 'mbFOCSmpc', 'solLP', 'solFOCS', 'solFOCSmpc', 'rtLP', 'rtFOCS', 'rtFOCSmpc'])
            #content
            if pf == 'partial':
                writer.writerows(np.array([self.mbLPtempred, self.mbFOCStempred, self.mbFOCSmpctempred, self.solLPtempred, self.solFOCStempred, self.solFOCSmpctempred, self.rtLPtempred, self.rtFOCStempred, self.rtFOCSmpctempred]).T.tolist())
            else: 
                writer.writerows(np.array([self.mbLPtemp, self.mbFOCStemp, self.mbFOCSmpctemp, self.solLPtemp, self.solFOCStemp, self.solFOCSmpctemp, self.rtLPtemp, self.rtFOCStemp, self.rtFOCSmpctemp]).T.tolist())
            return


    def write_timestep_to_csv(self, fc = None, pathOR = None):
        path = 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data/output/timesteps/'
        if pathOR is not None:
            path = pathOR
        #full
        pf = 'full'
        with open(path + 'timeStep_{}_{}.csv'.format(self.prefix[-1].rsplit('_', 2)[0], pf), 'w', newline='' ) as f:
            writer = csv.writer(f, delimiter=';')
            if fc is not None:
                n = len(fc)
                #header
                writer.writerow(['n', 'mbLP', 'mbFOCS', 'mbFOCSmpc', 'solLP', 'solFOCS', 'solFOCSmpc', 'rtLP', 'rtFOCS', 'rtFOCSmpc'])
                #content
                writer.writerows(np.array([fc, self.mbLP_median[-n:], self.mbFOCS_median[-n:], self.mbFOCSmpc_median[-n:], self.solLP_median[-n:], self.solFOCS_median[-n:], self.solFOCSmpc_median[-n:], self.rtLP_median[-n:], self.rtFOCS_median[-n:], self.rtFOCSmpc_median[-n:]]).T.tolist())
            else:
                #header
                writer.writerow(['mbLP', 'mbFOCS', 'mbFOCSmpc', 'solLP', 'solFOCS', 'solFOCSmpc', 'rtLP', 'rtFOCS', 'rtFOCSmpc'])
                #content
                writer.writerows(np.array([self.mbLP_median[-n:], self.mbFOCS_median[-n:], self.mbFOCSmpc_median[-n:], self.solLP_median[-n:], self.solFOCS_median[-n:], self.solFOCSmpc_median[-n:], self.rtLP_median[-n:], self.rtFOCS_median[-n:], self.rtFOCSmpc_median[-n:]]).T.tolist())
        #partial
        pf = 'partial'
        with open(path + 'timeStep_{}_{}.csv'.format(self.prefix[-1].rsplit('_', 2)[0], pf), 'w', newline='' ) as f:
            writer = csv.writer(f, delimiter=';')
            if fc is not None:
                n = len(fc)
                #header
                writer.writerow(['n', 'mbLP', 'mbFOCS', 'mbFOCSmpc', 'solLP', 'solFOCS', 'solFOCSmpc', 'rtLP', 'rtFOCS', 'rtFOCSmpc'])
                #content
                writer.writerows(np.array([fc, self.mbLPred_median[-n:], self.mbFOCSred_median[-n:], self.mbFOCSmpcred_median[-n:], self.solLPred_median[-n:], self.solFOCSred_median[-n:], self.solFOCSmpcred_median[-n:], self.rtLPred_median[-n:], self.rtFOCSred_median[-n:], self.rtFOCSmpcred_median[-n:]]).T.tolist())
            else:
                #header
                writer.writerow(['mbLP', 'mbFOCS', 'mbFOCSmpc', 'solLP', 'solFOCS', 'solFOCSmpc', 'rtLP', 'rtFOCS', 'rtFOCSmpc'])
                #content
                writer.writerows(np.array([self.mbLPred_median[-n:], self.mbFOCSred_median[-n:], self.mbFOCSmpcred_median[-n:], self.solLPred_median[-n:], self.solFOCSred_median[-n:], self.solFOCSmpcred_median[-n:], self.rtLPred_median[-n:], self.rtFOCSred_median[-n:], self.rtFOCSmpcred_median[-n:]]).T.tolist())
        return
    
    def write_flowmethod_to_csv(self, fc = None, sc = None):
        path = 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data/output/flowmethod/'
        #full
        pf = 'full'
        with open(path + 'flowmethod_{}_{}.csv'.format(self.prefix[-1].rsplit('_')[0], pf), 'w', newline='' ) as f:
            writer = csv.writer(f, delimiter=';')
            if (fc is not None) & (sc is not None):
                n = len(fc)*len(sc)
                fc2 = fc*len(sc)
                sc2 = [j for j in sc for i in range(0,len(fc))]
                #header
                writer.writerow(['n', 'timestep', 'mbLP', 'mbFOCS', 'mbFOCSmpc', 'solLP', 'solFOCS', 'solFOCSmpc', 'rtLP', 'rtFOCS', 'rtFOCSmpc'])
                #content
                writer.writerows(np.array([fc2, sc2, self.mbLP_median[-n:], self.mbFOCS_median[-n:], self.mbFOCSmpc_median[-n:], self.solLP_median[-n:], self.solFOCS_median[-n:], self.solFOCSmpc_median[-n:], self.rtLP_median[-n:], self.rtFOCS_median[-n:], self.rtFOCSmpc_median[-n:]]).T.tolist())
            else:
                #header
                writer.writerow(['mbLP', 'mbFOCS', 'mbFOCSmpc', 'solLP', 'solFOCS', 'solFOCSmpc', 'rtLP', 'rtFOCS', 'rtFOCSmpc'])
                #content
                writer.writerows(np.array([self.mbLP_median[-n:], self.mbFOCS_median[-n:], self.mbFOCSmpc_median[-n:], self.solLP_median[-n:], self.solFOCS_median[-n:], self.solFOCSmpc_median[-n:], self.rtLP_median[-n:], self.rtFOCS_median[-n:], self.rtFOCSmpc_median[-n:]]).T.tolist())
        #partial
        pf = 'partial'
        with open(path + 'flowmethod_{}_{}.csv'.format(self.prefix[-1].rsplit('_')[0], pf), 'w', newline='' ) as f:
            writer = csv.writer(f, delimiter=';')
            if (fc is not None) & (sc is not None):
                n = len(fc)*len(sc)
                fc2 = fc*len(sc)
                sc2 = [j for j in sc for i in range(0,len(fc))]
                #header
                writer.writerow(['n', 'timestep', 'mbLP', 'mbFOCS', 'mbFOCSmpc', 'solLP', 'solFOCS', 'solFOCSmpc', 'rtLP', 'rtFOCS', 'rtFOCSmpc'])
                #content
                writer.writerows(np.array([fc2, sc2, self.mbLPred_median[-n:], self.mbFOCSred_median[-n:], self.mbFOCSmpcred_median[-n:], self.solLPred_median[-n:], self.solFOCSred_median[-n:], self.solFOCSmpcred_median[-n:], self.rtLPred_median[-n:], self.rtFOCSred_median[-n:], self.rtFOCSmpcred_median[-n:]]).T.tolist())
            else:
                #header
                writer.writerow(['mbLP', 'mbFOCS', 'mbFOCSmpc', 'solLP', 'solFOCS', 'solFOCSmpc', 'rtLP', 'rtFOCS', 'rtFOCSmpc'])
                #content
                writer.writerows(np.array([self.mbLPred_median[-n:], self.mbFOCSred_median[-n:], self.mbFOCSmpcred_median[-n:], self.solLPred_median[-n:], self.solFOCSred_median[-n:], self.solFOCSmpcred_median[-n:], self.rtLPred_median[-n:], self.rtFOCSred_median[-n:], self.rtFOCSmpcred_median[-n:]]).T.tolist())
        return
