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
        self.rtFULPES1 = []
        self.rtFULPES2 = []

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
        self.rtFULPES1_median = []
        self.rtFULPES2_median = []

        # objective values
        self.ovLP = []
        self.ovFOCS = []
        self.ovAVR = []
        self.ovOA = []
        self.ovGreedy = []
        self.ovFOCSmpc = []
        self.ovRatio = []
        self.ovRatioAVR = []
        self.ovRatioOA = []
        self.ovRatioGreedy = []

        # power profiles
        self.power_profiles = pd.DataFrame({'i': []})
        self.agg_power_profiles = []

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

        '''for schedule experiments to record qos and qoe measures'''
        # job specific measures
        self.schedHeadersJobs = ['idx', 'weighting_method', 'day', 'count so far', 'ens_abs_exact', 'ens_rel_exact', 'ens_abs', 'ens_rel', 'qos1', 'qos2_real', 'qos2_plan', 'qos3_real', 'qos3_plan', 'qoe_real', 'qoe_real_pot', 'qoe_plan', 'qoe_plan_pot']
        # global measures
        self.schedHeadersGlobal = ['weighting_method', 'day', 'qoe_real_total_exact', 'qoe_real_pot_total_exact', 'qoe_plan_total_exact', 'qoe_plan_pot_total_exact', 'qoe_real_total_rel', 'qoe_real_pot_total_rel', 'qoe_plan_total_rel', 'qoe_plan_pot_total_rel', 'ens_abs_max', 'ens_rel_max', 'jain_ens_rel', 'jain_ens_rel_exact', 'hossfeld_ens_rel', 'hossfeld_ens_rel_exact', 'es_exact', 'es', 'ens_abs_exact', 'ens_abs', 'ens_rel_exact_avg', 'ens_rel_avg', 'qos_1_min', 'qos_2_plan_min', 'qos_2_real_min', 'qos_3_plan_min', 'qos_3_real_min', 'jain_qos_1', 'jain_qos_2_plan', 'jain_qos_2_real', 'jain_qos_3_plan', 'jain_qos_3_real', 'hossfeld_qos_1', 'hossfeld_qos_2_plan', 'hossfeld_qos_2_real', 'hossfeld_qos_3_plan', 'hossfeld_qos_3_real']
        
        self.schedJobs = []
        self.schedGlobal = []

    def write_objectives_schedule(self, schedule, weighting_method, day, idx, n_idx):
        self.schedJobs += [[idx, [weighting_method for i in idx], [day for i in idx], n_idx, schedule.jobs_ens_abs_exact, schedule.jobs_ens_rel_exact, schedule.jobs_ens_abs, schedule.jobs_ens_rel, schedule.jobs_es_rel, schedule.jobs_qos2_waiting_real, schedule.jobs_qos2_waiting_plan, schedule.jobs_qos3_powervar_real, schedule.jobs_qos3_powervar_plan, schedule.jobs_qoe_real_ed, schedule.jobs_qoe_real_pot, schedule.jobs_qoe_plan_ed, schedule.jobs_qoe_plan_pot]]
        self.schedGlobal += [[weighting_method, day, schedule.qoe_real_ed_total_exact, schedule.qoe_real_pot_total_exact, schedule.qoe_plan_ed_total_exact, schedule.qoe_plan_pot_total_exact, schedule.qoe_real_ed_total_rel, schedule.qoe_real_pot_total_rel, schedule.qoe_plan_ed_total_rel, schedule.qoe_plan_pot_total_rel, schedule.jobs_ens_abs_max, schedule.jobs_ens_rel_max, schedule.jain_ens_rel, schedule.jain_ens_rel_exact, schedule.hossfeld_ens_rel, schedule.hossfeld_ens_rel_exact, schedule.es_exact, schedule.es, schedule.ens_abs_exact, schedule.ens_abs, schedule.ens_rel_exact_avg, schedule.ens_rel_avg, schedule.qos_1_min, schedule.qos_2_plan_min, schedule.qos_2_real_min, schedule.qos_3_plan_min, schedule.qos_3_real_min, schedule.jain_qos_1, schedule.jain_qos_2_plan, schedule.jain_qos_2_real, schedule.jain_qos_3_plan, schedule.jain_qos_3_real, schedule.hossfeld_qos_1, schedule.hossfeld_qos_2_plan, schedule.hossfeld_qos_2_real, schedule.hossfeld_qos_3_plan, schedule.hossfeld_qos_3_real]]
        return
    
    def write_objectives_schedule_to_csv(self, data, weighting_methods, days, jobs, suffix = ''):
        # make dataframe of schedglobal and schedjobs
        dfSchedJobs = pd.concat([pd.DataFrame(np.array(x_list).T.tolist(), columns = self.schedHeadersJobs) for x_list in self.schedJobs])
        dfSchedJobs.reset_index(inplace=True)
        
        dfSchedGlobal = pd.DataFrame(self.schedGlobal, columns=self.schedHeadersGlobal)

        # bk save global results
        dfSchedGlobal.to_csv("data/output/schedules/qosqoe/schedule_global_objectives.csv", sep=';')

        # bk concat bk.schedJobs to instanceData
        # create idx column in instanceData
        data['idx'] = data.index.astype(str)
        # merge bk.schedJobs into instanceData for each weighting method
        header = data.columns.values.tolist()
        scheds = ['s','s_lin','s_q','s_linp','s_qp','s_prob','s_probfull']
        for sched in scheds:
            temp_j = copy.deepcopy(dfSchedJobs[dfSchedJobs['weighting_method']==sched])
            temp_j = temp_j.set_index('idx')
            data = data.merge(temp_j, on = 'idx', how = 'left', suffixes = ['' ,'_'+sched])
            data = data.drop(['index', 'weighting_method'], axis=1)
            header += [x+'_'+sched for x in ['day','count so far_y','ens_abs_exact','ens_rel_exact','ens_abs','ens_rel','qos1','qos2_real','qos2_plan','qos3_real','qos3_plan','qoe_real','qoe_real_pot','qoe_plan','qoe_plan_pot']]
        # rename column headers
        data.columns = header

        # drop day entries
        data = data.drop(['day_'+x for x in scheds[1:]], axis=1)
        data = data.rename({'day_s':'day'}, axis=1)
        data_filtered = data.dropna(axis=0, subset = 'day')
        # bk save instanceData
        data.to_csv("data/output/schedules/qosqoe/instanceData_j_obj.csv", sep=';')
        data_filtered.to_csv("data/output/schedules/qosqoe/instanceData_j_obj_filtered.csv", sep=';')
        
        # # bk save results per car
        #FIXME maybe drop data that has not been used in sim?
        ids = data_filtered['card_id'].unique()
        for idx, id in enumerate(ids):
            # filter instanceData for job
            temp = copy.deepcopy(data_filtered[data_filtered['card_id']==id])
            # save to csv the last x entries and count so far
            temp[temp.columns[-14*8+12:]].to_csv("data/output/schedules/qosqoe/perEV/instanceData_EV{}_obj.csv".format(idx), sep=';')

        # bk save results per schedule per day
        for sched in weighting_methods:
            temp_s = dfSchedJobs[dfSchedJobs['weighting_method']==sched]
            temp_s.to_csv("data/output/schedules/qosqoe/qosqoe_{}{}.csv".format(sched,suffix), sep = ';')
            dfSchedGlobal[dfSchedGlobal['weighting_method']==sched].to_csv("data/output/schedules/qosqoe/schedule_global_objectives_{}{}.csv".format(sched,suffix), sep = ';')
            for day in days:
                # filter results for sched and day
                temp = temp_s[temp_s['day']==str(day)]
                # save csv with proper prefix
                temp.to_csv("data/output/schedules/qosqoe/qosqoe_{}_{}{}.csv".format(sched, day,suffix), sep=';')

        return

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
        self.rtFULPES1temp = []
        self.rtFULPES2temp = []

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
        self.rtFULPES1 += [self.rtFULPES1temp]
        self.rtFULPES2 += [self.rtFULPES2temp]

        if len(self.mbLPtemp) > 0:
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
            self.rtFOCSmpc_median += [statistics.median(self.rtFOCSmpctemp)]
        self.rtFOCS_median += [statistics.median(self.rtFOCStemp)]
        if len(self.rtFULPES1temp)>0:
            self.rtFULPES1_median += [statistics.median(self.rtFULPES1temp)]
            self.rtFULPES2_median += [statistics.median(self.rtFULPES2temp)]

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
    
    def write_cf(self, focs, avr, oa, greedy = None):
        if avr is not None:
            self.ovAVR += [avr.objective()]

        if oa is not None:
            self.ovOA += [oa.objective()]

        if focs is not None:
            self.ovFOCS += [focs.objective()]

        if greedy is not None:
            self.ovGreedy += [greedy.objective(greedy.f_cost)]
        
        if (focs is not None) and (avr is not None):
            self.ovRatioAVR += [avr.objective()/focs.objective()]

        if (focs is not None) and (oa is not None):
            self.ovRatioOA += [oa.objective()/focs.objective()]

        if (focs is not None) and (greedy is not None):
            self.ovRatioGreedy += [greedy.objective(greedy.f_cost)/focs.objective()]

        #load duration is later only a sorted list of this. Since x1 and x2 are 'independent' in the sense that we don't want them in the same plot.
        return
    
    def write_power(self,focs,avr,oa,greedy=None):
        if len(self.prefix) > 0:
            if greedy is None:
                header = ['i', 'focs_'+str(self.prefix[-1]), 'avr_'+str(self.prefix[-1]), 'oa_'+str(self.prefix[-1])]
            else:
                header = ['i', 'focs_'+str(self.prefix[-1]), 'avr_'+str(self.prefix[-1]), 'oa_'+str(self.prefix[-1]), 'greedy_'+str(self.prefix[-1])]
        else:
            if greedy is None:
                header = ['i', 'focs_', 'avr_', 'oa_']
            else:
                header = ['i', 'focs_', 'avr_', 'oa_', 'greedy_']

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
        if greedy is not None: # Do same for greedy
            temp = []
            # determine powers
            p_i = [((greedy.f_cost["i"+str(i)]['t']/(greedy.instance.len_i[i]))*greedy.instance.timeBase) for i in range(0,greedy.instance.m)]
            for i in range(0,greedy.instance.m):
                # extrapolate
                temp += [p_i[i] for t in range(0,int(greedy.instance.len_i[i]/greedy.instance.timeStep))]
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
            writer.writerow(['index', 'ovFOCS', 'ovAVR', 'ovOA', 'ovGreedy', 'ovRatioAVR', 'ovRatioOA', 'ovRatioGreedy'])
            #content
            writer.writerows(np.array([index, self.ovFOCS, self.ovAVR, self.ovOA, self.ovGreedy, self.ovRatioAVR, self.ovRatioOA, self.ovRatioGreedy]).T.tolist())
        
        #write competitive ratios sorted
        with open(path + 'cf_{}_sorted.csv'.format(self.prefix[-1].rsplit('_', 1)[0]), 'w', newline='' ) as f:
            writer = csv.writer(f, delimiter=';')
            #header
            writer.writerow(['index', 'ovFOCS', 'ovAVR', 'ovOA', 'ovGreedy', 'ovRatioAVR', 'ovRatioOA', 'ovRatioGreedy'])
            #content
            writer.writerows(np.array([index, sorted(self.ovFOCS), sorted(self.ovAVR), sorted(self.ovOA), sorted(self.ovGreedy), sorted(self.ovRatioAVR), sorted(self.ovRatioOA), sorted(self.ovRatioGreedy)]).T.tolist())
        
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

    def write_fulpes_instance_to_csv(self, pf = ''):
        path = 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data/output/fulpes/'
        with open(path + 'fulpes_instance_{}.csv'.format(pf), 'w', newline='' ) as f:
            writer = csv.writer(f, delimiter=';')
            #header
            writer.writerow(['rtFOCS', 'rtFULPES1', 'rtFULPES2'])
            #content
            writer.writerows(np.array([self.rtFOCStemp, self.rtFULPES1temp, self.rtFULPES2temp]).T.tolist())
            return
        
    def write_fulpes_to_csv(self,pf = '', pathOR = None, fc = None):
        path = 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data/output/fulpes/'
        if pathOR is not None:
            path = pathOR
        with open(path + 'fulpes_median_{}.csv'.format(pf), 'w', newline='' ) as f:
            writer = csv.writer(f, delimiter=';')
            if fc is not None:
                n = len(fc)
                #header
                writer.writerow(['n', 'rtFOCS', 'rtFULPES1', 'rtFULPES2'])
                #content
                writer.writerows(np.array([fc, self.rtFOCS_median[-n:], self.rtFULPES1_median[-n:], self.rtFULPES2_median[-n:], ]).T.tolist())
            else:
                #header
                writer.writerow(['rtFOCS', 'rtFULPES1', 'rtFULPES2'])
                #content
                writer.writerows(np.array([self.rtFOCS_median[-n:], self.rtFULPES1_median[-n:], self.rtFULPES2_median[-n:], ]).T.tolist())
        

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

    def write_schedule_objectives(self, schedule, prefix = ''):
        path = 'C:/Users/WinschermannL/OneDrive - University of Twente/Documenten/Gridshield/Criticalintervals/FOCS_code/data/output/schedules/qosqoe'

        # metrics per job
        with open(path + 'qosqoe_j_{}.csv'.format(prefix), 'w', newline='' ) as f:
            writer = csv.writer(f,delimiter=';')
            #header
            writer.writerow(['jobs_ens_abs_exact', 'jobs_ens_rel_exact'])
            #content
            writer.writerows(np.array([schedule.self.jobs_ens_abs_exact,schedule.self.jobs_ens_rel_exact]))

        # metrics global
        with open(path + 'qosqoe_{}.csv'.format(prefix), 'w', newline='' ) as f:
            writer = csv.writer(f,delimiter=';')
            #header
            writer.writerow(['ens_abs_exact', 'ens_rel_exact'])
            #content
            writer.writerows(np.array([schedule.self.ens_abs_exact,schedule.self.ens_rel_exact]))
        return
    