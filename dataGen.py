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

import pandas as pd
from config_path import Config
import datetime
import copy
import math

Config.initializePath()
Config.path("VoI_code")
from filter_data_process import FilteringData
Config.path("FOCS_data")

timeSteps = [1, 10, 60, 900, 1800, 3600]
data = FilteringData()

# Real Training data
instanceData = data.filter_data(startTimeFilter = True, 
                                    afterStartDate = datetime.datetime(2022, 9, 1, 0, 0), 
                                    beforeEndDate = datetime.datetime(2023, 8, 31, 23, 59), 
                                    energyFilter = True, 
                                    energyCutOff = 1,
                                    maxDwellTime = 24,
                                    minDwellTime = 10/60,
                                    overnightFilter = True,
                                    defaultCapacity = 100000,
                                    defaultPower = 11000,
                                    managersFilter = None, # True: analysis without manangers | False: analysis ONLY managers | None: Analysis considering all together! 
                                    listmanagersFilter = ['1000019032','1000019086','1000019033','1000019034','1000019030','1000019031','1000019036','1000019029','1000019035','1000019037','1000019040','1000019038','1000019039','1000011271','1000011271','1000011272','1000011272','1000011273','1000011273','1000011274','1000011255','1000011255','1000011275','1000011275','1000011276','1000011276','1000011277','1000011277','1000011278','1000011278','1000011317'],
                                    idFilter = None,#['Plug & charge', '1 (chip)', '2 (chip)', '3 (chip)', '4 (chip)', '5 (chip)', '6 (chip)', '7 (chip)', '8 (chip)', '9 (chip)', '10 (chip)', '11 (chip)', '12 (chip)', '13 (chip)', '14 (chip)', '15 (chip)', '16 (chip)', '17 (chip)', '18 (chip)', '19 (chip)', '20 (chip)', '21 (chip)', '22 (chip)', 'Anoniem'],
                                    discretization = timeSteps)

#determine cumulative probabilities for departure times given historical information till current time
instanceData['original_index'] = instanceData.index
instanceData.sort_values(by='start_datetime_utc')
instanceData.reset_index(drop=True,inplace=True)
# add column: count so far
cols = [[0]*len(instanceData)]
col_heads = ['count so far']
# add 96 columns: 
for i in range(1,97):
    temp_col = [0]*len(instanceData)
    # add column with prob that still avail throughout interval given historic data so far
    for ev in pd.unique(instanceData['card_id']):
        temp_data = copy.deepcopy(instanceData[instanceData['card_id']==ev])
        temp_n = len(temp_data)
        # start on top
        for occ in range(0,temp_n-1):
            temp_index = temp_data.index[0]
            # drop top row
            temp_data = temp_data.iloc[1:]
            # count occurrances so far
            cols[0][temp_index] = len(temp_data)
            # determine chance still avail based on past
            p = len(temp_data[temp_data['t1_900']>=i]) / max(1,len(temp_data))
            # save result
            temp_col[temp_index] = p 
    col_heads += ['cumprob_i'+str(i)]
    cols += [temp_col]
    print('[MESSAGE]: i = {} complete.'.format(i))
# determine mean departure so far
temp_mean = [0]*len(instanceData) 
for ev in pd.unique(instanceData['card_id']):
    temp_data = copy.deepcopy(instanceData[instanceData['card_id']==ev])
    temp_n = len(temp_data)
    # start on top
    for occ in range(0,temp_n-1):
        temp_index = temp_data.index[0]
        # drop top row
        temp_data = temp_data.iloc[1:]
        # determine average departure
        avg = sum(temp_data['end_datetime_seconds'])/len(temp_data)
        # save result
        temp_mean[temp_index] = avg
cols += [temp_mean]
col_heads += ['end_datetime_seconds_mean']

cols_df = pd.DataFrame(cols).T
cols_df.columns = col_heads
instanceData = pd.concat([instanceData, cols_df], axis=1)
# correct prediction for feasibility given energy and power
instanceData['assumed_power_W'] = [22000 if instanceData["average_power_W"].iloc[j]>instanceData["maxPower"].iloc[j] else instanceData["maxPower"].iloc[j] for j in range(0,len(instanceData))] # in W
instanceData['end_datetime_seconds_req'] = 3600 * instanceData['total_energy_Wh']/instanceData['assumed_power_W']

# discretize
for dis in timeSteps:
    instanceData['t1_' + str(dis) + '_mean'] = (instanceData['end_datetime_seconds_mean']/dis).apply(lambda x: math.ceil(x))
    instanceData['t1_' + str(dis) + '_req'] = instanceData['t0_'+ str(dis)] + (instanceData['end_datetime_seconds_req']/dis).apply(lambda x: math.ceil(x))
    instanceData['t1_' + str(dis) + '_mean_feas'] = instanceData['t1_' + str(dis) + '_mean']
    instanceData['t1_' + str(dis) + '_mean_feas'][instanceData['t1_' + str(dis) + '_mean_feas'] < instanceData['t1_' + str(dis) + '_req']] = instanceData['t1_' + str(dis) + '_req']

instanceData.to_excel("filteredData_cumprob.xlsx")

# data set for Marleen
# instanceData[['total_energy_Wh', 'average_power_W', 'maxPower'] + ['t0_'+str(timeStep) for timeStep in timeSteps] + ['t1_'+str(timeStep) for timeStep in timeSteps]][:200].to_excel("DEMSdata_FOCS_v1.xlsx",header=['total_energy_Wh', 'average_power_W', 'maxPower'] + ['t0_'+str(timeStep) for timeStep in timeSteps] + ['t1_'+str(timeStep) for timeStep in timeSteps])
# data set for DEMS group 4
# instanceData = data.anonymize()
# instanceData[['EV_id', 'start_datetime_utc', 'end_datetime_utc', 'total_energy_Wh', 'average_power_W', 'maxPower'] + ['t0_'+str(timeStep) for timeStep in timeSteps] + ['t1_'+str(timeStep) for timeStep in timeSteps]].to_excel("DEMSdata_FOCS_v2.xlsx",header=['EV_id', 'arrival', 'departure', 'total_energy_Wh', 'average_power_W', 'maxPower'] + ['t0_'+str(timeStep) for timeStep in timeSteps] + ['t1_'+str(timeStep) for timeStep in timeSteps])

# data feitjes
# print('n = ', len(instanceData))
# print('max nr sessions per day = ', data.max_sessions_per_day(afterStartDate = datetime.datetime(2022, 9, 1, 0, 0), beforeEndDate = datetime.datetime(2023, 8, 31, 23, 59)))
