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

# data feitjes
# print('n = ', len(instanceData))
# print('max nr sessions per day = ', data.max_sessions_per_day(afterStartDate = datetime.datetime(2022, 9, 1, 0, 0), beforeEndDate = datetime.datetime(2023, 8, 31, 23, 59)))
