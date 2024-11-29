#    Data pre-processing. 
#    Flow-based Offline Charging Scheduler (FOCS)
#    For scheduling of large groups of EVs
#    Part of SmoothEMS met GridShield project - code developed by 
#    Leander van der Bijl, University of Twente, l.c.vanderbijl@utwente.nl
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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import math
import networkx as nx
import random
import datetime


#Functions

#Identify busses, note that indices muse be reset here
#algorithm: Look for when EndSoc != 100 but the next StartSoc == 100
def IdentifyBusses(df):
    # A dataframe of all busses that must be charged, Need is kWh to be ingested, Capacity is battery capacity, StartTime is the time the bus is plugged in (We assume universal plugin time for now)
    ToCharge = pd.DataFrame(columns=['Need','Capacity', 'StartTime'])
    for index, row in df.iterrows():
        endsoc = df.iloc[index]['EndSoc']
        if index != len(df.index) - 1:
            #Here I did 98 because I found a value of 98.3 before and I am not sure whether this bus is used again later as a full bus, should be insignificant to include.
            if endsoc < 98 and df.iloc[index + 1]['StartSoc'] == 100:
                capacity = df.iloc[index]['Capacity']
                ToCharge.loc[len(ToCharge.index)] = [capacity - endsoc, capacity, df.iloc[index]['Endtime']]
        else:
                capacity = df.iloc[index]['Capacity']
                ToCharge.loc[len(ToCharge.index)] = [capacity - endsoc, capacity, df.iloc[index]['Endtime']]
    return ToCharge



#This will give all start times of busses along with capacities required for the next day
def FindStartTimes(df):
    StartTimes = pd.DataFrame(columns=['Capacity', 'StartTime'])
    for index, row in df.iterrows():
        if index != 0:
            if df.iloc[index - 1]['EndSoc'] < 98 and df.iloc[index]['StartSoc'] == 100:
                capacity = df.iloc[index]['Capacity']
                StartTimes.loc[len(StartTimes.index)] = [capacity, df.iloc[index]['StartTime']]
        else:
            capacity = df.iloc[index]['Capacity']
            StartTimes.loc[len(StartTimes.index)] = [capacity, df.iloc[index]['StartTime']]
    return StartTimes



def PlotSchedule(MoThu):
    MoThu['StartTime'] = pd.to_datetime(MoThu['StartTime'], format='%H:%M:%S')
    MoThu['EndTime'] = pd.to_datetime(MoThu['EndTime'], format='%H:%M:%S')

    # Create a new figure and axis
    fig, ax = plt.subplots()

    # For each row in the DataFrame, plot a line from the start time to the end time
    for i, row in MoThu.iterrows():
        ax.plot([row['StartTime'], row['EndTime']], [i, i], color=plt.cm.Reds(row['Capacity']/MoThu['Capacity'].max()))

    # Set the x-axis to display time only
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    # Set the y-axis to display the indices of the DataFrame
    ax.set_yticks(range(len(MoThu)))
    ax.set_yticklabels(MoThu.index)

    # Set the labels and title
    ax.set_xlabel('Time')
    ax.set_ylabel('Index')
    ax.set_title('Start and End Times')

    # Show the plot
    plt.show()

#Create schedule for certain day. The groups are: Mo-Thu, Thu-Fri, Fri-Sat, Sat-Sun, Sun-Mon
def CreateSchedule(day):
    # Read the excel file
    activities = pd.read_excel("data/input/BusData.xlsx")




    #Clean the data
    #Drop activities without SOC
    activities = activities.dropna(subset=['StartSoc'])

    #Convert capacities to float
    activities['EndSoc'] = activities['EndSoc'].str.rstrip('%').astype('float')
    activities['StartSoc'] = activities['StartSoc'].str.rstrip('%').astype('float')

    #Convert times to datetime
    activities['StartTime'] = pd.to_datetime(activities['StartTime'], format='%H:%M:%S')
    activities['Endtime'] = pd.to_datetime(activities['Endtime'], format='%H:%M:%S')

    #Add timedifference
    activities['TimeDifference'] = (activities['Endtime'] - activities['StartTime']).dt.total_seconds() / 60.0

    #Add a column that represents capacity
    activities['Capacity'] = ((100/(activities['EndSoc'] - activities['StartSoc'])) * activities['IngestedKwh']).round()

    #Convert times to time
    activities['StartTime'] = pd.to_datetime(activities['StartTime']).dt.time
    activities['Endtime'] = pd.to_datetime(activities['Endtime']).dt.time

    #Reset indices
    activities = activities.reset_index(drop=True)

    #Add capacity of end of line busses to the activities
    capacity = 0
    for index, row in activities.iterrows():
        if index != len(activities.index) - 1:
            if activities.iloc[index]['Capacity'] > 0:
                capacity = activities.iloc[index]['Capacity'].round()
            if activities.iloc[index]['EndSoc'] < 98 and activities.iloc[index + 1]['StartSoc'] == 100:
                activities.at[index, 'Capacity'] = capacity
                capacity = 0
        else:
            if activities.iloc[index]['Capacity'] > 0:
                capacity = activities.iloc[index]['Capacity'].round()
            activities.at[index, 'Capacity'] = capacity
            capacity = 0

    #Add this capacity to all of the busses in the data
    indexlist = []
    for index, row in activities.iterrows():
        if pd.isna(activities.iloc[index]['Capacity']):
            indexlist.append(index)
        else:
            for indx in indexlist:
                activities.at[indx, 'Capacity'] = activities.iloc[index]['Capacity']
            indexlist = []


    #Start by splitting data into days of week
    if day == 'Mon-Thu':
        Incoming = activities[activities['DaysOfWeek'] == '1234*']
        Leaves = activities[activities['DaysOfWeek'] == '1234*']
    elif day == 'Thu-Fri':
        Incoming = activities[activities['DaysOfWeek'] == '1234*']
        Leaves = activities[activities['DaysOfWeek'] == '****5']
    elif day == 'Fri-Sat':
        Incoming = activities[activities['DaysOfWeek'] == '****5']
        Leaves = activities[activities['DaysOfWeek'] == 6]
    elif day == 'Sat-Sun':
        Incoming = activities[activities['DaysOfWeek'] == 6]
        Leaves = activities[activities['DaysOfWeek'] == 7]
    elif day == 'Sun-Mon':
        Incoming = activities[activities['DaysOfWeek'] == 7]
        Leaves = activities[activities['DaysOfWeek'] == '1234*']


    #Reset indices
    Incoming = Incoming.reset_index(drop=True)
    Leaves = Leaves.reset_index(drop=True)

    #Identify busses
    ToCharge = IdentifyBusses(Incoming)

    #Find start times
    StartTimes = FindStartTimes(Leaves)

    #Convert the ToCharge and StartTimes time to 900 second intervals
    ToCharge['t0_900'] = pd.to_datetime(ToCharge['StartTime'], format='%H:%M:%S').dt.hour * 4 + pd.to_datetime(ToCharge['StartTime'], format='%H:%M:%S').dt.minute * (1/15) + pd.to_datetime(ToCharge['StartTime'], format='%H:%M:%S').dt.second*(1/(60*15)) +1
    ToCharge['t0_900'] = ToCharge['t0_900'].astype(int)

    StartTimes['t1_900'] = pd.to_datetime(StartTimes['StartTime'], format='%H:%M:%S').dt.hour * 4 + pd.to_datetime(StartTimes['StartTime'], format='%H:%M:%S').dt.minute * (1/15) + pd.to_datetime(StartTimes['StartTime'], format='%H:%M:%S').dt.second*(1/(60*15)) +1
    StartTimes['t1_900'] = StartTimes['t1_900'].astype(int)

    for index, row in StartTimes.iterrows():
        #if ToCharge.iloc[index]['t0_900'] > StartTimes.iloc[index]['t1_900']:
        StartTimes.at[index, 't1_900']= StartTimes.iloc[index]['t1_900'] + 96

    #This is the optimal matching
    G = nx.Graph()
    for index, row in ToCharge.iterrows():
        for index2, row in StartTimes.iterrows():
            if (StartTimes.iloc[index2]['t1_900'] - ToCharge.iloc[index]['t0_900'])*30/4 - ToCharge.iloc[index]['Need'] >= 0:
                if StartTimes.iloc[index2]['Capacity'] == ToCharge.iloc[index]['Capacity']:
                    G.add_edge(index, index2+100)


    matching = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)

    addlist = []
    for index, row in ToCharge.iterrows():
        found = 0
        for edge in matching:
            if edge[0] < 100:
                edge0 = edge[0]
            else:
                edge0 = edge[1]
            if index == edge0:
                found = 1
        if found == 0:
            addlist.append(index)



    #Make the new dataframe
    MoThu = pd.DataFrame(columns=['Capacity', 't0_900', 't1_900', 'average_power_W'])
    for edge in matching:
        index0 = edge[0]
        index1 = edge[1]
        if index0 < 100:
            MoThu.loc[len(MoThu.index)] = [ToCharge.iloc[index0]['Need']*1000, ToCharge.iloc[index0]['t0_900'], StartTimes.iloc[index1-100]['t1_900'], 30000]
        else:
            MoThu.loc[len(MoThu.index)] = [ToCharge.iloc[index1]['Need']*1000, ToCharge.iloc[index1]['t0_900'], StartTimes.iloc[index0-100]['t1_900'], 30000]

    for index in addlist:
        MoThu.loc[len(MoThu.index)] = [ToCharge.iloc[index]['Need']*1000, ToCharge.iloc[index]['t0_900'], 193, 30000]


    # instanceData['maxpossible'] = (instanceData['t1_900'] - instanceData['t0_900'])*instanceData["maxPower"]/4
    # instanceData['plsbepositive'] = instanceData['maxpossible'] - instanceData['total_energy_Wh']



    MoThu['maxPower'] = 30000
    MoThu['average_power_W'] = 30000
    MoThu = MoThu.rename(columns={'StartTime': 'arrival', 'EndTime': 'departure', 'Capacity': 'total_energy_Wh'})


    return(MoThu)


def generate_baseload_random(num_timeslots, minval, maxval):
    #15 minutes time intervals
    baseload = []

    for t in range(0,num_timeslots):
        baseload.append(random.randint(minval, maxval))

    return baseload

def add_baseload(data, baseload):
    #We add the baseload by adding jobs. These jobs have to be done in a single interval and therefore we can simply give it a connection to only the time interval in which the load is present.
    #Find highest index job
    #Time starts now at 1 and ends at 192 (2 day spannen over 15 minute intervals)
    chargcap = max(baseload)*5 #assuming time intervals <= 1 hour

    timeslot = 1
    for job in baseload:
        data = data._append({'total_energy_Wh': job,'t0_900': timeslot, 't1_900': timeslot + 1, 'average_power_W': chargcap, 'maxPower': chargcap}, ignore_index=True)
        timeslot +=1
    
    return data

def generateDataWeek():
    #Add the dataframes of several days together
    data0 = CreateSchedule('Sun-Mon')  

    data1 = CreateSchedule('Mon-Thu')
    for index, row in data1.iterrows():
        data1.at[index, 't0_900'] += 96
        data1.at[index, 't1_900'] += 96

    data2 = CreateSchedule('Mon-Thu')
    for index, row in data2.iterrows():
        data2.at[index, 't0_900'] += 96*2
        data2.at[index, 't1_900'] += 96*2

    data3 = CreateSchedule('Mon-Thu')
    for index, row in data3.iterrows():
        data3.at[index, 't0_900'] += 96*3
        data3.at[index, 't1_900'] += 96*3

    data4 = CreateSchedule('Thu-Fri')
    for index, row in data4.iterrows():
        data4.at[index, 't0_900'] += 96*4
        data4.at[index, 't1_900'] += 96*4

    data5 = CreateSchedule('Fri-Sat')
    for index, row in data5.iterrows():
        data5.at[index, 't0_900'] += 96*5
        data5.at[index, 't1_900'] += 96*5

    data6 = CreateSchedule('Sat-Sun')
    for index, row in data6.iterrows():
        data6.at[index, 't0_900'] += 96*6
        data6.at[index, 't1_900'] += 96*6

    data7 = CreateSchedule('Sun-Mon')
    for index, row in data7.iterrows():
        data7.at[index, 't0_900'] += 96*7
        data7.at[index, 't1_900'] += 96*7

    data8 = CreateSchedule('Mon-Thu')
    for index, row in data8.iterrows():
        data8.at[index, 't0_900'] += 96*8
        data8.at[index, 't1_900'] += 96*8
    

    #total_time_horizon = 96*10 + 1


    #combine data
    data = data0
    data = data._append(data1)
    data = data._append(data2)
    data = data._append(data3)
    data = data._append(data4)
    data = data._append(data5)
    data = data._append(data6)
    data = data._append(data7)
    data = data._append(data8)


    return data


def import_baseload_data():
    #Import the baseload data
    baseload = pd.read_excel("data/input/baseloadData.xlsx")
    baseload = baseload.loc[14779:15737]
    baseload = baseload.reset_index(drop=True)
    baseloadlist = []
    for index, row in baseload.iterrows():
        baseloadlist.append(baseload.iloc[index]['Verbruik totaal']*1000)
    return baseloadlist
    
#There are 34 30kW chargers on site