#!/usr/bin/env python

"""Script to extract information for benefit-benefit, cost-cost 
and cost-benefit decision-making condition from raw data"""

# Import statements
import os
import re
import warnings

from tkinter import filedialog
from tkinter import *

from scipy import arange

import numpy as np
import pandas as pd
import datetime

# Sometimes animals manage to press the lever a second time when it retracts.
# In this case a file may contain data on more than 20 trials. If this is 
# the case, this function determines and returns the number of left, right 
# and omitted lever presses as well as reaction times for left and right
# lever presses. First, any second time lever presses are excluded.
# Because reaction time is recorded from the beginning of the trial rather
# than from the time of lever retraction this can easily be done by
# excluding any lever presses with a reaction time shorter than 10 secs.
# However, if a trial is omitted the reaction time is zero and, hence,
# these trials need to be kept. After cleaning up the data the first 20 trials
# are selected and the number of left, right and omitted lever presses as well
# as the corresponding reaction times are determined and returned
def count_lever_presses(text_data, length=0):        
    lever_data = re.sub(r'\s+[0-9]+:\s+|\s+', ' ', " ".join(text_data[38:])).split("D:")
            
    lever_presses = np.array(lever_data[0].strip().split(" "), dtype='float32')
    reaction_time = np.array(lever_data[1].strip().split(" "), dtype='float32')
            
    lever_presses_tidy = lever_presses[np.logical_or(reaction_time == 0.0, np.logical_and(reaction_time >=10.0, reaction_time <= 20.0))]
    reaction_time_tidy = reaction_time[np.logical_or(reaction_time == 0.0, np.logical_and(reaction_time >=10.0, reaction_time <= 20.0))]

    reaction_time_tidy[np.logical_and(reaction_time_tidy >=10.0, reaction_time_tidy <= 20.0)] -= 10

    if (length == 0):
        lever_presses_count = np.bincount(lever_presses_tidy.astype(int), minlength = 3)
        rt_left = reaction_time_tidy[lever_presses_tidy == 1.0]
        rt_right = reaction_time_tidy[lever_presses_tidy == 2.0]
        
    else:
        lever_presses_count = np.bincount(lever_presses_tidy.astype(int)[0:length], minlength = 3)
        rt_left = reaction_time_tidy[0:length][lever_presses_tidy[0:length] == 1.0]
        rt_right = reaction_time_tidy[0:length][lever_presses_tidy[0:length] == 2.0]
    
    return lever_presses_count[1], lever_presses_count[2], lever_presses_count[0], rt_left, rt_right

# For the case that a file contains excatly 20 trials, this function will
# extracts the reaction time for lever presses on the left or right lever
def extract_reaction_time(text_data):
    lever_data = re.sub(r'\s+[0-9]+:\s+|\s+', ' ', " ".join(text_data[38:])).split("D:")
            
    lever_presses = np.array(lever_data[0].strip().split(" "), dtype='float32')
    reaction_time = np.array(lever_data[1].strip().split(" "), dtype='float32')

    reaction_time[np.logical_and(reaction_time >=10.0, reaction_time <= 20.0)] -= 10
       
    rt_left = reaction_time[lever_presses == 1.0]
    rt_right = reaction_time[lever_presses == 2.0]
    
    return rt_left, rt_right 

# Store the subject ID, the date, the decision-making paradigm.
# the optogenetic condition, the group, high, low and omitted lever presses, 
# and reaction time for presses on the high and low lever. 
# The number of high and low lever presses either correspond to the 
# number of presses on the left or right lever depending on the group of the 
# subject, i.e. if the high cost-high benefit lever was the left or right lever 
# for that specific animal. The same holds true for the reaction time.
def append_all(subject, date, condition, opto, group, left, right, omit, rt_left, rt_right):
    
    subjects.append(subject)
    dates.append(date)
    conditions.append(condition)
    optogenetics.append(opto)

    if group == "LL":
        groups.append(1)
    elif group == "RL":
        groups.append(2)
    else:
        while True:
            group = input("Please input group_id for subject " + subject + "(1 for LL or 2 for RL): ")
            if (group != "1" and group != "2"):
                print("Please input 1 for LL and 2 for RL.")
                continue
            else:
                groups.append(int(group))
                break
        
    if groups[-1] == 1:
        high.append(left)
        low.append(right)
        rt_high.append(np.mean(rt_left))
        rt_low.append(np.mean(rt_right))
    if groups[-1] == 2:
        high.append(right)
        low.append(left)
        rt_high.append(np.mean(rt_right))
        rt_low.append(np.mean(rt_left))
        
    omitted.append(omit)

# Initiate empty lists and an empty dictionary to append/use later on
subjects = []
dates = []
conditions = []
optogenetics = []
groups = []
high = []
low = []
omitted = []
rt_high = []
rt_low = []
cache = dict()   

# Open a dialog to select the folder that the data are stored in
root = Tk()
root.withdraw()
path =  filedialog.askdirectory(title = "Select directory...")
root.destroy()         

# For each file in the selected directory open the file, extract the subject ID,
# experimental condition as specified in the file (condition), if light to
# induce optogenetic inhibition was delivered (opto), if the high cost-high 
# benefit lever was the left or right lever (group), the number of high, low 
# and omitted lever presses, and the reaction times for high and low lever
# presses, and store them in lists       
for file in os.listdir(path):
    with open(os.path.join(path, file), 'r') as openfile:
        # Read all lines of the open file
        raw_data = openfile.readlines()
        # Extract the date of the experiment
        date = datetime.datetime.strptime(raw_data[4][12:20], "%m/%d/%y").date()
        # Extract the subject ID
        subject = re.sub(r'Subject:|\s', '', raw_data[6])
        
        # Extract the condition, i.e. BB-DM for benefit-benefit decision-making,
        # BB-DM-OP for benefit-benefit decision-making with delivery of the 
        # light to induce optogenetic inhibition, CC-DM for cost-cost
        # decision-making, CC-DM-OP for cost-cost decision-making with 
        # delivery of the light to induce optogenetic inhibition, CB-DM for 
        # cost-benefit decision-making, and CB-DM-OP for cost-benefit
        # decision-making with delivery of the light to induce 
        # optogenetic inhibition. The condition is then coded as "1" for
        # benefit-benefit decision-making, "2" for cost-cost decision-making
        # and "3" for cost-benefit decision-making, and as "0" if no light
        # is delivered of "1" if light to induce optogenetic inhibition
        # was delivered. If input of the condition was not correct when 
        # the data was recorded, i.e. it does not match any of the short 
        # notations above, a value error is raised
        condition = re.sub(r'Experiment:|-[0-9]+|\s', '', raw_data[7])
        
        if condition == "BB-DM":
            condition = 1
            opto = 0
        elif condition == "BB-DM-OP":
            condition = 1
            opto = 1
        elif condition == "CC-DM":
            condition = 2
            opto = 0
        elif condition == "CC-DM-OP":
            condition = 2
            opto = 1
        elif condition == "CB-DM":
            condition = 3
            opto = 0
        elif condition == "CB-DM-OP":
            condition = 3
            opto = 1
        else:
            raise ValueError("Please make sure that Experiment in " + path + "/" + file + " specifies a valid condition, i.e. BB-DM, BB-DM-OP, CC-DM, CC-DM-OP, CB-DM or CB-DM-OP.") 
        
        # Extract the group, i.e. if the high cost-high benefit lever was the 
        # left or right lever for this specific animal as provided in the file
        group = raw_data[8][7:9]
        
        # Extract the number of presses on the left and right lever,
        # and the number of omitted lever presses
        left = int(float(re.sub(r'P:|\s', '', raw_data[26])))
        right = int(float(re.sub(r'Q:|\s', '', raw_data[27])))        
        omit = int(float(re.sub(r'R:|\s', '', raw_data[28])))
        
        # Calculate the total number of lever presses to check 
        # if the file contains 20 trials
        total = left + right + omit
        
        # The file may contain less than 20 trials, if the transponder attached
        # to the implanted LED fiber to remotely turn on and off the LED was
        # disconnected during the experiment, the experiment had to be stopped
        # and the remaining trials were recorded after the transponder was
        # reconnected to the LED. In this case, the number of presses on
        # the left and right lever, the number of omitted lever presses, and 
        # reaction times for left and right lever presses is extracted
        # and stored in a dictionary using a unique key for that specific
        # animal, date and condition. Whenever another file with data
        # recorded for the same animal on the same date and condition is
        # encountered the current dictionary entry is appended till a data of
        # a total of 20 trials is available. Once data for 20 or more trials
        # is available, values for subject ID, date, condition, opto, group,
        # number of high, low and omitted lever presses and reaction times
        # for high and low lever presses will be stored in lists
        if (total < 20):
            cache_key = subject + "/" + str(date) + "/" + str(condition) + "/" + str(opto) + "/" + group
            if cache_key in cache:
                values = cache[cache_key]
                if (total + sum(values[0:3]) < 20):
                    left, right, omit, rt_left, rt_right = count_lever_presses(raw_data, 0)
                    cache[cache_key] = [values[0] + left, values[1] + right, values[2] + omit, np.concatenate([values[3], rt_left]), np.concatenate([values[4], rt_right])]
                if (total + sum(values[0:3]) == 20):
                    rt_left, rt_right = extract_reaction_time(raw_data)                    
                    append_all(subject, date, condition, opto, group, values[0] + left, values[1] + right, values[2] + omit, np.concatenate([values[3], rt_left]), np.concatenate([values[4], rt_right]))
                    del cache[cache_key]
                if (total + sum(values[0:3]) > 20):
                    left, right, omit, rt_left, rt_right = count_lever_presses(raw_data, 20 - sum(values[0:3]))   
                    append_all(subject, date, condition, opto, group, values[0] + left, values[1] + right, values[2] + omit, np.concatenate([values[3], rt_left]), np.concatenate([values[4], rt_right]))
                    del cache[cache_key]
            else:
                rt_left, rt_right = extract_reaction_time(raw_data)
                cache[cache_key] = [left, right, omit, rt_left, rt_right]
        
        # If the file contains exactly 20 trials, extract the values for 
        # subject ID, date, condition, opto, group, number of high, low and 
        # omitted lever presses and reaction times for high and low 
        # lever presses and store them using the append_all function    
        if (total == 20):
            cache_key = subject + "/" + str(date) + "/" + str(condition) + "/" + str(opto) + "/" + group
            if cache_key in cache:
                values = cache[cache_key]
                left, right, omit, rt_left, rt_right = count_lever_presses(raw_data, 20 - sum(values[0:3]))          
                append_all(subject, date, condition, opto, group, values[0] + left, values[1] + right, values[2] + omit, np.concatenate([values[3], rt_left]), np.concatenate([values[4], rt_right]))
                del cache[cache_key]
            else:
                rt_left, rt_right = extract_reaction_time(raw_data)
                append_all(subject, date, condition, opto, group, left, right, omit, rt_left, rt_right)
        
        # Sometimes animals manage to press the lever a second time when
        # it retracts. In this case a file may contain data on more than
        # 20 trials. If this is the case, use the count_lever_presses function
        # to extract values for subject ID, date, condition, opto, group, number of high,
        # low and omitted lever presses and reaction times for high and low 
        # lever presses and store them using the append_all function 
        if (total > 20):
            left, right, omit, rt_left, rt_right = count_lever_presses(raw_data, 20)
            append_all(subject, date, condition, opto, group, left, right, omit, rt_left, rt_right)

# If there is any data remaining in the dictionary add it to the lists.
# This may be the case if there are any missing trials            
if cache:
    for key, value in cache.items():
        sdcog = key.split("/")
        append_all(sdcog[0], datetime.datetime.strptime(sdcog[1], "%Y-%m-%d").date(), int(sdcog[2]), int(sdcog[3]), sdcog[4], value[0], value[1], value[2], value[3], value[4])

rt_high = np.nan_to_num(rt_high)
rt_low = np.nan_to_num(rt_low)

# Store data in pandas data frame, aggregate data from the same condition,
# calculate the total number of trials and the 'performance',
# i.e. the number of high benefit/high cost/high cost-high benefit 
# lever presses as compared to the overall number of trials that
# were not omitted
data = pd.DataFrame({"subject": subjects, "condition": conditions, "opto": optogenetics, "group": groups, "high": high, "low": low, "omitted": omitted, "rt_high": rt_high, "rt_low": rt_low})
data = data.groupby(['subject', 'condition', 'opto', 'group'], as_index=False).agg({'high': np.sum, 'low': np.sum, 'omitted': np.sum, 'rt_high': np.mean, 'rt_low': np.mean})
data['total'] = data['high'] + data['low'] + data['omitted']
data['performance'] = (data['high'] / (data['low'] + data['high'])) * 100
data['performance'].fillna(0, inplace=True)
data['rt_diff'] = data['rt_high'] - data['rt_low']

# Provide a warning if there are any missing trials
if not all(x == 60 for x in data['total']):
    warnings.warn("Missing trials detected! Please check the data!")

# Store the data as ./tidy_data.csv
data.to_csv("./tidy_data.csv", index = False)

# Store data in a pandas data frame and sort by subject ID and date.
# For each subject each date should be present twice, once for
# opto = "0" and once for opto = "1". Conditions should be sorted based
# on the date, since animals were always first tested on the benefit-benefit,
# then on the cost-cost and last on the cost-benefit decision-making paradigm
# For each paradigm animals performed 20 trials on the opto = "0" and
# 20 trials on the opto = "1" condition on three consecutive days.
# Calculate the 'performance', i.e. the number of high benefit/high cost/
# high cost-high benefit lever presses as compared to the overall number of 
# trials that were not omitted. To discriminate data from the three
# consecutive days of testing for each condition, introduce a variable 
# day that codes data from day 1 as "1", from day 2 as "2" and from day 3 
# as "3". Drop the date column and store the data as ./tidy_data_day.csv
data_days = pd.DataFrame({"subject": subjects, "condition": conditions, "opto": optogenetics, "group": groups, "high": high, "low": low, "omitted": omitted, "rt_high": rt_high, "rt_low": rt_low, "date": dates})
data_days = data_days.sort_values(by=['subject','date'])
data_days['performance'] = (data_days['high'] / (data_days['low'] + data_days['high'])) * 100
data_days['performance'].fillna(0, inplace=True)
data_days['rt_diff'] = data_days['rt_high'] - data_days['rt_low']
data_days['day'] = int(len(data_days)/6) * [1,1,2,2,3,3]
data_days = data_days.drop(columns=['date'])

data_days.to_csv("./tidy_data_day.csv", index = False)                          