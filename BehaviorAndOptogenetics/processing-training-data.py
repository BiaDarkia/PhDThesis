#!/usr/bin/env python

"""Script to extract information from lever press training and
from behavioral training on the benefit-benefit, cost-cost 
and cost-benefit decision-making condition from raw data"""

# Import statements
import os
import re
import warnings

from tkinter import filedialog
from tkinter import *

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
def count_lever_presses(text_data):
    if (int(float(re.sub(r'P:|\s', '', text_data[26]))) + int(float(re.sub(r'Q:|\s', '', text_data[27])))) == 0:
        return 0, 0, 20, 0, 0, 0, 0, 20, 0, 0
    
    else:
        lever_data = re.sub(r'\s+[0-9]+:\s+|\s+', ' ', " ".join(text_data[38:])).split("D:")
            
        lever_presses = np.array(lever_data[0].strip().split(" "), dtype='float32')
        reaction_time = np.array(lever_data[1].strip().split(" "), dtype='float32')
            
        lever_presses_tidy = lever_presses[np.logical_or(reaction_time == 0.0, np.logical_and(reaction_time >=10.0, reaction_time <= 20.0))]
        lever_presses_count_block1 = np.bincount(lever_presses_tidy.astype(int)[0:20], minlength = 3)
        lever_presses_count_block2 = np.bincount(lever_presses_tidy.astype(int)[20:40], minlength = 3)
    
        reaction_time_tidy = reaction_time[np.logical_or(reaction_time == 0.0, np.logical_and(reaction_time >=10.0, reaction_time <= 20.0))]
        reaction_time_tidy[np.logical_and(reaction_time_tidy >=10.0, reaction_time_tidy <= 20.0)] -= 10

        rt_left_block1 = reaction_time_tidy[0:20][lever_presses_tidy[0:20] == 1.0]
        rt_right_block1 = reaction_time_tidy[0:20][lever_presses_tidy[0:20] == 2.0]
        rt_left_block2 = reaction_time_tidy[20:40][lever_presses_tidy[20:40] == 1.0]
        rt_right_block2 = reaction_time_tidy[20:40][lever_presses_tidy[20:40] == 2.0]
    
        return lever_presses_count_block1[1], lever_presses_count_block1[2], lever_presses_count_block1[0], rt_left_block1, rt_right_block1, lever_presses_count_block2[1], lever_presses_count_block2[2], lever_presses_count_block2[0], rt_left_block2, rt_right_block2

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

# For the benefit-benefit and cost-cost decision-making paradigm,
# store the subject ID, the date, the decision-making paradigm.
# the block (first or last 20 trials), the group, high, low and omitted lever presses, 
# and reaction time for presses on the high and low lever. 
# The number of high and low lever presses either correspond to the 
# number of presses on the left or right lever depending on the group of the 
# subject, i.e. if the high cost-high benefit lever was the left or right lever 
# for that specific animal. The same holds true for the reaction time.
def append_all(subject, date, condition, block, group, left, right, omit, rt_left, rt_right):
    
    subjects.append(subject)
    dates.append(date)
    conditions.append(condition)
    blocks.append(block)

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

# For the cost-benefit decision-making paradigm,
# store the subject ID, the date, the decision-making paradigm.
# the group, high, low and omitted lever presses, 
# and reaction time for presses on the high and low lever. 
# The number of high and low lever presses either correspond to the 
# number of presses on the left or right lever depending on the group of the 
# subject, i.e. if the high cost-high benefit lever was the left or right lever 
# for that specific animal. The same holds true for the reaction time.
def append_all_cb(subject, date, condition, group, left, right, omit):
    
    subjects_cb.append(subject)
    dates_cb.append(date)
    conditions_cb.append(condition)

    if group == "LL":
        groups_cb.append(1)
    elif group == "RL":
        groups_cb.append(2)
    else:
        while True:
            group = input("Please input group_id for subject " + subject + "(1 for LL or 2 for RL): ")
            if (group != "1" and group != "2"):
                print("Please input 1 for LL and 2 for RL.")
                continue
            else:
                groups_cb.append(int(group))
                break
 
    if groups_cb[-1] == 1:
        high_cb.append(left)
        low_cb.append(right)
    if groups_cb[-1] == 2:
        high_cb.append(right)
        low_cb.append(left)
        
    omitted_cb.append(omit)

# For l,ever press training, store the subject ID, the date, 
# the condition (left lever or right lever as high benefit lever),
# the group and the number of lever presses for each animal.
def append_all_lp(subject, date, condition, group, press):

    subjects_lp.append(subject)
    dates_lp.append(date)
    conditions_lp.append(condition)

    if group == "LL":
        groups_lp.append(1)
    elif group == "RL":
        groups_lp.append(2)
    else:
        while True:
            group = input("Please input group_id for subject " + subject + "(1 for LL or 2 for RL): ")
            if (group != "1" and group != "2"):
                print("Please input 1 for LL and 2 for RL.")
                continue
            else:
                groups_lp.append(int(group))
                break
        
    presses.append(press)

# Initiate empty lists and an empty dictionary to append/use later on 
subjects = []
dates = []
conditions = []
blocks = []
groups = []
high = []
low = []
omitted = []
rt_high = []
rt_low = []

subjects_cb = []
dates_cb = []
conditions_cb = []
groups_cb = []
high_cb = []
low_cb = []
omitted_cb = []

subjects_lp = []
dates_lp = []
conditions_lp = []
groups_lp = []
presses = []

cache = dict()
    
# Open a dialog to select the folder that the data are stored in   
root = Tk()
root.withdraw()
path =  filedialog.askdirectory(title = "Select directory...")
root.destroy()

# For each file in the selected directory open the file, extract the subject ID,
# experimental condition as specified in the file (condition), the block 
# (first or last 20 trials), if the high cost-high benefit lever was 
# the left or right lever (group), the number of high, low and omitted 
# lever presses, and the reaction times for high and low lever
# presses, and store them in lists. For the cost-benefit decision-making
# paradigm and for lever press training only a subset of these variables
# is extracted and stored, and for lever pressing training the overall number
# of lever presses is extracted and stored.   
for file in os.listdir(path):
    with open(os.path.join(path, file), 'r') as openfile:
        # Read all lines of the open file
        raw_data = openfile.readlines()
        # Extract the date of the experiment
        date = datetime.datetime.strptime(raw_data[4][12:20], "%m/%d/%y").date()
        # Extract the subject ID
        subject = re.sub(r'Subject:|\s', '', raw_data[6])
        
        # Extract the condition, i.e. BB-DM for benefit-benefit decision-making,
        # CC-DM for cost-cost decision-making, CB-DM for cost-benefit 
        # decision-making, LP-LL for lever pressing training on the
        # left lever and LP-RL for lever pressing training on the right lever. 
        # Conditions are coded by numbers. If input of the condition was not
        # correct when the data was recorded, i.e. it does not match any of 
        # the short notations above, a value error is raised
        condition = re.sub(r'Experiment:|-[0-9]+|\s', '', raw_data[7])
      
        if condition == "BB-DM":
            condition = 1
        elif condition == "CC-DM":
            condition = 2
        elif condition == "CB-DM":
            condition = 3
        elif condition == "LP-LL":
            condition = 11
        elif condition == "LP-RL":
            condition = 12
        else:
            raise ValueError("Please make sure that Experiment in " + path + "/" + file + " specifies a valid condition, i.e. BB-DM, BB-DM-OP, CC-DM, CC-DM-OP, CB-DM or CB-DM-OP.") 
        
        # Extract the group, i.e. if the high cost-high benefit lever was the 
        # left or right lever for this specific animal as provided in the file
        group = raw_data[8][7:9]
        
        # If the condition if either BB-DM or CC-DM extract the lever presses
        # on the left and right lever on the first 20 and last 20 trials using
        # the function count_lever_presses and use append_all to add the
        # extracted data to lists of previously extracted data
        if condition in [1, 2]:
            l1, r1, o1, rl1, rr1, l2, r2, o2, rl2, rr2  = count_lever_presses(raw_data)
            
            block = 0
            append_all(subject, date, condition, block, group, l1, r1, o1, rl1, rr1)
            
            block = 1
            append_all(subject, date, condition, block, group, l2, r2, o2, rl2, rr2)
        
        # If the condition is CB extract left, right and omitted lever presses
        # directly from the data and append the information using append all.
        # Information on lever presses is not actually used in later data
        # analysis and, hence, I did not account for occurrences when animals
        # manage to press the lever a second time while it retracts
        elif condition == 3:
            left = int(float(re.sub(r'P:|\s', '', raw_data[26])))
            right = int(float(re.sub(r'Q:|\s', '', raw_data[27])))        
            omit = int(float(re.sub(r'R:|\s', '', raw_data[28])))
            
            append_all_cb(subject, date, condition, group, left, right, omit)
        
        # If the condition is LP-LL or LP-RL extract and store the number of
        # lever presses
        elif condition in [11, 12]:
            press = int(float(re.sub(r'P:|\s', '', raw_data[27])))
            append_all_lp(subject, date, condition, group, press)

rt_high = np.nan_to_num(rt_high)
rt_low = np.nan_to_num(rt_low)

# Store data from the BB-DM and CC-DM condition in a pandas data frame, 
# calculate the total number of trials and the 'performance', 
# i.e. the number of high benefit/high cost/high cost-high benefit 
# lever presses as compared to the overall number of trials 
# that were not omitted as well as the difference in reaction times on
# high benefit/high cost/high benefit-high cost as compared to
# low benefit/low cost/low benefit-low cost choices
data = pd.DataFrame({"subject": subjects, "date": dates, "condition": conditions, "block": blocks, "group": groups, "high": high, "low": low, "omitted": omitted, "rt_high": rt_high, "rt_low": rt_low})
data['performance'] = (data['high'] / (data['low'] + data['high'])) * 100
data['performance'].fillna(0, inplace=True)
data['rt_diff'] = data['rt_high'] - data['rt_low']

# Based on the date calculate the day of behavioral training and store it in
# a variable day that codes the last day of behavioral testing as "-1", the
# second to last day as "-2" and so on
data['day'] = ""
data_days = []

for ele in list(set(subjects)):
    for cond in [1, 2]:
        data_day = data[(data['subject']==ele) & (data['condition']==cond)]
        data_day = data_day.sort_values('date', ascending=False)
        days = np.negative(np.arange(1, (len(data_day.index)/2)+1, 0.5).astype(int))
        data_day['day'] = days
        data_days.append(data_day)

# Add the day column, drop the date column and store data in a csv file
data = pd.concat(data_days)
data = data.drop(columns=['date'])
data.to_csv("./tidy_data_training.csv", index = False)

# Store data from lever pressing training, i.e. from the LP-LL and LP-RL
# condition, in a pandas data frame and store the data in a csv file
data_lp = pd.DataFrame({"subject": subjects_lp, "date": dates_lp, "condition": conditions_lp, "group": groups_lp, "presses": presses})
data_lp.to_csv("./tidy_data_training_lp.csv", index = False)

# Store data from the CB-DM condition in a pandas data frame,
# calculate the 'performance', i.e. the number of high benefit/high cost/
# high cost-high benefit lever presses as compared to the overall number 
# of trials that were not omitted, and store the data in a csv file
data_cb = pd.DataFrame({"subject": subjects_cb, "date": dates_cb, "condition": conditions_cb, "group": groups_cb, "high": high_cb, "low": low_cb, "omitted": omitted_cb})
data_cb['performance'] = (data_cb['high'] / (data_cb['low'] + data_cb['high'])) * 100
data_cb['performance'].fillna(0, inplace=True)
data_cb.to_csv("./tidy_data_training_cb.csv", index = False)
