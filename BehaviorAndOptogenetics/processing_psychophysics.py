#!/usr/bin/env python

"""Processing of data from training on the cost-benefit decision-making
paradigm that were acquired to plot a psychometric function. 

For a folder that contains files for each subject that store the data 
from training on the cost-benefit decision-making paradigm that were 
acquired to plot a psychometric function each file is opened and for
each subject a .csv file that contains the different dilutions of
sweetened condensed milk in one row and the corresponding number of
lever presses on the high cost-high benefit (high) and low cost-
low benefit (low) lever is saved in the folder ./data_psychophysics 

The psychometric function can then be plotted using the script
analysis_psyochophysics.py"""

# Import statements
import os
import re

from tkinter import filedialog
from tkinter import *

import numpy as np
import pandas as pd

# Sometimes animals manage to press the lever a second time when it retracts.
# In this case a file may contain data on more than 20 trials. If this is the 
# case, this function determines and returns the number of left, right and 
# omitted lever presses. First, any second lever presses are excluded.
# Because reaction time is recorded from the beginning of the trial rather
# than from the time of lever retraction this can easily be done by
# excluding any lever presses with a reaction time shorter than 10 secs.
# However, if a trial is omitted the reaction time is zero and, hence,
# these trials need to be kept. After cleaning up the data the first 20 trials
# are selected and the number of left, right and omitted lever presses is
# determined and returned
def count_lever_presses(text_data):
    lever_data = re.sub(r'\s+[0-9]+:\s+|\s+', ' ', " ".join(text_data[38:])).split("D:")
            
    lever_presses = np.array(lever_data[0].strip().split(" "), dtype='float32')
    reaction_time = np.array(lever_data[1].strip().split(" "), dtype='float32')
            
    lever_presses_tidy = lever_presses[np.logical_or(reaction_time == 0.0, np.logical_and(reaction_time >=10.0, reaction_time <= 20.0))]
    lever_presses_count = np.bincount(lever_presses_tidy.astype(int)[0:20], minlength = 3)
    
    return lever_presses_count[1], lever_presses_count[2], lever_presses_count[0]

# Store the subject ID, dilution of sweetened condensed milk, group, high, low
# and omitted lever presses. The number of high and low lever presses either
# correspond to the number of presses on the left or right lever depending
# on the group of the subject, i.e. if the high cost-high benefit lever was
# the left or right lever for that specific animal
def append_all(subject, dilution, group, left, right, omit):
        subjects.append(subject)
        dilutions.append(dilution)

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
        if groups[-1] == 2:
            high.append(right)
            low.append(left)
        
        omitted.append(omit)

# Initiate empty lists to append later on
subjects = []
dilutions = []
groups = []
high = []
low = []
omitted = []

# Open a dialog to select the folder that the data are stored in
root = Tk()
root.withdraw()
path =  filedialog.askdirectory(title = "Select directory...")
root.destroy()

# For each file in the selected directory open the file, extract the subject ID,
# dilution of the sweetened condensed milk as provided in the data file, 
# if the high cost-high benefit lever was the left or right lever (group) and
# the number of high, low and omitted lever presses, and store them in lists
for file in os.listdir(path):
    with open(os.path.join(path, file), 'r') as openfile:
        # Read all lines of the open file
        raw_data = openfile.readlines()
        # Extract the subject ID
        subject = re.sub(r'Subject:|\s', '', raw_data[6])

        # Extract the dilution of sweetened condensed milk
        dilution = int(re.sub(r'Experiment: CB-DM-|\s', '', raw_data[7]))
        
        # Extract if the high cost-high benefit lever was the left or right
        # lever for this specific animal as provided in the file (group)
        group = raw_data[8][7:9]
        
        # Extract the number of presses on the left and right lever,
        # and the number of omitted lever presses
        left = int(float(re.sub(r'P:|\s', '', raw_data[26])))        
        right = int(float(re.sub(r'Q:|\s', '', raw_data[27])))        
        omit = int(float(re.sub(r'R:|\s', '', raw_data[28])))
        
        # Calculate the total number of lever presses to check 
        # if the file contains 20 trials
        total = left + right + omit
        
        # If the file contains less than 20 trials raise a value error
        if (total < 20):
            raise ValueError("The file " + path + "/" + file + "contains less than 20 trials.")
         
        # If the file contains exactly 20 trials store the subject ID,
        # the dilution of sweetened condensed milk, the group and the number
        # of high, low and omitted lever presses using the append_all function
        if (total == 20):
            append_all(subject, dilution, group, left, right, omit)
        
        # Sometimes animals manage to press the lever a second time when
        # it retracts. In this case a file may contain data on more than
        # 20 trials. If this is the case, use the count_lever_presses function
        # to extract the number of high, low and omitted trials and
        # use the append_all function after to store the data          
        if (total > 20):
            left, right, omit = count_lever_presses(raw_data)
            append_all(subject, dilution, group, left, right, omit)

# Store the subject ID, dilution of sweetened condensed milk and number of
# high and low lever presses in a data frame and sort by dilution                          
data = pd.DataFrame({"subject": subjects, "dilution": dilutions, "high": high, "low": low})
data = data.sort_values(by=['dilution'])

# For each subject write the data to a .csv file and store it in the folder
# ./data_psychophysics
for subject in subjects:
    data_subject = data.loc[data['subject'] == subject, 'dilution':]
    data_subject.to_csv("./data_psychophysics/psychophysics_" + subject + ".csv", index = False)
                          