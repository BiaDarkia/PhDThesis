#!/usr/bin/env python version 2.7

"""Script to fit a psychophysical function (Weibull function)
to the data collected while animals were trained on the
cost-benefit decision-making paradigm"""

# Import statements
import Tkinter, tkFileDialog
import pandas as pd
import numpy as np
from psychopy.data import FitWeibull
import pylab
import matplotlib.pyplot as plt

# Provide a guess for the threshold and slope of the psychophysical function 
# (Weibull function) that will be fitted to the data
threshold = 5
slope = 2

# Open a file dialog to select a file from an animal that contains that
# animal's behavioral data from cost-benefit decision-making training
root = Tkinter.Tk()
root.withdraw()

file = tkFileDialog.askopenfilename()

# Read the data from the file
cb = pd.read_csv(file, index_col=False)

dilution = cb.loc[:,'dilution'].values.astype(float)
low = cb.loc[:,'low'].values.astype(float)
high = cb.loc[:,'high'].values.astype(float)

total = high+low

# Account for dilutions of sweetened condensed milk that the animal
# omitted all trials by deleting the condition from the data
zero_elements, = np.where(total == 0)

if zero_elements is not 0:
	low = np.delete(low, np.where(total == 0))
	high = np.delete(high, np.where(total == 0))
	dilution = np.delete(dilution, np.where(total == 0))

# Calculate the percentage of lox benefit-low cost choices out of all
# trials that were not omitted	
percentage = low/(high+low)

# Fit the data to that percentage
fit = FitWeibull(dilution, percentage, guess=[threshold, slope], expectedMin=0)
dilution_fit = pylab.arange(0, 20, 0.001)
estimated_fit = fit.eval(dilution_fit)

# Predict the dilution of sweetened condensed milk that will have the animal
# choose the low benefit-low cost option and the high benefit-high cost
# option in about 50% of non-omitted trials and print it
prediction = fit.inverse(0.5)
print(prediction)

# Plot data and psychophysical curve (Weibull function)
plt.suptitle('50 percent level = %f' %prediction, fontsize = 10)
plt.title('Threshold =%6.2f, Slope =%6.2f' %(threshold, slope), fontsize = 8)
plt.axis([-0.1,20,0,1.01])
plt.plot(dilution, percentage, 'ro')
plt.plot(dilution_fit, estimated_fit)
plt.show()

# The plot must be saved manually!!!
