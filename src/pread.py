"""Pandas csv file reading wrapper"""
import pandas as pd
import os
from numpy import exp

def convert_units(age_col, rad_col, logL_col):
    ages = []
    radii = []
    lumin = []

    # convert to years
    for i in range(1, len(age_col)):
        ages.append(float(age_col[i])*1e9)

    # convert to m
    for i in range(1, len(rad_col)):
        radii.append(float(rad_col[i])*6.9911e7)

    # convert to watts
    for i in range(1, len(logL_col)):
        lumin.append((10**(float(logL_col[i])))*3.828e26)

    return ages, radii, lumin

class Pread:

    def __init__(self, dir, subdir, filename):
        df = pd.read_csv(os.path.join(dir, subdir, filename), sep="\s+",
                        engine='python', names=['M', 'age', 'RvRJ', 'logL'])
        age_col = df.age
        rad_col = df.RvRJ
        logL_col = df.logL
        self.ages, self.radii, self.lumin = convert_units(age_col, rad_col, logL_col)

    def find_time_index(self, time):
        for i in range(1, len(self.ages)):
            if self.ages[i-1] <= time <= self.ages[i]:
                return i-1
        return -1
