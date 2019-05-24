"""Pandas csv file reading wrapper"""
import pandas as pd
import os
from numpy import exp
import numpy as np
from scipy.interpolate import interp1d
 
Rjup = 6.9911e7
Lsun = 3.828e26
Mearth = 5.97e24

def convert_units(age_col, rad_col, logL_col):
    ages = []
    radii = []
    lumin = []

    # convert to years
    for i in range(1, len(age_col)):
        ages.append(float(age_col[i])*1e9)

    # convert to m
    for i in range(1, len(rad_col)):
        radii.append(float(rad_col[i])*Rjup)

    # convert to watts
    for i in range(1, len(logL_col)):
        lumin.append((10**(float(logL_col[i])))*Lsun)

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

class BaraffeModelFixedTime():
    def __init__(self, Z, time):
        ''' Provides interpolation functions that allow you to extract planet luminosity and radius from Baraffe 2008 models as a function of planet mass.
        This class provides L and R functions at a fixed time set at initialization.

        === Arguments ===
        Z: fraction of heavy metals. Must be one of '002' = 0.02, '010' = 0.1, '050' = 0.5, '090' = 0.9
        time: Time since formation (in years)
        '''
        if Z=='002':
            self.masses = np.array([20, 100, 318, 636, 1590, 3180])
        elif Z=='010':
            self.masses = np.array([20, 100, 318, 636])
        elif Z=='050':
            self.masses = np.array([20, 100, 318])
        elif Z=='090':
            self.masses = np.array([20, 100, 318])
        else:
            raise AttributeError("Z not in Baraffe list. Must be '002', '010', '050' or '090'")

        Ls, Rs = np.zeros(len(self.masses)), np.zeros(len(self.masses))
        for i, mass in enumerate(self.masses):
            data = Pread('data', 'Z'+Z, 'pltlum_M'+str(mass)+'Z'+Z+'.csv') 
            if time < data.ages[0]:
                raise Warning("Asked for planet model at a time earlier than available in Baraffe 2008 grid. Outputting values at time = {0}".format(data.ages[0]))
                time = data.ages[0]
            if time > data.ages[-1]:
                raise Warning("Asked for planet model at a time later than available in Baraffe 2008 grid. Outputting values at time = {0}".format(data.ages[-1]))
                time = data.ages[-1]

            index = data.find_time_index(time)          
            Ls[i] = data.lumin[index]
            Rs[i] = data.radii[index]

        self.logL = interp1d(np.log10(self.masses*Mearth), np.log10(Ls), kind='linear') # interp going from log M (in kg) to log L (in SI)
        self.logR = interp1d(np.log10(self.masses*Mearth), np.log10(Rs), kind='linear') # interp going from log M (in kg) to log R (in m)

    def L(self, M):
        ''' Takes a mass or array of masses (in kg) and returns a luminosity (in W) or array of luminosities for the time and heavy metals fraction Z that the class was initialized with.'''
        try:
            return 10**self.logL(np.log10(M))
        except:
            raise AttributeError("mass value outside Baraffe grid. Valid masses between {0} and {1} Earth masses".format(self.masses[0], self.masses[-1]))

    def R(self, M):
        ''' Takes a mass or array of masses (in kg) and returns a radius (in m) or array of radii for the time and heavy metals fraction Z that the class was initialized with.'''
        try:
            return 10**self.logR(np.log10(M))
        except:
            raise AttributeError("mass value outside Baraffe grid. Valid masses between {0} and {1} Earth masses".format(self.masses[0], self.masses[-1]))

class BaraffeModelFixedMass():
    def __init__(self, Z, mass):
        ''' Provides interpolation functions that allow you to extract planet luminosity and radius from Baraffe 2008 models as a function of time.
        This class provides L and R functions at a fixed mass set at initialization.

        === Arguments ===
        Z: fraction of heavy metals. Must be one of '002' = 0.02, '010' = 0.1, '050' = 0.5, '090' = 0.9
        mass: planet mass (in Earth masses!). Allowed masses depend on Z. For 0.02, m = [20, 100, 318, 636, 1590, 3180], for 0.1 m= [20, 100, 318, 636], for 0.5 and 0.9, m=[20, 100, 318]
        '''
        try:
            data = Pread('data', 'Z'+Z, 'pltlum_M'+str(int(mass))+'Z'+Z+'.csv') 
        except:
            raise AttributeError("Could not find data file for the passed Z or planet mass. See docstring for allowed values")

        self.logL = interp1d(np.log10(data.ages), np.log10(data.lumin), kind='linear') # interp going from log t (in yrs) to log L (in SI)
        self.logR = interp1d(np.log10(data.ages), np.log10(data.radii), kind='linear') # interp going from log t (in yrs) to log R (in m)

    def L(self, t):
        ''' Takes a time or array of times (in years) and returns a luminosity (in W) or array of luminosities'''
        try:
            return 10**self.logL(np.log10(t))
        except:
            raise AttributeError("mass value outside Baraffe grid. Valid masses between {0} and {1} Earth masses".format(self.masses[0], self.masses[-1]))

    def R(self, t):
        ''' Takes a time or array of times (in years) and returns a radius (in m) or array of radii'''
        try:
            return 10**self.logR(np.log10(t))
        except:
            raise AttributeError("mass value outside Baraffe grid. Valid masses between {0} and {1} Earth masses".format(self.masses[0], self.masses[-1]))
