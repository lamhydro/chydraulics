#!/usr/bin/python3

# -*- coding: utf-8 -*-

import json
import pandas as pd
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from tabulate import tabulate

from . import clib

class GVF():
  """
  Class to estimate the gradually varied flow in a open channel flow
  """
  def __init__(self):

    # Read input data
    with open(sys.argv[1]) as f:
      self._data = json.load(f)
    
    print(self._data)

    # Set gravity
    self._setGravity()

    # Set unit convertion factor
    self._setConFac()

    # Convert degrees angles in radians
    self._conAnglToRad()

    # Set the channel slope angle
    #self._setSoAngle() 

    # Estimate the critical depth
    self._criticalDepth()

    # Estimate the normal depth
    self._normalDepth()

    # Estimate the Froude number
    self._froudeNumber()

    # Classify the profile type
    self._profileType()

    # Compute profile
    self._computeProfile()

    self._plotFlowProfile()

  def _setGravity(self):
    """
    Set the gravity constant
    """
    self._g = clib.gravity(self._data['US'])

  def _setConFac(self):
    """
    Set the convertion factor
    """
    self._conf = clib.convertionFactor(self._data['US'])

  def _conAnglToRad(self):
    """
    Convert angles from degrees to radians
    """

    ns = 0 # number of homogeneous sections
    # Loop through the outer dictionary
    for outer_key, inner_dict in self._data.items():
        #print(f"Outer key: {outer_key}")
        # Loop through the inner dictionary
        if isinstance(inner_dict, dict):  # Check if the value is a dictionary
            ns += 1
            for inner_key, inner_value in inner_dict.items():
                if (inner_key == 'theta1') or (inner_key == 'theta2'):
                    if inner_value != "":
                        #print(f"  {inner_key}: {inner_value}")
                        self._data[outer_key][inner_key] = math.radians(inner_value)  
                        
    self._data["ns"] = ns

#    def _setSoAngle(self):
    #  """
    #  Set a the channel slope angle
    #  """
    #  # Loop through the outer dictionary
    #  for outer_key, inner_dict in self._data.items():
        #  #print(f"Outer key: {outer_key}")
        #  # Loop through the inner dictionary
        #  if isinstance(inner_dict, dict):  # Check if the value is a dictionary
            #  for inner_key, inner_value in inner_dict.items():
                #  if (inner_key == 'So'):
                    #  #print(f"  {inner_key}: {inner_value}")
                    #  self._data[outer_key][inner_key] = clib.slopeAngle(inner_value)

  def _criticalDepth(self):
    """
    Estimate the critical depth
    """
    # Loop through the outer dictionary
    for outer_key, inner_dict in self._data.items():
        # Loop through the inner dictionary
        if isinstance(inner_dict, dict):  # Check if the value is a dictionary
            if self._data[outer_key]["r"] != "":
                self._data[outer_key]["yc"] = clib.criticalDepthC(self._data[outer_key]['So'], self._data[outer_key]['alpha'], self._data[outer_key]['Q'], self._g, self._data[outer_key]['r'])
            else:
                self._data[outer_key]["yc"] = clib.criticalDepth(self._data[outer_key]['So'], self._data[outer_key]['alpha'], self._data[outer_key]['Q'], self._g, self._data[outer_key]['b'], self._data[outer_key]['theta1'], self._data[outer_key]['theta2'])


    print(self._data)


  def _normalDepth(self):
    """
    Estimate the normal depth
    """
    
    # Loop through the outer dictionary
    for outer_key, inner_dict in self._data.items():
        # Loop through the inner dictionary
        if isinstance(inner_dict, dict):  # Check if the value is a dictionary
            if self._data[outer_key]['So'] > 0.0:
                i = 1
                if self._data[outer_key]["r"] != "": # Cicular channel
                    self._a = clib.TI
                    self._b = clib.TF
                    while i<= clib.NMAX:
                        c = (self._a + self._b)/2.0
                        fa = self._data[outer_key]['Q'] - clib.QmanningC(self._conf, self._data[outer_key]['So'], self._data[outer_key]['n'], 2*self._a, self._data[outer_key]['r'])   
                        fb = self._data[outer_key]['Q'] - clib.QmanningC(self._conf, self._data[outer_key]['So'], self._data[outer_key]['n'], 2*self._b, self._data[outer_key]['r'])   
                        fc = self._data[outer_key]['Q'] - clib.QmanningC(self._conf, self._data[outer_key]['So'], self._data[outer_key]['n'], 2*c, self._data[outer_key]['r'])   
                        if abs(fc) < clib.ERROR or (self._b-self._a)*0.5<clib.ERROR:
                          break
                  
                        if np.sign(fc)== np.sign(fa):
                          self._a = c
                        else:
                          self._b = c
                        i += 1  
                    if c==90.0:
                        self._data[outer_key]["yn"] = self._data[outer_key]['r']
                    if c>0.0 and c<90.0:
                        self._data[outer_key]["yn"] = self._data[outer_key]['r']*(1.0 - math.cos(math.radians(c)))
                    if c>90.0 and c<180.0:
                        #self._data[outer_key]["yn"] = self._data[outer_key]['r']*(1.0 + math.sin(math.radians(c)))
                        self._data[outer_key]["yn"] = self._data[outer_key]['r']*(1.0 + math.cos(math.radians(180-c)))
                else:
                    self._a = clib.YI
                    self._b = clib.YF
                    while i<= clib.NMAX:
                        c = (self._a + self._b)/2.0
                        fa = self._data[outer_key]['Q'] - clib.Qmanning(self._conf, self._data[outer_key]['So'], self._data[outer_key]['n'], self._data[outer_key]['b'], self._data[outer_key]['theta1'], self._data[outer_key]['theta2'], self._a)   
                        fb = self._data[outer_key]['Q'] - clib.Qmanning(self._conf, self._data[outer_key]['So'], self._data[outer_key]['n'], self._data[outer_key]['b'], self._data[outer_key]['theta1'], self._data[outer_key]['theta2'], self._b)   
                        fc = self._data[outer_key]['Q'] - clib.Qmanning(self._conf, self._data[outer_key]['So'], self._data[outer_key]['n'], self._data[outer_key]['b'], self._data[outer_key]['theta1'], self._data[outer_key]['theta2'], c)   
                        if abs(fc) < clib.ERROR or (self._b-self._a)*0.5<clib.ERROR:
                          break
                        if np.sign(fc)== np.sign(fa):
                          self._a = c
                        else:
                          self._b = c
                        i += 1  
 
                    self._data[outer_key]["yn"] = c
            else:
                self._data[outer_key]["yn"] = math.nan

    print(self._data)

  def _froudeNumber(self):
    """
    Estimate the Froude number
    """
    # Loop through the outer dictionary
    for outer_key, inner_dict in self._data.items():
        # Loop through the inner dictionary
        if isinstance(inner_dict, dict):  # Check if the value is a dictionary
            if self._data[outer_key]['So'] > 0.0:
                if self._data[outer_key]["r"] != "": # Cicular channel
                    self._data[outer_key]["Fr"] = clib.froudeNumberC(self._data[outer_key]['So'], self._data[outer_key]['alpha'], self._data[outer_key]['Q'], self._g, self._data[outer_key]['r'], self._data[outer_key]['yn'])
                else:
                    self._data[outer_key]["Fr"] = clib.froudeNumber(self._data[outer_key]['So'], self._data[outer_key]['alpha'], self._data[outer_key]['Q'], self._g, self._data[outer_key]['b'], self._data[outer_key]['theta1'], self._data[outer_key]['theta2'], self._data[outer_key]['yn'])
            else:
                self._data[outer_key]["Fr"] = math.nan

    print(self._data)

  def _profileType(self):
    """
    Get the profile type
    """
    # Loop through the outer dictionary
    for outer_key, inner_dict in self._data.items():
        # Loop through the inner dictionary
        if isinstance(inner_dict, dict):  # Check if the value is a dictionary
            if self._data[outer_key]["y1"] == "":
                self._data[outer_key]["y1"] = self._data[outer_key]["yc"]

            self._data[outer_key]["PT"] = clib.classify_gvf_profiles(self._data[outer_key]['y1'], self._data[outer_key]['yn'], self._data[outer_key]['yc'], self._data[outer_key]['So'])

    print(self._data)
    #sys.exit()

  def _computeProfile(self):
    """
    Compute gradually varied flow profile
    """

    def computeZ(So, xs):
        i = 0
        z = 0
        zs = []
        zs.append(z)
        while (i<(len(xs)-1)):
            dx = xs[i]-xs[i+1]
            z = z + So*dx
            zs.append(z)
            i+=1
        return(zs)

    if self._data['MT'] == 'DS': # Direct step method
        # Loop through the outer dictionary
        for outer_key, inner_dict in self._data.items():
            # Loop through the inner dictionary
            if isinstance(inner_dict, dict):  # Check if the value is a dictionary
                
                if self._data[outer_key]["r"] != "": # Cicular channel
                    [ys, xs, As, Rs, Vs, Sfs, Es]= clib.direct_stepC(self._conf, self._data[outer_key]['So'], self._data[outer_key]['alpha'], self._data[outer_key]['Q'], self._g, self._data[outer_key]['r'], self._data[outer_key]['y1'], self._data[outer_key]['x1'], self._data[outer_key]['yn'],self._data[outer_key]['n'], self._data[outer_key]['step'], self._data[outer_key]['PT'], self._data[outer_key]['yc'])
                else:
                    [ys, xs, As, Rs, Vs, Sfs, Es]= clib.direct_step(self._conf, self._data[outer_key]['So'], self._data[outer_key]['alpha'], self._data[outer_key]['Q'], self._g, self._data[outer_key]['b'], self._data[outer_key]['theta1'], self._data[outer_key]['theta2'], self._data[outer_key]['y1'], self._data[outer_key]['x1'], self._data[outer_key]['yn'],self._data[outer_key]['n'], self._data[outer_key]['step'], self._data[outer_key]['PT'], self._data[outer_key]['yc'])
                zs = computeZ(self._data[outer_key]['So'], xs)

                yns = [self._data[outer_key]['yn']] * len(zs)
                ycs = [self._data[outer_key]['yc']] * len(zs)

    elif self._data['MT'] == 'SS': # Standard step method
        
        # Loop through the outer dictionary
        for outer_key, inner_dict in self._data.items():
            # Loop through the inner dictionary
            if isinstance(inner_dict, dict):  # Check if the value is a dictionary
                
                dz = self._data[outer_key]["So"] * self._data[outer_key]["step"]

                if self._data[outer_key]["r"] != "": # Cicular channel
                    print('here')
                else:
                    [ys, xs, zs, As, Rs, Vs, Sfs, Es] = clib.standard_step(self._conf, self._data[outer_key]['Q'], self._g, self._data[outer_key]['theta1'], self._data[outer_key]['theta2'], self._data[outer_key]['n'], self._data[outer_key]['So'], self._data[outer_key]['alpha'], self._data[outer_key]['b'], self._data[outer_key]['y1'], self._data[outer_key]['x1'], self._data[outer_key]['z1'], self._data[outer_key]['step'], dz, self._data[outer_key]['yn'], self._data[outer_key]['yc'], self._data[outer_key]['PT'], self._data[outer_key]['xf'])

                    #zs = computeZ(self._data[outer_key]['So'], xs)

                yns = [self._data[outer_key]['yn']] * len(zs)
                ycs = [self._data[outer_key]['yc']] * len(zs)

                if self._data[outer_key]['PT'] in ['M1','M2','S1','H2']:
                    xs.reverse()
                if self._data[outer_key]['PT'] in ['A2']:
                    xs.reverse()
                    zs.reverse()
    #xs = [-ii for ii in xs]
    #ys.reverse()
    #zs.reverse()
    #print(ys, xs, As, Rs, Vs, Sfs, Es)
    self._df = pd.DataFrame({'y': ys, 'Area': As, 'Rh': Rs, 'Vel':Vs, 'Sf':Sfs, 'Es':Es, 'x':xs, 'z':zs, 'yn':yns, 'yc':ycs})
    # Print the DataFrame nicely
    print(tabulate(self._df, headers='keys', tablefmt='fancy_grid'))

    #sys.exit()


  def _plotFlowProfile(self):
    
    if self._data['US'] == 'IS':
        ylab = 'm'
        xlab = 'm'
    elif self._data['US'] == 'BG':
        ylab = 'ft'
        xlab = 'ft'
     
    x = self._df['x'].tolist()
    y = self._df['y'].tolist()
    z = self._df['z'].tolist()
    yn = self._df['yn'].tolist()
    yc = self._df['yc'].tolist()
    Es = self._df['Es'].tolist()
    
    yr = [i + j for i, j in zip(y, z)]
    ynr = [i + j for i, j in zip(yn, z)]
    ycr = [i + j for i, j in zip(yc, z)]
    Et = [i + j for i, j in zip(Es, z)]

    #print(self._df.dtypes)
    plt.figure(figsize=(15, 5))

    # Plot abline (y = slope * x + intercept)
    #plt.plot(x, y_abline, color='red', linestyle='--', label="Abline", lw=0.5)
    plt.plot(x, yr, color='blue', label="Water surface", lw=1.5)
    plt.plot(x, z, color='black', label="Bottom", lw=1.0)
    plt.plot(x, ynr, color='green', linestyle='--', label="Normal depth", lw=0.5)
    plt.plot(x, ycr, color='red', linestyle='-.', label="Critical depth", lw=0.5)
    #plt.plot(x, Et, color='orange', linestyle='-.', label="Energy line", lw=0.5)
    #plt.plot(self._Es, y_curve, color='black')
    #plt.plot(self._Ec, self._yc, '.', color='red', markersize=10)
   
    # Set square aspect ratio
    #plt.gca().set_aspect('equal', adjustable='box')

    #plt.axhline(0, color='black', lw=0.5, ls='--')  # Reference line at y=0
    #plt.axvline(0, color='black', lw=0.5, ls='--')  # Reference line at x=0
    plt.grid()
    #plt.xlim(0, 1.5*np.max(y_curve))  # Set x-axis limits from 0 to 10
    # Add labels and legend
    plt.xlabel('Horizontal distance (%s)' % xlab)
    plt.ylabel('(%s)' % ylab)
    plt.title("Longitudinal channel flow profile", fontsize=14, color='blue', loc='center')
    plt.legend()
    
    # Show plot
    plt.show()


