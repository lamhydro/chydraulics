#!/usr/bin/python3

# -*- coding: utf-8 -*-

import json
import pandas as pd
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from . import clib

class SpecificEnergy():
  """
  Class to estimate the normal depth in a open channel flow
  """
  def __init__(self):

    # Read input data
    with open(sys.argv[1]) as f:
      self._data = json.load(f)
    
    # Set gravity
    self._setGravity()

    # Set unit convertion factor
    #self._setConFac()

    # Set the cross sectio type 
    #self._setSectionType()

    # Convert degrees angles in radians
    self._conAnglToRad()

    # Set type of problem and solve it
    #self._problemType()

    # Set the channel slope angle
    self._setSoAngle() 

    # Set a list of water depth values
    self._setyr()
    
    if len(self._data['b']) == 1 or self._data['b']=="":

        # Plot section
        self._plotSection()
    
        # Estimate specific energy curve
        self._specificEnergyCurve()
    
        # Estimate specific energy curve
        self._specificForceCurve()
    
    else:
        if self._data['ST'] == 2:
            # Estimate specific energy curve
            self._specificEnergyCurve2()
    
            # Estimate specific energy curve
            self._specificForceCurve2()
    

    # Plote specific energy curve
    self._plotSpecificEnergyCurve()

    # Estimate specific energy curve
    self._plotSpecificForceCurve()

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
    if self._data['theta1'] !="":
      self._data['theta1'] = math.radians(self._data['theta1'])  
    if self._data['theta2'] !="":
      self._data['theta2'] = math.radians(self._data['theta2'])  

  def _setSoAngle(self):
    """
    Set a the channel slope angle
    """
    self._thetaS = clib.slopeAngle(self._data['So'])

  def _setyr(self):
    """
    Set a list of water depth values
    """
    step = self._data['step']
    self._yr = np.arange(step, 0.99*self._data['y'], step).tolist()

  def _specificForceCurve(self):
    """
    Estimate the speficic force curve
    """
    self._Fs = list()
    for Q in self._data['Q']:
        if self._data['ST'] == 1:
          self._Fs.append(clib.specificForceCurveC(self._yr, self._thetaS, self._data['beta'], Q, self._g, self._data['r']))
        elif self._data['ST'] == 2:
          self._Fs.append(clib.specificForceCurve(self._yr, self._thetaS, self._data['beta'], Q, self._g, self._data['b'][0], self._data['theta1'], self._data['theta2']))
        elif self._data['ST'] == 3:
          self._Fs.append(clib.specificForceCurveI(self._yr, self._thetaS, self._data['beta'], Q, self._g, self._data['xs'], self._data['ys']))
 
    #print(self._Fs)

  def _specificForceCurve2(self):
    """
    Estimate the speficic force curve for variable width
    """
    self._Fs2 = list()
    for b in self._data['b']:
        self._Fs2.append(clib.specificForceCurve(self._yr, self._thetaS, self._data['beta'], self._data['Q'], self._g, b, self._data['theta1'], self._data['theta2']))
          
       
  def _specificEnergyCurve(self):
    """
    Estimate the speficic energy curve
    """

    self._Es = list()
    for Q in self._data['Q']:
        if self._data['ST'] == 1:
          self._Es.append(clib.specificEnergyCurveC(self._yr, self._thetaS, self._data['alpha'], Q, self._g, self._data['r']))
        elif self._data['ST'] == 2:
          self._Es.append(clib.specificEnergyCurve(self._yr, self._thetaS, self._data['alpha'], Q, self._g, self._data['b'][0], self._data['theta1'], self._data['theta2']))
        elif self._data['ST'] == 3:
          self._Es.append(clib.specificEnergyCurveI(self._yr, self._thetaS, self._data['alpha'], Q, self._g, self._data['xs'], self._data['ys']))
    #print(self._Es)

  def _specificEnergyCurve2(self):
    """
    Estimate the speficic energy curve for variable width
    """

    self._Es2 = list()
    for b in self._data['b']:
          self._Es2.append(clib.specificEnergyCurve(self._yr, self._thetaS, self._data['alpha'], self._data['Q'], self._g, b, self._data['theta1'], self._data['theta2']))

  def _plotSpecificForceCurve(self):
    
    if self._data['US'] == 'IS':
        if len(self._data['b'])==1 or self._data['b']=="":
            lab = 'm3/s'
        else:
            lab = 'm'
        ylab = 'm'
        xlab = 'm^3'
    elif self._data['US'] == 'BG':
        if len(self._data['b'])==1 or self._data['b']=="":
            lab = 'ft3/s'
        else:
            lab = 'm'
        ylab = 'ft'
        xlab = 'ft^3'
     
#      # Define abline parameters
    #  slope = self._thetaS
    #  intercept = 0

    y_curve = np.array(self._yr) #np.sin(x)  # Example curve
    #  #y_abline = math.cos(math.radians(slope)) * x + intercept
    #  #x_abline = (y_curve-intercept)/math.cos(math.radians(slope))
    #  x_abline = y_curve*math.cos(math.radians(slope)) + intercept

    # Create a color map
    if len(self._data['b'])==1 or self._data['b']=="":
        colors = cm.viridis(np.linspace(0, 1, len(self._data['Q'])))  # Use 'viridis' colormap for n colors
    else:
        colors = cm.viridis(np.linspace(0, 1, len(self._data['b'])))  # Use 'viridis' colormap for n colors

    plt.figure(figsize=(10, 6))

    # Plot abline (y = slope * x + intercept)
    #plt.plot(x, y_abline, color='red', linestyle='--', label="Abline", lw=0.5)
   # plt.plot(x_abline, y_curve, color='red', linestyle='--', label="Abline", lw=0.5)
    if  len(self._data['b'])==1 or self._data['b']=="":
        fss = self._Fs
        vari = self._data['Q']
        lenv = len(self._data['Q'])
    else:
        fss = self._Fs2
        vari = self._data['b']
        lenv = len(self._data['b'])

    for Fs,var,i in zip(fss, vari,range(lenv)):
        # Generate example data
        x = np.array(Fs) #np.linspace(0, 10, 100)
        
        # Plot the curve
        plt.plot(x, y_curve, label=str(var)+' ('+lab+')', color=colors[i])
    
    # Set square aspect ratio
    #plt.gca().set_aspect('equal', adjustable='box')

    plt.axhline(0, color='black', lw=0.5, ls='--')  # Reference line at y=0
    plt.axvline(0, color='black', lw=0.5, ls='--')  # Reference line at x=0
    plt.grid()
    #plt.xlim(0, 2.5*np.max(y_curve))  # Set x-axis limits from 0 to 10
    plt.xlim(0, 2000)  # Set x-axis limits from 0 to 10
    # Add labels and legend
    plt.xlabel('Specific Force (%s)' % xlab)
    plt.ylabel('Water depth (%s)' % ylab)
    plt.legend()
    
    # Show plot
    plt.show()

     
  def _plotSpecificEnergyCurve(self):
    
    if self._data['US'] == 'IS':
        if len(self._data['b'])==1 or self._data['b']=="":
            lab = 'm3/s'
        else:
            lab = 'm'
        ylab = 'm'
        xlab = 'm'
    elif self._data['US'] == 'BG':
        if len(self._data['b'])==1 or self._data['b']=="":
            lab = 'ft3/s'
        else:
            lab = 'ft'
        ylab = 'ft'
        xlab = 'ft'
     
    # Define abline parameters
    slope = self._thetaS
    intercept = 0

    y_curve = np.array(self._yr) #np.sin(x)  # Example curve
    #y_abline = math.cos(math.radians(slope)) * x + intercept
    #x_abline = (y_curve-intercept)/math.cos(math.radians(slope))
    x_abline = y_curve*math.cos(math.radians(slope)) + intercept

    # Create a color map
    if len(self._data['b'])==1 or self._data['b']=="":
        colors = cm.viridis(np.linspace(0, 1, len(self._data['Q'])))  # Use 'viridis' colormap for n colors
    else:
        colors = cm.viridis(np.linspace(0, 1, len(self._data['b'])))  # Use 'viridis' colormap for n colors

    plt.figure(figsize=(10, 6))

    # Plot abline (y = slope * x + intercept)
    #plt.plot(x, y_abline, color='red', linestyle='--', label="Abline", lw=0.5)
    plt.plot(x_abline, y_curve, color='red', linestyle='--', label="Abline", lw=0.5)
    if  len(self._data['b'])==1 or self._data['b']=="":
        ess = self._Es
        vari = self._data['Q']
        lenv = len(self._data['Q'])
    else:
        ess = self._Es2
        vari = self._data['b']
        lenv = len(self._data['b'])

    for Es,var,i in zip(ess, vari,range(lenv)):
        # Generate example data
        x = np.array(Es) #np.linspace(0, 10, 100)
        
        # Plot the curve
        plt.plot(x, y_curve, label=str(var)+' ('+lab+')', color=colors[i])
    
    # Set square aspect ratio
    plt.gca().set_aspect('equal', adjustable='box')

    plt.axhline(0, color='black', lw=0.5, ls='--')  # Reference line at y=0
    plt.axvline(0, color='black', lw=0.5, ls='--')  # Reference line at x=0
    plt.grid()
    plt.xlim(0, 1.5*np.max(y_curve))  # Set x-axis limits from 0 to 10
    # Add labels and legend
    plt.xlabel('Specific energy (%s)' % xlab)
    plt.ylabel('Water depth (%s)' % ylab)
    plt.legend()
    
    # Show plot
    plt.show()

  def _plotSection(self):
    """
    Estimate the coord x and y of the centroid of a polygon 
    Where:
      vertices: Tuple with the coord x and y of the polygon area
    """
    if self._data['US'] == 'IS':
        lab = 'm'
    elif self._data['US'] == 'BG':
        lab = 'ft'
  
    if self._data['ST'] == 1:
        center_x, center_y = 0, 0  # Center of the circle
        radius = self._data['r']  # Radius of the circle
        # Generate points on the circle
        t = np.linspace(0, 2 * np.pi, 100)
        x = center_x + radius * np.cos(t)
        y = center_y + radius * np.sin(t)
    elif self._data['ST'] == 2:
        x, y = clib.getXYcoordSec(self._data['b'][0], self._data['theta1'], self._data['theta2'], self._data['y'])
        x = np.array(x)
        y = np.array(y)
    elif self._data['ST'] == 3:
        x = np.array(self._data['xs'])
        y = np.array(self._data['ys'])
    
    plt.plot(x, y, marker='o', linestyle='-', color='b')
  
    # Add labels and title
    plt.xlabel('X'+'('+lab+')')
    plt.ylabel('Y'+'('+lab+')')
    plt.title('Channel cross section')
    plt.grid(True)  # Add a grid for better readability
    #plt.gca().set_aspect('equal', adjustable='box')
    
    # Show plot
    plt.show()


#    #def _setSectionType(self):
  #  #  """
  #  #  Define the section type: circular or non circular
  #  #  """
  #  #  if self._data['r']!="":
  #  #    self._data['ST'] = 1 # Circular cross section
  #  #  else:
  #  #    self._data['ST'] = 2 # No circular cross section

  #  def _setIterDomain(self):
    #  """
    #  Set the domain of iteration
    #  """
    #  if self._data['ST']==1:
      #  self._a = clib.TI
      #  self._b = clib.TF
    #  elif self._data['ST']==2:
      #  self._a = clib.YI
      #  self._b = clib.YF


  #  def _problemType(self):
    #  if self._data['y'] == "":
      #  print('')
      #  print('Find the normal depth (y_n)')
      #  print('')
      #  self._getNormalDepth()
    #  if self._data['Q'] == "":
      #  print('')
      #  print('Find the normal discharge (Q)')
      #  print('')
      #  self._getNormalDisch()  
    #  if self._data['So'] == "":
      #  print('')
      #  print('Find the channel slope (So)')
      #  print('')
      #  self._getChSlope()  
    #  if self._data['n'] == "":
      #  print('')
      #  print('Find Manning roughness factor (n)')
      #  print('')
      #  self._getNmanning()  

  #  def _getNormalDisch(self):
    #  """
    #  Estimate the normal discharge
    #  """
    #  if self._data['ST'] == 1:
      #  theta = clib.thetaInC(self._data['r'], self._data['y'])
      #  self._data['Q'] = clib.QmanningC(self._conf, self._data['So'], self._data['n'], theta, self._data['r'])   
    #  elif self._data['ST'] == 2:
      #  self._data['Q'] = clib.Qmanning(self._conf, self._data['So'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], self._data['y'])   

    #  # Printing results
    #  if self._data['US'] == 'IS':
      #  print('Normal discharge (Q) = %8.3f m³/s' % self._data['Q']) 
    #  elif self._data['US'] == 'BG':
      #  print('Normal discharge (Q) = %8.3f ft³/s' % self._data['Q']) 


  #  def _getChSlope(self):
    #  """
    #  Estimate the channel slope
    #  """
    #  if self._data['ST'] == 1:
      #  theta = clib.thetaInC(self._data['r'], self._data['y'])
      #  self._data['So'] = clib.SomanningC(self._conf, self._data['Q'], self._data['n'], theta, self._data['r'])   
    #  elif self._data['ST'] == 2:
      #  self._data['So'] = clib.Somanning(self._conf, self._data['Q'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], self._data['y'])   

    #  # Printing results
    #  if self._data['US'] == 'IS':
      #  print('Channel slope (So) = %8.5f m/m' % self._data['So']) 
    #  elif self._data['US'] == 'BG':
      #  print('Channel slope (So) = %8.5f ft/ft' % self._data['So']) 


  #  def _getNmanning(self):
    #  """
    #  Estimate the Manning roughness coefficient
    #  """
    #  if self._data['ST'] == 1:
      #  theta = clib.thetaInC(self._data['r'], self._data['y'])
      #  self._data['n'] = clib.NmanningC(self._conf, self._data['Q'], self._data['So'], theta, self._data['r'])   
    #  elif self._data['ST'] == 2:
      #  self._data['n'] = clib.Nmanning(self._conf, self._data['Q'], self._data['So'], self._data['b'], self._data['theta1'], self._data['theta2'], self._data['y'])   

    #  # Printing results
    #  if self._data['US'] == 'IS':
      #  print('Manning roughness (n) = %8.5f' % self._data['n']) 
    #  elif self._data['US'] == 'BG':
      #  print('Manning roughness (n) = %8.5f' % self._data['n']) 



  #  def _getNormalDepth(self):
    #  """
    #  Estimate the normal depth
    #  """
    
    #  # Set iteration limits
    #  self._setIterDomain()

    #  i = 1
    #  while i<= clib.NMAX:
      #  c = (self._a + self._b)/2.0
      #  print(c)

      #  if self._data['ST'] == 1:
        #  fa = self._data['Q'] - clib.QmanningC(self._conf, self._data['So'], self._data['n'], self._a, self._data['r'])   
        #  fb = self._data['Q'] - clib.QmanningC(self._conf, self._data['So'], self._data['n'], self._b, self._data['r'])   
        #  fc = self._data['Q'] - clib.QmanningC(self._conf, self._data['So'], self._data['n'], c, self._data['r'])   
      #  elif self._data['ST'] == 2:
        #  fa = self._data['Q'] - clib.Qmanning(self._conf, self._data['So'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], self._a)   
        #  fb = self._data['Q'] - clib.Qmanning(self._conf, self._data['So'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], self._b)   
        #  fc = self._data['Q'] - clib.Qmanning(self._conf, self._data['So'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], c)   

      #  if abs(fc) < clib.ERROR or (self._b-self._a)*0.5<clib.ERROR:
        #  break

      #  if np.sign(fc)== np.sign(fa):
        #  self._a = c
      #  else:
        #  self._b = c
      #  i += 1  

    #  if self._data['ST'] == 1: # Estimating y for circular channel
      #  if c==90.0:
        #  yn = self._data['r']
      #  if c>0.0 and c<90.0:
        #  yn = self._data['r']*(1.0 - math.cos(math.radians(c)))
      #  if c>90.0 and c<180.0:
        #  yn = self._data['r']*(1.0 + math.sin(math.radians(c)))
    #  elif self._data['ST'] == 2:
        #  yn = c
  
    #  # Printing results
    #  if self._data['US'] == 'IS':
      #  print('Normal depth (y_n) = %8.3f m' % yn) 
    #  elif self._data['US'] == 'BG':
      #  print('Normal depth (y_n) = %8.3f ft' % yn) 

