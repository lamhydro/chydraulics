#!/usr/bin/python3

# -*- coding: utf-8 -*-

import json
import pandas as pd
import sys
import numpy as np
import math

from . import clib

class UniformFlow():
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
    self._setConFac()

    # Set the cross sectio type 
    #self._setSectionType()

    # Convert degrees angles in radians
    self._conAnglToRad()

    # Set type of problem and solve it
    self._problemType()


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

  #def _setSectionType(self):
  #  """
  #  Define the section type: circular or non circular
  #  """
  #  if self._data['r']!="":
  #    self._data['ST'] = 1 # Circular cross section
  #  else:
  #    self._data['ST'] = 2 # No circular cross section

  def _setIterDomain(self):
    """
    Set the domain of iteration
    """
    if self._data['ST']==1:
      self._a = clib.TI
      self._b = clib.TF
    elif self._data['ST']==2:
      self._a = clib.YI
      self._b = clib.YF


  def _problemType(self):
    if self._data['y'] == "":
      print('')
      print('Find the normal depth (y_n)')
      print('')
      self._getNormalDepth()
    if self._data['Q'] == "":
      print('')
      print('Find the normal discharge (Q)')
      print('')
      self._getNormalDisch()  
    if self._data['So'] == "":
      print('')
      print('Find the channel slope (So)')
      print('')
      self._getChSlope()  
    if self._data['n'] == "":
      print('')
      print('Find Manning roughness factor (n)')
      print('')
      self._getNmanning()  

  def _getNormalDisch(self):
    """
    Estimate the normal discharge
    """
    if self._data['ST'] == 1:
      theta = clib.thetaInC(self._data['r'], self._data['y'])
      self._data['Q'] = clib.QmanningC(self._conf, self._data['So'], self._data['n'], theta, self._data['r'])   
    elif self._data['ST'] == 2:
      self._data['Q'] = clib.Qmanning(self._conf, self._data['So'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], self._data['y'])   

    # Printing results
    if self._data['US'] == 'IS':
      print('Normal discharge (Q) = %8.3f m³/s' % self._data['Q']) 
    elif self._data['US'] == 'BG':
      print('Normal discharge (Q) = %8.3f ft³/s' % self._data['Q']) 


  def _getChSlope(self):
    """
    Estimate the channel slope
    """
    if self._data['ST'] == 1:
      theta = clib.thetaInC(self._data['r'], self._data['y'])
      self._data['So'] = clib.SomanningC(self._conf, self._data['Q'], self._data['n'], theta, self._data['r'])   
    elif self._data['ST'] == 2:
      self._data['So'] = clib.Somanning(self._conf, self._data['Q'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], self._data['y'])   

    # Printing results
    if self._data['US'] == 'IS':
      print('Channel slope (So) = %8.5f m/m' % self._data['So']) 
    elif self._data['US'] == 'BG':
      print('Channel slope (So) = %8.5f ft/ft' % self._data['So']) 


  def _getNmanning(self):
    """
    Estimate the Manning roughness coefficient
    """
    if self._data['ST'] == 1:
      theta = clib.thetaInC(self._data['r'], self._data['y'])
      self._data['n'] = clib.NmanningC(self._conf, self._data['Q'], self._data['So'], theta, self._data['r'])   
    elif self._data['ST'] == 2:
      self._data['n'] = clib.Nmanning(self._conf, self._data['Q'], self._data['So'], self._data['b'], self._data['theta1'], self._data['theta2'], self._data['y'])   

    # Printing results
    if self._data['US'] == 'IS':
      print('Manning roughness (n) = %8.5f' % self._data['n']) 
    elif self._data['US'] == 'BG':
      print('Manning roughness (n) = %8.5f' % self._data['n']) 



  def _getNormalDepth(self):
    """
    Estimate the normal depth
    """
    
    # Set iteration limits
    self._setIterDomain()

    i = 1
    while i<= clib.NMAX:
      c = (self._a + self._b)/2.0
      print(c)

      if self._data['ST'] == 1:
        fa = self._data['Q'] - clib.QmanningC(self._conf, self._data['So'], self._data['n'], 2*self._a, self._data['r'])   
        fb = self._data['Q'] - clib.QmanningC(self._conf, self._data['So'], self._data['n'], 2*self._b, self._data['r'])   
        fc = self._data['Q'] - clib.QmanningC(self._conf, self._data['So'], self._data['n'], 2*c, self._data['r'])   
      elif self._data['ST'] == 2:
        fa = self._data['Q'] - clib.Qmanning(self._conf, self._data['So'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], self._a)   
        fb = self._data['Q'] - clib.Qmanning(self._conf, self._data['So'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], self._b)   
        fc = self._data['Q'] - clib.Qmanning(self._conf, self._data['So'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], c)   

      if abs(fc) < clib.ERROR or (self._b-self._a)*0.5<clib.ERROR:
        break

      if np.sign(fc)== np.sign(fa):
        self._a = c
      else:
        self._b = c
      i += 1  

    if self._data['ST'] == 1: # Estimating y for circular channel
      if c==90.0:
        yn = self._data['r']
      if c>0.0 and c<90.0:
        yn = self._data['r']*(1.0 - math.cos(math.radians(c)))
      if c>90.0 and c<180.0:
        yn = self._data['r']*(1.0 + math.cos(math.radians(180-c)))
    elif self._data['ST'] == 2:
        yn = c
  
    # Printing results
    if self._data['US'] == 'IS':
      print('Normal depth (y_n) = %8.3f m' % yn) 
    elif self._data['US'] == 'BG':
      print('Normal depth (y_n) = %8.3f ft' % yn) 

