#!/usr/bin/python3

# -*- coding: utf-8 -*-

import json
import pandas as pd
import sys
import numpy as np
import math

from . import clib

class UniformFlowI():
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

    # Set type of problem and solve it
    #self._problemType()

    self._getNormalDepth()

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

#  def _getYmax(self):
#    """
#    Get the maximun y coord.
#    """
#    self._ymax = max(self._data['ys'])
#
#  def _getYmin(self):
#    """
#    Get the minimun y coord.
#    """
#    self._ymin = min(self._data['ys'])
#
#  def _getXmax(self):
#    """
#    Get the maximun X coord.
#    """
#    self._xmax = max(self._data['xs'])
#
#  def _getXmin(self):
#    """
#    Get the minimun x coord.
#    """
#    self._xmin = min(self._data['xs'])


  def _getNormalDepth(self):
    """
    Estimate the normal depth
    """
    
    # Set iteration limits
    self._a = min(self._data['ys']) + clib.YI
    self._b = min(self._data['ys']) + clib.YF

    i = 1
    while i<= clib.NMAX:
      c = (self._a + self._b)/2.0
      #print(i, c, self._a, self._b, 'hereeeeeee')

      xsn, ysn, nsn = clib.interpAnewSection(self._data['xs'], self._data['ys'], self._data['ns'], self._a)
      fa = self._data['Q'] - clib.QmanningI(self._conf, self._data['So'], nsn, xsn, ysn)   
      xsn, ysn, nsn = clib.interpAnewSection(self._data['xs'], self._data['ys'], self._data['ns'], self._b)
      fb = self._data['Q'] - clib.QmanningI(self._conf, self._data['So'], nsn, xsn, ysn)   
      xsn, ysn, nsn = clib.interpAnewSection(self._data['xs'], self._data['ys'], self._data['ns'], c)
      fc = self._data['Q'] - clib.QmanningI(self._conf, self._data['So'], nsn, xsn, ysn)   

      if abs(fc) < clib.ERROR or (self._b-self._a)*0.5<clib.ERROR:
        break

      if np.sign(fc)== np.sign(fa):
        self._a = c
      else:
        self._b = c
      i += 1  
 
    # Printing results
    if self._data['US'] == 'IS':
      print('Normal depth (y_n) = %8.3f m' % c) 
    elif self._data['US'] == 'BG':
      print('Normal depth (y_n) = %8.3f ft' % c) 


#  def _problemType(self):
#    if self._data['y'] == "":
#      print('')
#      print('Find the normal depth (y_n)')
#      print('')
#      self._getNormalDepth()
#    if self._data['Q'] == "":
#      print('')
#      print('Find the normal discharge (Q)')
#      print('')
#      self._getNormalDisch()  
#    if self._data['So'] == "":
#      print('')
#      print('Find the channel slope (So)')
#      print('')
#      self._getChSlope()  
#    if self._data['n'] == "":
#      print('')
#      print('Find Manning roughness factor (n)')
#      print('')
#      self._getNmanning()  
#
#  def _getNormalDisch(self):
#    """
#    Estimate the normal discharge
#    """
#    if self._data['ST'] == 1:
#      theta = clib.thetaInC(self._data['r'], self._data['y'])
#      self._data['Q'] = clib.QmanningC(self._conf, self._data['So'], self._data['n'], theta, self._data['r'])   
#    elif self._data['ST'] == 2:
#      self._data['Q'] = clib.Qmanning(self._conf, self._data['So'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], self._data['y'])   
#
#    # Printing results
#    if self._data['US'] == 'IS':
#      print('Normal discharge (Q) = %8.3f m³/s' % self._data['Q']) 
#    elif self._data['US'] == 'BG':
#      print('Normal discharge (Q) = %f8.3f ft³/s' % self._data['Q']) 
#
#
#  def _getChSlope(self):
#    """
#    Estimate the channel slope
#    """
#    if self._data['ST'] == 1:
#      theta = clib.thetaInC(self._data['r'], self._data['y'])
#      self._data['So'] = clib.SomanningC(self._conf, self._data['Q'], self._data['n'], theta, self._data['r'])   
#    elif self._data['ST'] == 2:
#      self._data['So'] = clib.Somanning(self._conf, self._data['Q'], self._data['n'], self._data['b'], self._data['theta1'], self._data['theta2'], self._data['y'])   
#
#    # Printing results
#    if self._data['US'] == 'IS':
#      print('Channel slope (So) = %8.5f m/m' % self._data['So']) 
#    elif self._data['US'] == 'BG':
#      print('Channel slope (So) = %8.5f ft/ft' % self._data['So']) 
#
#
#  def _getNmanning(self):
#    """
#    Estimate the Manning roughness coefficient
#    """
#    if self._data['ST'] == 1:
#      theta = clib.thetaInC(self._data['r'], self._data['y'])
#      self._data['n'] = clib.NmanningC(self._conf, self._data['Q'], self._data['So'], theta, self._data['r'])   
#    elif self._data['ST'] == 2:
#      self._data['n'] = clib.Nmanning(self._conf, self._data['Q'], self._data['So'], self._data['b'], self._data['theta1'], self._data['theta2'], self._data['y'])   
#
#    # Printing results
#    if self._data['US'] == 'IS':
#      print('Manning roughness (n) = %8.5f' % self._data['n']) 
#    elif self._data['US'] == 'BG':
#      print('Manning roughness (n) = %8.5f' % self._data['n']) 




