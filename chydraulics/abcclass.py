#!/usr/bin/python3

# -*- coding: utf-8 -*-

import json
import pandas as pd
import sys
import numpy as np
import math

from . import clib

class abCoeff():
  """
  Class to compute numerically the energy (alpha) and the momentum (beta) coefficients in channels.
  """
  def __init__(self):

    # Read input data
    with open(sys.argv[1]) as f:
      self._data = json.load(f)

    #print(self._data)

    # Get the wet area
    self._getWetArea()

    # Get the discharge
    self._getDischarge()

    # Get the energy coefficient
    self._getEnergyCoeff()
    
    # Get the momentum coefficient
    self._getMomentCoeff()

  def _getWetArea(self):
    """
    Estimate the wet area and the subsection areas
    """

    self._area, self._subareas = clib.irregCrossSecA(self._data['xs'], self._data['ys'])

    # Printing results
    if self._data['US'] == 'IS':
      print('Wet area (A) = %8.2f m^2' % self._area) 
    elif self._data['US'] == 'BG':
      print('Wet area (A) = %8.2f ft^2' % self._area) 

  def _getDischarge(self):
    """
    Get the discharge
    """
    self._Q = clib.QI(self._data['xs'], self._data['ys'], self._data['v'])

    # Printing results
    if self._data['US'] == 'IS':
      print('Discharge (Q) = %8.2f m^3/s' % self._Q) 
    elif self._data['US'] == 'BG':
      print('Discharge (Q) = %8.2f ft^3/s' % self._Q) 

  def _getEnergyCoeff(self):
    """
    Get the energy coefficient
    """
    self._alpha = clib.energyCoeff(self._data['xs'], self._data['ys'], self._data['v'])
    print('Energy coefficient (alpha) = %8.5f m^3/s' % self._alpha) 

  def _getMomentCoeff(self):
    """
    Get the momentum coefficient
    """
    self._beta = clib.momentCoeff(self._data['xs'], self._data['ys'], self._data['v'])
    print('Energy momentum (beta) = %8.5f m^3/s' % self._beta) 


