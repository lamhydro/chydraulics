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

class CriticalFlow():
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
    if self._data['ST'] == 2:
        self._conAnglToRad()

    # Set type of problem and solve it
    #self._problemType()

    # Set the channel slope angle
    self._setSoAngle() 


    # Plot the cross section
    self._plotSection()

    # Estimate the critical depth
    self._criticalDepth()

    # Estimate the specific energy
    self._specificEnergy()

    # Estimate the specific energy
    self._specificForce()

    # Set a list of water depth values
    self._setyr()
        
    # Estimate specific energy curve
    self._specificEnergyCurve()
 
    # Plot specific energy curve
    self._plotSpecificEnergyCurve()

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

  def _criticalDepth(self):
    """
    Estimate the critical depth
    """
    if self._data['ST'] == 1:
       self._yc = clib.criticalDepthC(self._thetaS, self._data['alpha'], self._data['Q'], self._g, self._data['r'])
    elif self._data['ST'] == 2:
       self._yc = clib.criticalDepth(self._thetaS, self._data['alpha'], self._data['Q'], self._g, self._data['b'], self._data['theta1'], self._data['theta2'])
    elif self._data['ST'] == 3:
       self._yc = clib.criticalDepthI(self._thetaS, self._data['alpha'], self._data['Q'], self._g, self._data['xs'], self._data['ys'])
       #  self._Fs.append(clib.specificForceCurveI(self._yr, self._thetaS, self._data['beta'], Q, self._g, self._data['xs'], self._data['ys']))

    # Printing results
    if self._data['US'] == 'IS':
      print('Critical depth (y_c) = %8.3f m' % self._yc) 
    elif self._data['US'] == 'BG':
      print('Critical depth (y_c) = %8.3f ft' % self._yc) 

  def _specificEnergy(self):
    """
    Estimate the critical energy 
    """
    if self._data['ST'] == 1:
       theta = clib.thetaInC(self._data['r'], self._yc)
       A = clib.areaC(theta, self._data['r'])
    elif self._data['ST'] == 2:
       A = clib.area(self._data['b'], self._data['theta1'], self._data['theta2'], self._yc)
    elif self._data['ST'] == 3:
       xsn, ysn, dummy = clib.interpAnewSection(self._data['xs'], self._data['ys'], self._data['xs'], self._yc)
       A,_ = clib.irregCrossSecA(xsn, ysn)
       #  self._Fs.append(clib.specificForceCurveI(self._yr, self._thetaS, self._data['beta'], Q, self._g, self._data['xs'], self._data['ys']))
    self._Ec = clib.specificEnergy(self._yc, self._thetaS, self._data['alpha'], self._data['Q'], self._g, A)

    # Printing results
    if self._data['US'] == 'IS':
      print('--Critical energy (E_c) = %8.3f m' % self._Ec) 
    elif self._data['US'] == 'BG':
      print('--Critical energy (E_c) = %8.3f ft' % self._Ec) 

  def _specificForce(self):
    """
    Estimate the critical force 
    """
    if self._data['ST'] == 1:
       theta = clib.thetaInC(self._data['r'], self._yc)
       A = clib.areaC(theta, self._data['r'])
       z = clib.cirseg_centroid(self._data['r'],theta,self._yc)
       B = clib.surfaceWidthC(theta, self._data['r'])
       P = clib.perimeterC(theta, self._data['r'])
    elif self._data['ST'] == 2:
       A = clib.area(self._data['b'], self._data['theta1'], self._data['theta2'], self._yc)
       xs,ys = clib.getXYcoordSec(self._data['b'], self._data['theta1'], self._data['theta2'], self._yc)
       _,z = clib.polygon_centroid(xs,ys)
       B = clib.surfaceWidth(self._data['b'], self._data['theta1'], self._data['theta2'], self._yc)
       P = clib.perimeter(self._data['b'], self._data['theta1'], self._data['theta2'], self._yc)
    elif self._data['ST'] == 3:
       xsn, ysn, dummy = clib.interpAnewSection(self._data['xs'], self._data['ys'], self._data['xs'], self._yc)
       A,_ = clib.irregCrossSecA(xsn, ysn)
       _,z = clib.polygon_centroid(xsn,ysn)
       B = clib.surfaceWidthI(xsn, ysn)
       P,_ = clib.irregCrossSecP(xsn, ysn)
       #  self._Fs.append(clib.specificForceCurveI(self._yr, self._thetaS, self._data['beta'], Q, self._g, self._data['xs'], self._data['ys']))
    #self._Fsc = clib.specificEnergy(self._yc, self._thetaS, self._data['alpha'], self._data['Q'], self._g, A)
    self._Fsc = clib.specificForce(self._data['beta'], self._data['Q'], self._g, A, z)

    # Printing results
    if self._data['US'] == 'IS':
      print('--Critical force (F_sc) = %8.3f m^3' % self._Fsc) 
      print('--Wet area (A) = %8.3f m^2' % A) 
      print('--Wet perimeter (P) = %8.3f m' % P) 
      print('--Water surface width (B) = %8.3f m' % B) 
    elif self._data['US'] == 'BG':
      print('--Critical force (F_sc) = %8.3f ft^3' % self._Fsc) 
      print('--Wet area (A) = %8.3f ft^2' % A) 
      print('--Wet perimeter (P) = %8.3f ft' % P) 
      print('--Water surface width (B) = %8.3f m' % B) 

  def _setyr(self):
    """
    Set a list of water depth values
    """
    step = self._data['step']
    self._yr = np.arange(step, 0.99*self._data['y'], step).tolist()


  def _specificEnergyCurve(self):
    """
    Estimate the speficic energy curve
    """

    self._Es = list()
    if self._data['ST'] == 1:
      self._Es = clib.specificEnergyCurveC(self._yr, self._thetaS, self._data['alpha'], self._data['Q'], self._g, self._data['r'])
    elif self._data['ST'] == 2:
      self._Es = clib.specificEnergyCurve(self._yr, self._thetaS, self._data['alpha'], self._data['Q'], self._g, self._data['b'], self._data['theta1'], self._data['theta2'])
    elif self._data['ST'] == 3:
      self._Es = clib.specificEnergyCurveI(self._yr, self._thetaS, self._data['alpha'], self._data['Q'], self._g, self._data['xs'], self._data['ys'])
 
  def _plotSpecificEnergyCurve(self):
    
    if self._data['US'] == 'IS':
        ylab = 'm'
        xlab = 'm'
    elif self._data['US'] == 'BG':
        ylab = 'ft'
        xlab = 'ft'
     
    # Define abline parameters
    slope = self._thetaS
    intercept = 0

    y_curve = np.array(self._yr) #np.sin(x)  # Example curve
    #y_abline = math.cos(math.radians(slope)) * x + intercept
    #x_abline = (y_curve-intercept)/math.cos(math.radians(slope))
    x_abline = y_curve*math.cos(math.radians(slope)) + intercept

    plt.figure(figsize=(10, 6))

    # Plot abline (y = slope * x + intercept)
    #plt.plot(x, y_abline, color='red', linestyle='--', label="Abline", lw=0.5)
    plt.plot(x_abline, y_curve, color='blue', linestyle='--', label="Abline", lw=0.5)
    plt.plot(self._Es, y_curve, color='black')
    plt.plot(self._Ec, self._yc, '.', color='red', markersize=10)
   
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
        x, y = clib.getXYcoordSec(self._data['b'], self._data['theta1'], self._data['theta2'], self._data['y'])
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



#    def _specificForceCurve(self):
    #  """
    #  Estimate the speficic force curve
    #  """
    #  self._Fs = list()
    #  for Q in self._data['Q']:
        #  if self._data['ST'] == 1:
          #  self._Fs.append(clib.specificForceCurveC(self._yr, self._thetaS, self._data['beta'], Q, self._g, self._data['r']))
        #  elif self._data['ST'] == 2:
          #  self._Fs.append(clib.specificForceCurve(self._yr, self._thetaS, self._data['beta'], Q, self._g, self._data['b'][0], self._data['theta1'], self._data['theta2']))
        #  elif self._data['ST'] == 3:
          #  self._Fs.append(clib.specificForceCurveI(self._yr, self._thetaS, self._data['beta'], Q, self._g, self._data['xs'], self._data['ys']))
 
    #  #print(self._Fs)

  #  def _specificForceCurve2(self):
    #  """
    #  Estimate the speficic force curve for variable width
    #  """
    #  self._Fs2 = list()
    #  for b in self._data['b']:
        #  self._Fs2.append(clib.specificForceCurve(self._yr, self._thetaS, self._data['beta'], self._data['Q'], self._g, b, self._data['theta1'], self._data['theta2']))
          
       
  #  def _specificEnergyCurve(self):
    #  """
    #  Estimate the speficic energy curve
    #  """

    #  self._Es = list()
    #  for Q in self._data['Q']:
        #  if self._data['ST'] == 1:
          #  self._Es.append(clib.specificEnergyCurveC(self._yr, self._thetaS, self._data['alpha'], Q, self._g, self._data['r']))
        #  elif self._data['ST'] == 2:
          #  self._Es.append(clib.specificEnergyCurve(self._yr, self._thetaS, self._data['alpha'], Q, self._g, self._data['b'][0], self._data['theta1'], self._data['theta2']))
        #  elif self._data['ST'] == 3:
          #  self._Es.append(clib.specificEnergyCurveI(self._yr, self._thetaS, self._data['alpha'], Q, self._g, self._data['xs'], self._data['ys']))
    #  #print(self._Es)

  #  def _specificEnergyCurve2(self):
    #  """
    #  Estimate the speficic energy curve for variable width
    #  """

    #  self._Es2 = list()
    #  for b in self._data['b']:
          #  self._Es2.append(clib.specificEnergyCurve(self._yr, self._thetaS, self._data['alpha'], self._data['Q'], self._g, b, self._data['theta1'], self._data['theta2']))

