#!/usr/bin/python3

# -*- coding: utf-8 -*-

"""
This is repository with functions for open channel flow hydraulics
"""

# Defining some constants
ERROR = 1.0e-5
NMAX = 10000
YI = 0.01
YF = 1000
TI = 0.01
TF = 179.9

import math

def gravity(US):
  """
  Return the gravity acceleration depend on the unit system
  Where:
    US: Is the unit system. US='IS' or US='BG'
  """
  try:
    if US=='IS':
      return 9.81
    elif US=='BG':
      return 32.2
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def convertionFactor(US):
  """
  Return the convertion factor in the Manning equation
  Where:
    US: Is the unit system. US='IS' or US='BG'
  """
  try:
    if US=='IS':
      return 1.0
    elif US=='BG':
      return 1.49
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def area(b, theta1, theta2, y):
  """
  Area of the cross section of a rectagular, triangular o trapezoidal channel.
  Where:
    b = Channel width
    theta1 = Right diagonal angle with the horizontal in degrees
    theta2 = Left diagonal angle with the horizontal in degrees
    y = water depth
  """
  try:
    return (y/2.0)*(y*( (math.cos(theta1)/math.sin(theta1))+  (math.cos(theta2)/math.sin(theta2)) ) + 2*b) 
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def perimeter(b, theta1, theta2, y):
  """
  Perimeter of the cross section of a rectagular, triangular o trapezoidal channel.
  Where:
    b = Channel width
    theta1 = Right diagonal angle with the horizontal in degrees
    theta2 = Left diagonal angle with the horizontal in degrees
    y = water depth
  """
  try:
    return y*( 1./math.sin(theta1) +  1./math.sin(theta2) ) + b 
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def hydraulicRatio(b, theta1, theta2, y):
  """
  Hydraulic ratio of the cross section of a rectagular, triangular o trapezoidal channel.
  Where:
    b = Channel width
    theta1 = Right diagonal angle with the horizontal in degrees
    theta2 = Left diagonal angle with the horizontal in degrees
    y = water depth
  """
  try:
    return area(b, theta1, theta2, y)/perimeter(b, theta1, theta2, y)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def areaC(theta, r):
  """
  Area of the cross section of a circular channel.
  Where:
    r = Radius of cross section
    theta = Angle formed between the vertical axe and the radio to y in degrees
  """
  try:
    theta_r = math.radians(theta)
    return  (r**2.)*(theta_r-(math.sin(2.0*theta_r)/2.0))
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def perimeterC(theta, r):
  """
  Perimeter of the cross section of a circular channel.
  Where:
    r = Radius of cross section
    theta = Angle formed between the vertical axe and the radio to y in degrees
  """
  try:
    theta_r = math.radians(theta)
    return 2*r*theta_r
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def Qmanning(alpha, So, n, b, theta1, theta2, y):
  """
  Estimation of discharge for a non circular channel using Manning equation
  Where:
    alpha = Unit convertion factor    
    So = Channel slope  
    n = Manning roughness factor
    b = Channel width
    theta1 = Right diagonal angle with the horizontal in degrees
    theta2 = Left diagonal angle with the horizontal in degrees
    y = water depth
  """
  A = area(b, theta1, theta2, y)
  P = perimeter(b, theta1, theta2, y)
  Rh = A/P
  try:
    return (alpha/n)*A*(Rh**(2.0/3.0))*(So**0.5)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def QmanningC(alpha, So, n, theta, r):
  """
  Estimation of discharge for a non circular channel using Manning equation
  Where:
    alpha = Unit convertion factor    
    So = Channel slope  
    n = Manning roughness factor
    theta = Angle formed between the vertical axe and the radio to y in degrees
    r = Radius of cross section
  """
  A = areaC(theta, r)
  P = perimeterC(theta, r)
  Rh = A/P
  try:
    return (alpha/n)*A*(Rh**(2.0/3.0))*(So**0.5)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def Somanning(alpha, Q, n, b, theta1, theta2, y):  
  """
  Estimating channel slope for a non circular channel using Manning equation
  Where:
    alpha = Unit convertion factor    
    Q = Discharge
    n = Manning roughness factor
    b = Channel width
    theta1 = Right diagonal angle with the horizontal in degrees
    theta2 = Left diagonal angle with the horizontal in degrees
    y = water depth
  """
  A = area(b, theta1, theta2, y)
  P = perimeter(b, theta1, theta2, y)
  Rh = A/P
  try:
    return ((Q*n)/(alpha*A*(Rh**(2.0/3.0))))**2.
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def SomanningC(alpha, Q, n, theta, r):  
  """
  Estimating channel slope for a circular channel using Manning equation
  Where:
    alpha = Unit convertion factor    
    Q = Discharge
    n = Manning roughness factor
    theta = Angle formed between the vertical axe and the radio to y in degrees
    r = Radius of cross section
  """
  A = areaC(theta, r)
  P = perimeterC(theta, r)
  Rh = A/P
  try:
    return ((Q*n)/(alpha*A*(Rh**(2.0/3.0))))**2.
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def Nmanning(alpha, Q, So, b, theta1, theta2, y):  
  """
  Estimating Manning roughness factor for a non circular channel using Manning equation
  Where:
    alpha = Unit convertion factor    
    Q = Discharge
    So = Channel slope
    b = Channel width
    theta1 = Right diagonal angle with the horizontal in degrees
    theta2 = Left diagonal angle with the horizontal in degrees
    y = water depth
  """
  A = area(b, theta1, theta2, y)
  P = perimeter(b, theta1, theta2, y)
  Rh = A/P
  try:
    return (alpha/Q)*A*(Rh**(2.0/3.0))*(So**0.5)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def NmanningC(alpha, Q, So, theta, r):  
  """
  Estimating Manning roughness factor for a circular channel using Manning equation
  Where:
    alpha = Unit convertion factor    
    Q = Discharge
    So = Channel slope
    theta = Angle formed between the vertical axe and the radio to y in degrees
    r = Radius of cross section
  """
  A = areaC(theta, r)
  P = perimeterC(theta, r)
  Rh = A/P
  try:
    return (alpha/Q)*A*(Rh**(2.0/3.0))*(So**0.5)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def thetaInC(r, y):
  """
  Estimating the angle for a circular channel given y
  Where:
    r = Radius of cross section
    y = water depth
  """
  if y == r:
    return 90
  elif y < r:
    return (180.0/math.pi)*math.acos((r-y)/r)
  elif y > r:
    return (180.0/math.pi)*math.asin((y-r)/r)


