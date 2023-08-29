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

#def triangleArea(b, h):
#  """
#  Estimating Manning roughness factor for a circular channel using Manning equation
#  Where:
#    b = Triangle base
#    h = Triangle height
#  """
#  try:
#    return b*h*0.5
#  except ValueError:
#    print("Oops!  That was no valid number.  Try again...")
#
#def trapezoidArea(b1, b2, h):
#  """
#  Estimating Manning roughness factor for a circular channel using Manning equation
#  Where:
#    b1 = Trapezoid base 1
#    b2 = Trapezoid base 2
#    h =  Trapezoid height
#  """
#  try:
#    return (b1+b2)*h*0.5
#  except ValueError:
#    print("Oops!  That was no valid number.  Try again...")

def linearInterp(x1, y1, x2, y2, y):
  """
  Linear interpolation between pair of points. Return the x coord for y
  Where:
    x1 = x coord for point 1
    y1 = y coord for point 1
    x2 = x coord for point 2
    y2 = y coord for point 2
    y  = y coord. It is between y1 and y2.
  """
  dx = x2-x1
  try:
    if y == y1:
      return x1
    elif y == y2:
      return x2
    elif y<y1 and y>y2:
      return (y1-y)*dx/(y1-y2) + x1    
    elif y<y2 and y>y1:
      return -(y2-y)*dx/(y2-y1) + x2   
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def interpAnewSection(xs, ys, ns, y):
  """
  Get a new section given a y coord. Return new xs and ys
  Where:
    xs = List of x coords
    ys = List of y coords
    y = Given y coord.
  """
  try:
    n = len(xs)
    ymax = max(ys)
    xsn = xs[:]
    ysn = ys[:]
    nsn = ns[:]
    if y==ymax:
      return(xsn,ysn,nsn)
    elif y>ymax:
      xsn.insert(0,xs[0])
      ysn.insert(0,y)
      nsn.insert(0,ns[0])
      xsn.append(xs[n-1])
      ysn.append(y)
      nsn.append(ns[n-2])
      return(xsn,ysn,nsn)
    elif y<ymax:
      xsn = []
      ysn = []
      nsn = []
      sw = False
      for i in range(n-1):
        if ys[i]>y and ys[i+1]<y:
          sw = True
          xi = linearInterp(xs[i], ys[i], xs[i+1], ys[i+1], y)
          xsn.append(xi)
          ysn.append(y)
        elif ys[i+1]>y and ys[i]<y:
          xi = linearInterp(xs[i], ys[i], xs[i+1], ys[i+1], y)
          xsn.append(xi)
          ysn.append(y)
          nsn.append(ns[i])
          return(xsn,ysn,nsn)
        if sw:
          xsn.append(xs[i+1])
          ysn.append(ys[i+1])
          nsn.append(ns[i+1])
          
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def irregCrossSecA(xs, ys):
  """
  Area of a irregular channel cross section
  Where:
    xs = List of x coords
    ys = List of y coords
  """
  
#  try:
#    n = len(xs)
#    areaL = []
#    area_l = (xs[1]-xs[0])*(ys[0]-ys[1])*0.5
#    areaL.append(area_l)
#    area = 0
#    for i in range(1,n-2,1):
#      ai = ((ys[0]-ys[i])+(ys[0]-ys[i+1]))*(xs[i+1]-xs[i])*0.5
#      area += ai
#      areaL.append(ai)
#    area_r = (xs[n-1]-xs[n-2])*(ys[n-1]-ys[n-2])*0.5
#    areaL.append(area_r)
#    areaT = area_l + area_r + area
#    return(areaT, areaL)

  try:
    n = len(xs)
    areaL = []
    areaT = 0
    for i in range(n-1):
      ai = ((ys[0]-ys[i])+(ys[0]-ys[i+1]))*(xs[i+1]-xs[i])*0.5
      areaT += ai
      areaL.append(ai)
    return(areaT, areaL)

  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def irregCrossSecP(xs, ys):
  """
  Perimeter of a irregular channel cross section
  Where:
    xs = List of x coords
    ys = List of y coords
  """
  
  try:
    n = len(xs)
    peri = 0.
    periL = []
    for i in range(n-1):
      pi = math.sqrt( ((xs[i]-xs[i+1])**2.) + ((ys[i]-ys[i+1])**2.) )  
      peri += pi
      periL.append(pi)
    return(peri, periL)

  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def QmanningI(alpha, So, ns, xs, ys):
  """
  Estimation of discharge for an irregular cross section
  Where:
    alpha = Unit convertion factor    
    ns = List of Manning roughness coefficients
    So = Channel slope  
    xs = List of x coords
    ys = List of y coords
    """
  try:
    At,Al = irregCrossSecA(xs, ys)
    Pt,Pl = irregCrossSecP(xs, ys)
    n = len(Al)
    nt = 0.
    for i in range(n):
      nt += ns[i]*Al[i]
    nt = nt/At
    Rh = At/Pt
    return At*(alpha/nt)*(Rh**(2./3.))*(So**0.5)
    #print(Al, At)
    #print(Pl, Pt)
    #cte = alpha*(So**0.5)
    #Q = [0]*n
    #for i in range(n):
    #  #Q += (Al[i]**(5./3.))*(Pl[i]**(-2./3.))*(1./ns[i])
    #  Rh = Al[i]/Pl[i]
    #  Q[i]= Al[i]*(alpha/ns[i])*(Rh**(2./3.))*(So**0.5)
    #return sum(Q)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")


