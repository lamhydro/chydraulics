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
import matplotlib.pyplot as plt
import numpy as np

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
    theta1 = Right diagonal angle with the horizontal in radians 
    theta2 = Left diagonal angle with the horizontal in radians 
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
    theta1 = Right diagonal angle with the horizontal in radians
    theta2 = Left diagonal angle with the horizontal in radians
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
    theta1 = Right diagonal angle with the horizontal in radians
    theta2 = Left diagonal angle with the horizontal in radians 
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
    theta1 = Right diagonal angle with the horizontal in radians
    theta2 = Left diagonal angle with the horizontal in radians
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
    theta1 = Right diagonal angle with the horizontal in radians
    theta2 = Left diagonal angle with the horizontal in radians
    y = water depth
  """
  try:
    A = area(b, theta1, theta2, y)
    P = perimeter(b, theta1, theta2, y)
    Rh = A/P
    return ((Q*n)/(alpha*A*(Rh**(2.0/3.0))))**2.
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def SomanningI(alpha, Q, ns, xs, ys):  
  """
  Estimating channel slope for an irregular cross section using Manning equation
  Where:
    alpha = Unit convertion factor    
    Q  = Discharge  
    ns = List of Manning roughness coefficients
    xs = List of x coords
    ys = List of y coords
  """
  try:
    At,Al = irregCrossSecA(xs, ys)
    Pt,Pl = irregCrossSecP(xs, ys)
    na = weightedAve(Al, ns)
    Rh = At/Pt
    return ((Q*na)/(alpha*At*(Rh**(2.0/3.0))))**2.
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
  try:
    A = areaC(theta, r)
    P = perimeterC(theta, r)
    Rh = A/P
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
    theta1 = Right diagonal angle with the horizontal in radians
    theta2 = Left diagonal angle with the horizontal in radians
    y = water depth
  """
  try:
    A = area(b, theta1, theta2, y)
    P = perimeter(b, theta1, theta2, y)
    Rh = A/P
    return (alpha/Q)*A*(Rh**(2.0/3.0))*(So**0.5)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def NmanningI(alpha, Q, So, xs, ys):  
  """
  Estimating Manning roughness factor for an irregular cross section using Manning equation
  Where:
    alpha = Unit convertion factor    
    Q  = Discharge  
    So = Channel slope
    xs = List of x coords
    ys = List of y coords
  """
  try:
    At,Al = irregCrossSecA(xs, ys)
    Pt,Pl = irregCrossSecP(xs, ys)
    Rh = At/Pt
    return (alpha/Q)*At*(Rh**(2.0/3.0))*(So**0.5)
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
    return 180
  elif y < r:
    return math.degrees(2*math.acos((r-y)/r))
  elif y > r:
    return 360-math.degrees(2*math.acos((y-r)/r))
      #return (180.0/math.pi)*math.asin((y-r)/r)

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
        if ys[i]>=y and ys[i+1]<=y:
          sw = True
          xi = linearInterp(xs[i], ys[i], xs[i+1], ys[i+1], y)
          xsn.append(xi)
          ysn.append(y)
        elif ys[i+1]>=y and ys[i]<=y:
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


def weightedAve(wl, xl):
  """
  Estimation of weighted average
  Where:
    wl = List of weights 
    xl = List of values to be averaged
  """

  try:
    n = len(xl)
    xa = 0.
    for i in range(n):
      xa += wl[i]*xl[i]
    return xa/sum(wl)

  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def equiNmanning_EinsteinHorton(Ps, ns):
  """
  Estimation of equivalent Manning n, using the Einstein-Horton equation.
  Input:
  - Ps: List of perimeters for each segment
  - ns: List of Mannign n for each segment
  Output:
  - Equivalent Manning n
  """
  try:
    n = len(Ps)
    psns = 0.
    for i in range(n):
      psns += Ps[i]*((ns[i])**(3./2))
    return (psns/sum(Ps))**(2./3)

  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def equiNmanning_EinsteinHorton(Ps, ns):
  """
  Estimation of equivalent Manning n, using the Einstein-Horton equation.
  Input:
  - Ps: List of perimeters for each segment
  - ns: List of Mannign n for each segment
  Output:
  - Equivalent Manning n
  """
  try:
    n = len(Ps)
    psns = 0.
    for i in range(n):
      psns += Ps[i]*((ns[i])**(3./2))
    return (psns/sum(Ps))**(2./3)

  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def QI(xs, ys, vs):
  """
  Estimation of discharge for an irregular cross section 
  Where:
    xs = List of x coords
    ys = List of y coords
    vs = List of velocities en the area subsections
  """

  # Get the areas in the subsections
  _,Al = irregCrossSecA(xs, ys)

  return (sum(v*a for v,a in zip(vs,Al)))


def QmanningI(alpha, So, ns, xs, ys):
  """
  Estimation of discharge for an irregular cross section from the Mannin equation
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
    #na = weightedAve(Al, ns)
    na = equiNmanning_EinsteinHorton(Pl,ns)
    Rh = At/Pt
    return At*(alpha/na)*(Rh**(2./3.))*(So**0.5)
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

def energyCoeff(xs, ys, vs):
  """
  Estimation of the correction coefficient for the kinematic energy
  Where:
    xs = List of x coords
    ys = List of y coords
    vs = List of velocities en the area subsections
  """

  try:
    # Get the areas in the subsections
    At,Al = irregCrossSecA(xs, ys)

    # Compute the coeff.
    num1 = 0
    den = 0
    for v,a in zip(vs,Al):
        num1 += (v**3)*a
        den += v*a
    return(num1*(At**2)/(den**3)) 
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def momentCoeff(xs, ys, vs):
  """
  Estimation of the correction coefficient for the momentum
  Where:
    xs = List of x coords
    ys = List of y coords
    vs = List of velocities en the area subsections
  """

  try:
    # Get the areas in the subsections
    At,Al = irregCrossSecA(xs, ys)

    # Compute the coeff.
    num1 = 0
    den = 0
    for v,a in zip(vs,Al):
        num1 += (v**2)*a
        den += v*a
    return(num1*At/(den**2)) 
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def slopeAngle(So):
  """
  Estimate the channel slope angle base on the slope
  Where:
    So = Channel slope (L/L)
  """

  try:
    return(math.degrees(math.atan(So))) 
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def specificEnergy(y, thetaS, alpha, Q, g, A):
  """
  Estimate the speficic energy
  Where:
    y = Vertical water depth
    thetaS = Slope angle channel
    alpha = Energy coefficient
    Q = Discharge
    g = Gravity
    A = Cross section area
  """

  try:
    thetar = math.radians(thetaS)
    E = y*(math.cos(thetar)) +alpha*(Q**2)/(2*g*(A**2)) 
    return(E)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def specificEnergyCurveI(yr, thetaS, alpha, Q, g, xs, ys):
  """
  Estimate the speficic energy curve for a irregular channel
  Where:
    yr = List of vertical water depths
    thetaS = Slope angle channel
    alpha = Energy coefficient
    Q = Discharge
    g = Gravity
    xs = List of x coords
    ys = List of y coords
  """

  try:
    Es = list()
    for y in yr:
      xsn, ysn, dummy = interpAnewSection(xs, ys, xs, y)
      A,_ = irregCrossSecA(xsn, ysn)
      Es.append(specificEnergy(y, thetaS, alpha, Q, g, A))
    return(Es)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def specificEnergyCurve(yr, thetaS, alpha, Q, g, b, theta1, theta2):
  """
  Estimate the speficic energy curve for a prismatic channel
  Where:
    yr = List of vertical water depths
    theta = Slope angle channel
    alpha = Energy coefficient
    Q = Discharge
    g = Gravity
    b = Channel width
    theta1 = Right diagonal angle with the horizontal in radians
    theta2 = Left diagonal angle with the horizontal in radians
  """

  try:
    Es = list()
    for y in yr:
      A = area(b, theta1, theta2, y)
      Es.append(specificEnergy(y, thetaS, alpha, Q, g, A))
    return(Es)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")


def specificEnergyCurveC(yr, thetaS, alpha, Q, g, r):
  """
  Estimate the speficic energy curve for a circular channel
  Where:
    yr = List of vertical water depths
    theta = Slope angle channel
    alpha = Energy coefficient
    Q = Discharge
    g = Gravity
    r = Radio
  """

  try:
    Es = list()
    for y in yr:
      theta = thetaInC(r, y)
      #print(y, theta)
      A = areaC(theta, r)
      Es.append(specificEnergy(y, thetaS, alpha, Q, g, A))
    return(Es)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def specificForce(beta, Q, g, A, z):
  """
  Estimate the speficic energy
  Where:
    beta = Momentum coefficient
    Q = Discharge
    g = Gravity
    A = Cross section area 
    z = vertical depth to the area centroid
  """

  try:
    F = A*z + beta*(Q**2)/(g*A) 
    return(F)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def specificForceCurve(yr, thetaS, beta, Q, g, b, theta1, theta2):
  """
  Estimate the speficic force curve for a prismatic channel
  Where:
    yr = List of vertical water depths
    theta = Slope angle channel
    beta = Momentum coefficient
    Q = Discharge
    g = Gravity
    b = Channel width
    theta1 = Right diagonal angle with the horizontal in radians
    theta2 = Left diagonal angle with the horizontal in radians
  """

  try:
    Fs = list()
    for y in yr:
      A = area(b, theta1, theta2, y)
      xs,ys = getXYcoordSec(b, theta1, theta2, y)
      _,z = polygon_centroid(xs,ys)
      Fs.append(specificForce(beta, Q, g, A, z))
    return(Fs)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")


def specificForceCurveI(yr, thetaS, beta, Q, g, xs, ys):
  """
  Estimate the speficic force curve for a irregular channel
  Where:
    yr = List of vertical water depths
    thetaS = Slope angle channel
    beta = Momentum coefficient
    Q = Discharge
    g = Gravity
    xs = List of x coords
    ys = List of y coords
  """

  try:
    Fs = list()
    for y in yr:
      xsn, ysn, dummy = interpAnewSection(xs, ys, xs, y)
      A,_ = irregCrossSecA(xsn, ysn)
      _,z = polygon_centroid(xsn,ysn)
      Fs.append(specificForce(beta, Q, g, A, z))
    return(Fs)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def specificForceCurveC(yr, thetaS, beta, Q, g, r):
  """
  Estimate the speficic force curve for a circular channel
  Where:
    yr = List of vertical water depths
    thetaS = Slope angle channel
    beta = Energy momentum 
    Q = Discharge
    g = Gravity
    r = Radio
  """

  try:
    Fs = list()
    for y in yr:
      theta = thetaInC(r, y)
      #print(y, theta)
      A = areaC(theta, r)
      z = cirseg_centroid(r,theta,y)
      Fs.append(specificForce(beta, Q, g, A, z))
    return(Fs)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")


def getXYcoordSec(b, theta1, theta2, y):
  """
  Set the x and y coordinates of a prismatic section
  where:
    b = Channel width
    theta1 = Right diagonal angle with the horizontal in radians
    theta2 = Left diagonal angle with the horizontal in radians
    y = water depth
  """
  try:
      
    if theta1 == math.pi/2:
        a1 = 0.0
    elif theta1 < math.pi/2:
        a1 = y/math.tan(theta1)

    if theta2 == math.pi/2:
        a2 = 0.0
    elif theta2 < math.pi/2:
        a2 = y/math.tan(theta2)

    xs = list()
    ys = list()
    if b==0:
      xs.append(0); ys.append(y)
      xs.append(a1); ys.append(0)
      xs.append(a1+a2); ys.append(y)
    else:
      xs.append(0); ys.append(y)
      xs.append(a1); ys.append(0)
      xs.append(a1+b); ys.append(0)
      xs.append(a1+b+a2); ys.append(y)
    return(xs,ys)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def cirseg_centroid(r,theta,y):
  """
  Estimate the centroid y coordinate of a circular segment
  where:
    r: Circle radio
    theta: Angle formed between y and r in degrees
    y: height of the circular segment
  """
  try:
    theta_r = math.radians(theta)
    z = (4*r*((math.sin(theta_r/2.))**3.)/(3*(theta_r-math.sin(theta_r))))-(r-y)
    #print(y,z)
    return(z)
  except ValueError:
    print("Oops!  That was no valid number.  Try again...")

def polygon_centroid(x,y):
  """
  Estimate the coord x and y of the centroid of a polygon 
  Where:
    xs: a list with x coord of the section
    ys: a list with y coord of the section
  """

  try:
 
    # Vertices should be a list of (x, y) tuples
    x = np.array(x)
    y = np.array(y)
    
    # Close the polygon by adding the first vertex to the end
    x = np.append(x, x[0])
    y = np.append(y, y[0])
    
    # Calculate the area
    A = 0.5 * np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
    
    # Calculate centroid coordinates
    Cx = (1 / (6 * A)) * np.sum((x[:-1] + x[1:]) * (x[:-1] * y[1:] - x[1:] * y[:-1]))
    Cy = (1 / (6 * A)) * np.sum((y[:-1] + y[1:]) * (x[:-1] * y[1:] - x[1:] * y[:-1]))
    
    return Cx, Cy

  except ValueError:
    print("Oops!  That was no valid number.  Try again...")
