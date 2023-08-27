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
      print('Normal discharge (Q) = %f8.3f ft³/s' % self._data['Q']) 


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

      if self._data['ST'] == 1:
        fa = self._data['Q'] - clib.QmanningC(self._conf, self._data['So'], self._data['n'], self._a, self._data['r'])   
        fb = self._data['Q'] - clib.QmanningC(self._conf, self._data['So'], self._data['n'], self._b, self._data['r'])   
        fc = self._data['Q'] - clib.QmanningC(self._conf, self._data['So'], self._data['n'], c, self._data['r'])   
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
        yn = self._data['r']*(1.0 + math.sin(math.radians(c)))
    elif self._data['ST'] == 2:
        yn = c
  
    # Printing results
    if self._data['US'] == 'IS':
      print('Normal depth (y_n) = %f8.3 m' % yn) 
    elif self._data['US'] == 'BG':
      print('Normal depth (y_n) = %f8.3 ft' % yn) 

class SimplePipes():

  def __init__(self):

    # Set the input data
    self._data = InputData()._data
    self._typec = InputData()._typec
    #self._Ht = InputData()._Ht
    self._E1 = InputData()._E1
    self._E2 = InputData()._E2
    self._SK = InputData()._SK
    self._g = InputData()._g

    # Executing de calculation
    self.procedure()

  def procedure(self):
    if self._typec ==1:
      print('')
      print('Proving the system design')
      print('')
      self.designTest()
    elif self._typec ==2:
      print('')
      print('Estimating system power')
      print('')
      self.systemPower()
    elif self._typec ==3:
      print('')
      print('System design')
      print('')
      self.pipeDesign()
      
  def designTest(self): 
    """
    Estimate de velocity/discharge
    """
    if self._data['IM'] == 'fp':
      self._res = splib.f_fp2(self._g, self._data['ks'], self._data['rho'], self._data['mu'], self._data['D'], self._E1, self._E2, self._data['Pu']['h'], self._data['Tu']['h'], self._data['L'], self._SK)
    elif self._data['IM'] == 'nr': 
      self._res = splib.f_nr2(self._g, self._data['ks'], self._data['rho'], self._data['mu'], self._data['D'], self._E1, self._E2, self._data['Pu']['h'], self._data['Tu']['h'], self._data['L'], self._SK)
   
    # Set variables 
    self._V = self._res['V']
    self._f = self._res['f']
    self._hf = self._res['hfrl'][-1]

    # Estimate accesory losses
    self._het = splib.he(self._g,self._SK,self._V)

    # Discharge estimation
    self._Q =  splib.Qc(self._V, self._data['D']) 

    # Get flow regime
    self._fr = splib.flowRegime(self._data['rho'], self._data['mu'], self._data['D'], self._V)

    # Printing results
    self.printIter_designTest()
    self.print_designTest_results()

  def printIter_designTest(self):
    """
    Print iteration table for design test
    """
    print('   Print iteration table for design test   ')
    #print(pd.DataFrame(self._table, columns = ["hf", "he", "f", "V", "Dhf"]))
    self._table = list(zip(self._res['fl'], self._res['Vl'], self._res['hfrl']))
    print(pd.DataFrame(self._table, columns = ["f","V","hf"]))

  def print_designTest_results(self):
    """
    Print design test results
    """
    if self._data['US']=='IS':
      print("Q = %8.4f m³/s" % self._Q)
      print("V = %8.4f m/s" % self._V)
      print("f = %8.4f" % self._f)
      print("hf = %8.4f m" % self._hf)
      print("he = %8.4f m" % self._het)
      print("Flow regime = %s" % self._fr)
    elif self._data['US']=='ES':
      print("Q = %8.4f ft³/s" % self._Q)
      print("V = %8.4f ft/s" % self._V)
      print("f = %8.4f" % self._f)
      print("hf = %8.4f ft" % self._hf)
      print("he = %8.4f ft" % self._het)
      print("Flow regime = %s" % self._fr)

   
  def systemPower(self):
    """
    Estimate the system power
    """
    
    # Calculate de velocity
    self._V = splib.Vc(self._data['Q'], self._data['D'])
    
    # Estimate he
    self._het = splib.he(self._g, self._SK, self._V)

    # Estimate f
    if self._data['IM'] == 'fp':
      self._fdic = splib.f_fp(self._g, self._data['ks'], self._data['rho'], self._data['mu'], self._data['D'], self._V)
    elif self._data['IM'] == 'nr': 
      self._fdic = splib.f_nr(self._g, self._data['ks'], self._data['rho'], self._data['mu'], self._data['D'], self._V)
    
    # Estimate hf
    self._hf = splib.hf(self._g, self._fdic['f'], self._data['L'], self._data['D'], self._V)

    # Total head energy deliver by a pumb from Bernoulli equation
    self._Ht = splib.HPu(self._E1, self._E2, self._het, self._hf, self._data['Tu']['h'])

    # Nominal system power
    self._P = splib.pot(self._data['US'], self._data['Pu']['ef'], self._data['Q'], self._Ht, self._g, self._data['rho'])
    
    # Get flow regime
    self._fr = splib.flowRegime(self._data['rho'], self._data['mu'], self._data['D'], self._V)

    # Printing results
    self.printIter_systemPower()
    self.print_systemPower_results()


  def printIter_systemPower(self):
    """
    Print iteration table for system power
    """
    print('   Print iteration table for system power   ')
    self._table = list(zip(self._fdic['f_list'], self._fdic['df']))
    print(pd.DataFrame(self._table, columns = ["f","df"]))

  def print_systemPower_results(self):
    """
    Print system power results
    """
    if self._data['US']=='IS':
      print("P = %8.2f W" % self._P)
      print("V = %8.4f m/s" % self._V)
      print("f = %8.6f" % self._fdic['f'])
      print("hf = %8.4f m" % self._hf)
      print("he = %8.4f m" % self._het)
      print("Hp = %8.4f m" % self._Ht)
      print("Flow regime = %s" % self._fr)
    elif self._data['US']=='ES':
      print("P = %8.4f Hp" % self._P)
      print("V = %8.4f ft/s" % self._V)
      print("f = %8.6f" % self._f)
      print("hf = %8.4f ft" % self._hf)
      print("he = %8.4f ft" % self._het)
      print("Hp = %8.4f ft" % self._Ht)
      print("Flow regime = %s" % self._fr)


  def pipeDesign(self):
    """
    Estimation of pipe diameter
    """
     
    # Comertial diameter list in inches
    #CD = [1./8, 1./4, 3./8, 1./2, 3./4, 1., 1.+1./4, 1.+1./2, 2., 2.+1./2, 3., 3.+1./2, 4, 4.+1./2, 5., 6., 8., 10., 12., 14., 16., 18., 20., 24., 28., 32., 36., 40., 42., 44., 48., 52., 56., 60.] 
    #CD = [80.42, 103.42, 152.22, 198.48, 247.09, 293.07]
    CD={'3/4':23.63,'1':30.20,'1 1/4':38.14,'1 1/2':43.68,'2':54.58,'2 1/12':66.07,'3':80.42,'4':103.42,'6':152.22,'8':198.21,'10':247.09,'12':293.07,'14':321.76,'16':367.70,'18':413.66,'20':459.64,'24':551.54}   

    # Loop to find the optimal comertial diameter
    self._table = []
    #for D in CD:
    for self._D, self._data['D'] in CD.items():
      self._data['D'] *= 0.001 
      if self._data['IM'] == 'fp':
        self._res = splib.f_fp2(self._g, self._data['ks'], self._data['rho'], self._data['mu'], self._data['D'], self._E1, self._E2, self._data['Pu']['h'], self._data['Tu']['h'], self._data['L'], self._SK)
      elif self._data['IM'] == 'nr': 
        self._res = splib.f_nr2(self._g, self._data['ks'], self._data['rho'], self._data['mu'], self._data['D'], self._E1, self._E2, self._data['Pu']['h'], self._data['Tu']['h'], self._data['L'], self._SK)
      # Set variables 
      self._V = self._res['V']
      self._f = self._res['f']
      self._hf = self._res['hfrl'][-1]

      # Estimate accesory losses
      self._het = splib.he(self._g,self._SK,self._V)

      # Discharge estimation
      self._Q =  splib.Qc(self._V, self._data['D']) 
      testQ = self._Q >= self._data['Q']
       
      self._table.append([self._data['D']*1000., self._Q, testQ, self._f, self._het, self._hf])
      #if self._Q >= self._data['Q']:
      if testQ:
        break

    # Get flow regime
    self._fr = splib.flowRegime(self._data['rho'], self._data['mu'], self._data['D'], self._V)

    # Printing results
    self.printIter_pipeDesign()
    self.print_pipeDesign_results()


  def printIter_pipeDesign(self):
    """
    Print iteration table for system power
    """
    print('   Print iteration table for pipe design   ')
    print(pd.DataFrame(self._table, columns = ["D(mm)", "Q", "Q>=Qd", "f", "he", "hf"]))

  def print_pipeDesign_results(self):
    """
    Print system power results
    """
    if self._data['US']=='IS':
      print("D = %s pulg" % self._D)
      print("Q = %8.4f m/s" % self._Q)
      print("f = %8.6f" % self._f)
      print("hf = %8.4f m" % self._hf)
      print("he = %8.4f m" % self._het)
      print("Flow regime = %s" % self._fr)
    elif self._data['US']=='ES':
      print("D = %s pulg" % self._DC)
      print("Q = %8.4f ft/s" % self._Q)
      print("f = %8.6f" % self._f)
      print("hf = %8.4f ft" % self._hf)
      print("he = %8.4f ft" % self._het)
      print("Flow regime = %s" % self._fr)













#  def designTest2(self):
#    """
#    Estimate de velocity/discharge
#    """
#    # Main loop
#    self._hf = self._Ht
#    #while abs(diff)<=CONST['error']:
#    self._table=[]
#    while True:
#      # Estimate Vi
#      self._V = vel(self._g, self._data['ks'], self._data['rho'], self._data['mu'], self._data['D'], self._data['L'], self._hf)
#    
#      # Friction factor
#      self._f = f_dw(self._g, self._hf, self._data['L'], self._data['D'], self._V)
#    
#      # Estimate hf2
#      hf2 =hfb(self._g, self._Ht, self._data['zo'], self._data['K'], self._V)
#    
#      # Estimate he
#      self._het = he(self._g,sum(self._data['K']),self._V)
#    
#      self._diff = hf2-self._hf
#      self._table.append([self._hf, self._het, self._f, self._V, self._diff])
#      if abs(self._diff)<=1.e-5:
#        break
#      self._hf = hf2
#    
#    # Discharge estimation
#    self._Q =  Qc(self._V, self._data['D']) 
#
#    # Printing results
#    self.printIter_designTest()
#    self.print_designTest_results()



    # Suppose diameter (must be commertial and small
#  def CD(i):
#  """
#  Return comertial diameter based on the i.
#  """
#  if i in range(4):
#    return i*1./8
#  elif i in range(4,8):
#    return i*2*1./8
#  elif i in range(9,15):
#    return i*4*1./8
#  elif i == 16:
#    return 6.
#  elif i in range(17,23):
#    return i*16*1./8
#  elif i in range(24,2):

#  def pipeDesign(self):
#    """
#    Estimate of the pipe comercial diameter
#    """
#    
#    # Suppose hf  
#    self._hf = self._Ht
#     
#    # Suppose diameter (must be commertial and small
#    CD = [1./8, 1./4, 3./8, 1./2, 3./4, 1., 1.+1./4, 1.+1./2, 2., 2.+1./2, 3., 3.+1./2, 4, 4.+1./2, 5., 6., 8., 10., 12., 14., 16., 18., 20., 24., 28., 32., 36., 40., 42., 44., 48., 52., 56., 60.] 
#    self._data['D'] = CD[0]
#
#    i = 1
#    while True:
#      # Estimate V
#      self._V = vel(self._g, self._data['ks'], self._data['rho'], self._data['mu'], self._data['D'], self._data['L'], self._hf)
#    
#      # Estimate Q
#      self._Q =  Qc(self._V, self._D) 
#
#      if self._Q >= self._data['Q']:
#        while True:
#          # Estimate hf2
#          hf2 =hfb(self._g, self._Ht, self._data['zo'], self._data['K'], self._V)
#
#          if abs(hf2-self._hf)<=1.e-5:
#            break
#          else:
#            # Estimate V
#            self._V = vel(self._g, self._data['ks'], self._data['rho'], self._data['mu'], self._data['D'], self._data['L'], self._hf)
#
#      else:
#        self._data['D'] = CD[i]
#        i+=1
#    
#    # Estimate Q
#    self._Q =  Qc(self._V, self._D) 
    
   

#if __name__ == '__main__':
  
