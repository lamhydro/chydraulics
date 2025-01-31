#!/usr/bin/python3

# -*- coding: utf-8 -*-

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from chydraulics.sefclass import *

def main():

  # Class to estimate the specific energy in a channel and to plot the curve (E vs y)
  SpecificEnergyForce()

if __name__ == '__main__':
  
  main()  


