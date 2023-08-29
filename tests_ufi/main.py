#!/usr/bin/python3

# -*- coding: utf-8 -*-

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from chydraulics.uficlass import *

def main():

  # Class to execute de calculus for unifor flow
  UniformFlowI()


if __name__ == '__main__':
  
  main()  


