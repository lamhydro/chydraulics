#!/usr/bin/python3

# -*- coding: utf-8 -*-

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from chydraulics.cfclass import *

def main():

  # Class to estimate the critical flow depth
  CriticalFlow()

if __name__ == '__main__':
  
  main()  


