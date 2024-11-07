#!/usr/bin/python3

# -*- coding: utf-8 -*-

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from chydraulics.abcclass import *

def main():

  # Class to compute the energy coefficient (alpha) and the momentum coefficient (beta) in channels.
  abCoeff()


if __name__ == '__main__':
  
  main()  


