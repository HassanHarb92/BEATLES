#!/usr/bin/python

from __future__ import division
import sys
import math
import cmath
import numpy as np
from numpy import genfromtxt
import csv
from decimal import Decimal
import os
import random
from BEATLES import *


# Natural Orbitals Transformation
#
# Script that reads in the canonical orbitals from a formatted checkpoint file and performs the Lowdin Natural orbital transformtion
# outputs a new fchk file that can be opened to visualize the natural orbitals
# 

filename1 = sys.argv[1]
filename3 = "NaturalOrbitals-"+filename1

NBasis, NElem, Charge, Multiplicity, NAtoms, SCFEnergy = NBasGrab(filename1)
NOvecsA, NOvecsB, NOvalsA, NOvalsB = CalcNO(filename1,NBasis)
WriteMOs(filename1,filename3,NOvecsA,NOvecsB,NOvalsA,NOvalsB,NBasis) 



