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
from BEATLES import *


# MixandMatch:  Reads in two fchk files, generates two new fchk files new MO coefficients as follows
#               fchk1: Alpha_1, Alpha_2
#               fchk2: Beta_1, Beta_2


print "Mix and Match: Reads in two checkpoint files (chkpt file 1, chkpt file 2)."
print "               Generates two checkpoint files (chkpt file 3, chkpt file 4)as follows:\n"
print "Checkpoint file 3: New Alpha MOs = Alpha MOs from chkpt file 1. New Beta MOs = Alpha MOs from chkpt file 2."
print "Checkpoint file 4: New Alpha MOs = Beta MOs from chkpt file 1. New Beta MOs = Beta MOs from chkpt file 2."

filename1 = sys.argv[1]
filename2 = sys.argv[2]

NBasis, NElem, Charge, Multiplicity, NAtoms, SCFEnergy = NBasGrab(filename1)

C1_a = MatGrab(filename1,NBasis,1)
C1_b = MatGrab(filename1,NBasis,-1)
C2_a = MatGrab(filename2,NBasis,1)
C2_b = MatGrab(filename2,NBasis,-1)

C1_a = column2square(C1_a,NBasis)
C1_b = column2square(C1_b,NBasis)
C2_a = column2square(C2_a,NBasis)
C2_b = column2square(C2_b,NBasis)

AOE_a_1 = MatGrab(filename1,NBasis,3)
AOE_b_1 = MatGrab(filename1,NBasis,-3)
AOE_a_2 = MatGrab(filename2,NBasis,3)
AOE_b_2 = MatGrab(filename2,NBasis,-3)

# Write Alpha - Alpha Checkpoint file

filename3 = 'Alpha-Alpha-'+filename1
WriteMOs(filename1,filename3,C1_a,C2_a,AOE_a_1,AOE_a_2,NBasis)

# Write Beta - Beta Checkpoint file

filename4 = 'Beta-Beta-'+filename1
WriteMOs(filename1,filename4,C1_b,C2_b,AOE_b_1,AOE_b_2,NBasis)

PrintLyrics()

