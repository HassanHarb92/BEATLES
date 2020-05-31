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


# Distances Calculator
#
# Script that calculates the distances between all atoms in a molecule and outputs them into a new text file
# 

print "-------------------------------------"
print "Internuclear Distance Calculator Code"
print "-------------------------------------\n"
filename1 = sys.argv[1]
filename3 = os.path.splitext(filename1)[0]+"-Distance.txt"

print "Retrieving Molecule's info from: ", filename1

NBasis, NElem, Charge, Multiplicity, NAtoms, SCFEnergy = NBasGrab(filename1)
Distance_Matrix, Atomic_Symbol = DistanceMatrix(filename1)

print "Internuclear Distances:\n"

for i in range(0,NAtoms-1):
    print "Atom "+Atomic_Symbol[i]+" ("+str(i+1)+") :"
    for j in range(i+1,NAtoms):
        print "Dist. "+Atomic_Symbol[i]+" ("+str(i+1)+") - "+Atomic_Symbol[j]+"("+str(j+1)+") = "+str(Distance_Matrix[i,j])
    print "\n"

with open(filename3, 'w') as output:
    output.write("Internuclear distances for atoms in "+ filename1+"\n")
    output.write("Number of Atoms = "+str(NAtoms)+"\n\n")
    for i in range(0,NAtoms-1):
      output.write("Atom "+Atomic_Symbol[i]+" ("+str(i+1)+") :\n")
      for j in range(i+1,NAtoms):
          output.write("Dist. "+Atomic_Symbol[i]+" ("+str(i+1)+") - "+Atomic_Symbol[j]+"("+str(j+1)+") = "+str(Distance_Matrix[i,j])+"\n")
      output.write("\n")

print "Internuclear distances successfully written to ", filename3

PrintLyrics()

