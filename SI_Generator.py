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

arguments = len(sys.argv)
SI_file = sys.argv[1]
SI_array = list() 

print "------------------------------------- "
print "Ringo (SI)-Star Supporting Information Generator"
print "Hassan Harb & Hrant P. Hratchian"
print "University of California, Merced"
print "To cite this program, please use the following DOI: "
print "------------------------------------- \n"

print "Generating SI file for ", arguments-1, " files"
print "SI data will be printed to ", SI_file

for i in range (2,arguments):
    filename = PrintSI(sys.argv[i],1)
    SI_array.append(filename)

print "List of SI files =", SI_array

with open(SI_file, 'w') as output:
     for i in range(0,len(SI_array)):
         with open (SI_array[i]) as inputfile:
              output.write(inputfile.read())
         os.remove(SI_array[i])



