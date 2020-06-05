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
import requests


# Basis set exchange interface with BEATLES
# Uses requests package (make sure to have it installed!)

print "_____________________________________________________________"
print "Get Basis Interface of the BEATLES v 0.1"
print "Retrieves basis set information from www.BasisSetExchange.org"
print "\n A New Basis Set Exchange: An Open, Up-to-date Resource for the Molecular Sciences Community.\n Benjamin P. Pritchard, Doaa Altarawy, Brett Didier, Tara D. Gibson, Theresa L. Windus.\n J. Chem. Inf. Model. 2019, 59(11), 4814-4820"
print "_____________________________________________________________"
print "This is beta version."
print "it can interface with www.basissetexchange.org but can't figure out what to do if basis set is misspelled/missing"
print "Last updated by Hassan Harb, June 1, 2020"
print "NOTICE: Remember to add a backslash if the basisset has a special character, e.g. 6-311++G* should be called as 6-311++G\*"
print"\n"

arguments_length = len(sys.argv)
arguments = []*(arguments_length - 2)

basis = sys.argv[1]
formats = sys.argv[2]
elements = '?elements='
print "Number of elements", arguments_length -3
element_list = [" "]*(arguments_length-3)

filename = basis+'_'

for i in range (3,arguments_length):
   if (i < arguments_length-1):
         elements = elements+sys.argv[i]+','
         filename = filename+sys.argv[i]+'-'
         element_list[i-3] = sys.argv[i]
   else:
         elements = elements+sys.argv[i]
         filename = filename+sys.argv[i]
         element_list[i-3] = sys.argv[i]

print "elements =", element_list
print "Basis =", basis
print "Format =", formats

if (formats == 'Gaussian' or formats == 'GAUSSIAN' or formats == 'gaussian'):
    formats = 'gaussian94'

Base_URL = 'https://www.basissetexchange.org/api'

full_path = Base_URL+'/basis/'+basis+'/format/'+formats+'/'+elements +'&make_general=true'
r2 = requests.get(full_path)
basisset = r2.text

if (r2.status_code != 200):
    print "ERROR, Basis set not found! Make sure you typed in the correct input."
    exit()

filename = filename.replace(' ','-')
filename = filename+'.bsf'

with open (filename, 'w') as bsf:
    bsf.write(basisset)

print "Basis set written to file ", filename

PrintLyrics()

