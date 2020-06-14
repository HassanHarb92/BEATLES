#!/usr/bin/python

from __future__ import division
import sys
import math
import cmath
import numpy as np
from numpy import genfromtxt
import csv
from decimal import Decimal 
from BEATLES import *

filename1 = sys.argv[1]

NBasis, NElem, Charge, Multiplicity, NAtoms, SCFEnergy = NBasGrab(filename1)
NElec_1, NAlpha_1, NBeta_1 = NElec(filename1)

print "------------------------------------- "
print "Corresponding Orbitals"
print "version 1.1"
print "Hassan Harb & Hrant P. Hratchian"
print "University of California, Merced"
print "Last edited by Hassan Harb, May 29, 2020"
print "------------------------------------- \n"

print "Notice: This version works for both closed shell and open shell states"
print "        The current version only works over occupied orbitals"
print "Number of Basis functions =", NBasis
print "Number of Alpha electrons =", NAlpha_1
print "Number of Beta  electrons =", NBeta_1

#if (NAlpha_1 != NBeta_1):
#   print "Alpha and Beta electrons are not equal! Code will terminate."
#   exit()

Ca_1 = MatGrab(filename1,NBasis,1)
Cb_1 = MatGrab(filename1,NBasis,-1)

S = GetOverlap(Ca_1,NBasis)

Ca_1 = column2square(Ca_1,NBasis)
Cb_1 = column2square(Cb_1,NBasis)

Ca_1_O, Ca_1_V = OVParse(Ca_1,NBasis,NAlpha_1)
Cb_1_O, Cb_1_V  = OVParse(Cb_1,NBasis,NBeta_1)

Sigma_OO, U_OO, V_OO, D_OO = Biorthog(Ca_1_O,Cb_1_O,S,-1)
Sigma_VV, U_VV, V_VV, D_VV = Biorthog(Ca_1_V,Cb_1_V,S,-1)

Ca_1_O_Biorth = np.dot(Ca_1_O,np.transpose(V_OO))*np.linalg.det(np.transpose(V_OO))
Ca_1_V_Biorth = np.dot(Ca_1_V,np.transpose(V_VV))*np.linalg.det(np.transpose(V_VV))

Cb_1_O_Biorth = np.dot(Cb_1_O,U_OO)*np.linalg.det(np.transpose(U_VV))
Cb_1_V_Biorth = np.dot(Cb_1_V,U_VV)*np.linalg.det(np.transpose(U_VV))

Ca_Biorth = OVMerge(Ca_1_O_Biorth,Ca_1_V_Biorth,NAlpha_1,NBasis)
Cb_Biorth = OVMerge(Cb_1_O_Biorth,Cb_1_V_Biorth,NBeta_1,NBasis)

if (NAlpha_1 != NBeta_1):
   Sigma_alpha_nulls = np.zeros(NAlpha_1-NBeta_1)
   Sigma_alpha = np.concatenate((Sigma_OO,Sigma_alpha_nulls))
   Sigma_alpha = np.concatenate((Sigma_alpha,Sigma_VV))
   Sigma_beta = np.concatenate((Sigma_OO,Sigma_alpha_nulls))
   Sigma_beta = np.concatenate((Sigma_beta,Sigma_VV))

else:
   Sigma_alpha = np.concatenate((Sigma_OO,Sigma_VV))
   Sigma_beta = np.concatenate((Sigma_OO,Sigma_VV))

array_alpha = np.zeros(NBasis)
array_beta = np.zeros(NBasis)
thresh_eigval = 0.01

for i in range (0,NAlpha_1):
    if (Sigma_alpha[i] < thresh_eigval):
       array_alpha[i] = 1
       print "Alpha orbital ", i+1, " has an overlap less than the threshold"
for i in range (0,NBeta_1):
    if (Sigma_beta[i] < thresh_eigval):
       array_beta[i] = 1
       print "Beta orbital ", i+1, " has an overlap less than the threshold"

array_alpha_sum = np.sum(array_alpha)
array_beta_sum  = np.sum(array_beta)

Ca_nulls = np.zeros((NBasis,array_alpha_sum))
Cb_nulls = np.zeros((NBasis,array_beta_sum))
j = 0
for i in range(0,NBasis):
    if (array_alpha[i] == 1):
       Ca_nulls[:,j] = Ca_Biorth[:,i]
       j = j + 1

k = 0
for i in range(0,NBasis):
    if (array_beta[i] ==1):
       Cb_nulls[:,k] = Cb_Biorth[:,i]
       k = k + 1

Sigma_alpha = np.sqrt(Sigma_alpha)
Sigma_beta = np.sqrt(Sigma_beta)

Ca_Biorth_final = Ca_1
Cb_Biorth_final = Cb_1

sigma_alpha_final = Sigma_alpha
sigma_beta_final = Sigma_beta

for i in range(0,NAlpha_1):
          Ca_Biorth_final[:,i] = Ca_Biorth[:,i]

for i in range(0,NBeta_1):
          Cb_Biorth_final[:,i] = Cb_Biorth[:,i]


filename3 = "CorrespondingOrbitals-"+filename1

WriteMOsQChem(filename1,filename3,Ca_Biorth,Cb_Biorth,sigma_alpha_final,sigma_beta_final,NBasis)

listoftitles = ["Orbital Number", "Occupancy", "Occ. abOV"]

print "\nLegend:"
print ('{:<15s}{:<20s}'.format("Docc.", "Doubly Occupied"))
print ('{:<15s}{:<20s}'.format("Asing.", "Alpha Singly Occupied"))
print ('{:<15s}{:<20s}'.format("Bsing.", "Beta Singly Occupied"))
print ('{:<15s}{:<20s}'.format("abOV", "Alpha - Beta Overlap"))
print "\n"
print "-------------------------------------"
print "Final Corresponding orbitals results"
print "-------------------------------------"
print "Alpha Corresponding orbitals"
for i in range(0,NAlpha_1+1):
    if (i==0):
       print('{:<20s}{:^20s}{:^20s}'.format(listoftitles[0],listoftitles[1],listoftitles[2]))
    else:
     if (Sigma_alpha[i-1] < 0.1):
         print('{:<20d}{:^20s}{:^20.2f}'.format(i,"Asing.",Sigma_alpha[i-1]))
     else:
         print('{:<20d}{:^20s}{:^20.2f}'.format(i,"Docc.",Sigma_alpha[i-1]))

print "\nBeta Corresponding orbitals"
for i in range(0,NBeta_1+1):
    if (i==0):
       print('{:<20s}{:^20s}{:^20s}'.format(listoftitles[0],listoftitles[1],listoftitles[2]))
    else:
     if (Sigma_beta[i-1] < 0.1):
         print('{:<20d}{:^20s}{:^20.2f}'.format(i,"Bsing.",Sigma_beta[i-1]))
     else:
         print('{:<20d}{:^20s}{:^20.2f}'.format(i,"Docc.",Sigma_beta[i-1]))

print "\nCorresponding orbital analysis done! To visualize the orbitals please open", filename3,"\n"
