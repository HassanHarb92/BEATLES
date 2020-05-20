#!/usr/bin/python

#######################################################################################################################
#
#  NIOs  v. 1.0
#  By Hassan Harb
#  Program that reads in SCF Densities of ground and detached states and returns the eigenvalues and eigenvectors of
#  natural ionization orbitals
#
#  usage: python nio_beatles.py ground.fchk detached.fchk
#
#  For more details, check: J. Chem. Phys., 144, 204117 (2016)
#
#  Last edited by Hassan Harb, June 18, 2019
#
#######################################################################################################################

from __future__ import division
import sys
import math
import cmath
import numpy as np
from numpy import genfromtxt
import csv
from decimal import Decimal 
from BEATLES import *

# Function: Calculate Ncol, the column vector that shows the MO contributions of NIOs

def CalNCol(T,k,NBasis,S):
   Ti = np.zeros(NBasis)
   for i in range(0,NBasis):
      Ti[i] = T[i,k].real

   T_in = np.dot(Ti,Ti)
   N = np.outer(Ti,Ti)
   I = np.identity(NBasis)
   N  = np.multiply(N,I)
   N_cont = np.sum(N)
   print "Inner product =", T_in
   print "Contraction = ", N_cont
   print "trace N = ", np.trace(N)

   NCol = np.zeros(NBasis)
   for i in range(0,NBasis):
       for j in range(0,NBasis):
          NCol[i] = NCol[i] + N[j,i]
       if (NCol[i] < 0.0000000001 ):
          NCol[i] = 0
   NIOPops(NCol,NBasis)
   return NCol

#Function: Calculate the percentage contributions of MOs in each NIO

def NIOPops(NCol,NBasis):
    for i in range(0,NBasis):
       if (np.absolute(NCol[i]) >= 0.01):
           print "The Contribution from MO ", i+1 , " is ", np.around(NCol[i],decimals=2)*100, "percent \n"
    print "---------------\n"

# Funciton: Calculate the Occ and Virt percentages in each NCol

def NColPop(NCol,NOcc):
    Occ = 0.0
    Virt = 0.0
    for i in range(0,NOcc):
        Occ = Occ + NCol[i]
    for j in range(NOcc,NBasis):
        Virt = Virt + NCol[j]

    return Occ, Virt

#Part 1: Read in the matrix files from two checkpoint files

filename1 = sys.argv[1]
filename2 = sys.argv[2]
filename3 = "NIO-"+filename1
acc = 8 #accuracy to the acc's decimal place
minEigVal = 0.15

print "------------------------------------- "
print "Natural Ionization Orbitals"
print "Hassan Harb & Hrant P. Hratchian"
print "University of California, Merced"
print "Last updated by Hassan Harb, May 19, 2020"
print "------------------------------------- \n"

print "NIO: Calculates the density difference between two checkpoint files.\n"
print "Checkpoint 1:", filename1
print "Checkpoint 2:", filename2

NBasis, NElem, Charge, Multiplicity, NAtoms, SCFEnergy = NBasGrab(filename1)

P1_a, P1_b = MatGrab(filename1,NBasis,2)
P2_a, P2_b = MatGrab(filename2,NBasis,2)

Ca_1 = MatGrab(filename1,NBasis,1)
Cb_1 = MatGrab(filename1,NBasis,-1)

S = GetOverlap(Ca_1,NBasis)

Svals, Svecs = np.linalg.eig(S)
Sval_minhalf = (np.diag(Svals**(0.5)))
Shalf = np.dot(Svecs,np.dot(Sval_minhalf,np.transpose(Svecs)))
Shalf = Shalf.real

DP_a = P2_a - P1_a
DP_b = P2_b - P1_b

DP_a = np.dot(Shalf,np.dot(DP_a,Shalf))
DP_b = np.dot(Shalf,np.dot(DP_b,Shalf))

e1, U1 = np.linalg.eig(DP_a)
e2, U2 = np.linalg.eig(DP_b)

e1 = e1.real
e2 = e2.real
U1 = U1.real
U2 = U2.real

V1 = np.dot(np.linalg.inv(Shalf),U1)
V2 = np.dot(np.linalg.inv(Shalf),U2)

WriteMOs(filename1,filename3,V1,V2,e1,e2,NBasis)

Ca_1 = column2square(Ca_1,NBasis)
Cb_1 = column2square(Cb_1,NBasis)

NElec, NAlpha, NBeta = NElec(filename1)

Ca_1_O, Ca_1_V = OVParse(Ca_1,NBasis,NAlpha)
Cb_1_O, Cb_1_V = OVParse(Cb_1,NBasis,NBeta)

Ta = np.dot(np.transpose(Ca_1),np.dot(S,V1))
Tb = np.dot(np.transpose(Cb_1),np.dot(S,V2))

eigs_alpha = 0
eigs_beta  = 0

for m in range(0,NBasis):
    if (np.absolute(e1[m].real) >= minEigVal):
        eigs_alpha = eigs_alpha + 1
    if (np.absolute(e2[m].real) >= minEigVal):
        eigs_beta = eigs_beta + 1

eigs_alpha_pop = np.zeros((eigs_alpha,3))
eigs_beta_pop = np.zeros((eigs_beta,3))
l = 0
s = 0
for i in range(0,NBasis):
    if (np.absolute(e1[i].real) >= minEigVal):
        print "The ", i,"th element of e1 has an eigenvalue of ", np.around(e1[i].real,decimals=3) ," perform population analysis\n"
        eigs_alpha_pop[l,0] = np.around(e1[i].real,decimals=3)
        print "Population analysis on Alpha NIO ", i+1, "\n"
        print "----------------------------\n"
        NCol_i = CalNCol(Ta,i,NBasis,S)
        eigs_alpha_pop[l,1], eigs_alpha_pop[l,2] = NColPop(NCol_i,NAlpha)
        l = l + 1

    if (np.absolute(e2[i].real) >= minEigVal):
        print "The ", i,"th element of e2 has an eigenvalue of ", np.around(e2[i].real,decimals=3) ," perform population analysis\n"
        eigs_beta_pop[s,0] = np.around(e2[i].real,decimals=3)
        print "Population analysis on Beta NIO ", i+1, "\n"
        print "----------------------------\n"
        NCol_i = CalNCol(Tb,i,NBasis,S)
        eigs_beta_pop[s,1], eigs_beta_pop[s,2] = NColPop(NCol_i,NBeta)
        s = s + 1

print "Number of alpha electrons =", NAlpha
print "Number of  beta electrons =",  NBeta
print "Number of basis functions =", NBasis
listoftitles = ["NIO Eigenvalue", "Percent Occ. (%)", "Percent Vir. (%)"]

print "-------------------"
print "Final NIO Results"
print "-------------------"
print "\nAlpha NIOs:\n"
for i in range(0,eigs_alpha+1):
    if (i==0): 
        print('{:<20s}{:^20s}{:^20s}{:^20s}'.format("NIO Number",listoftitles[0],listoftitles[1],listoftitles[2]))
    else:
       j = i - 1
       print('{:<20d}{:^20.2f}{:^20.2f}{:^20.2f}'.format(i,eigs_alpha_pop[j,0],eigs_alpha_pop[j,1]*100,eigs_alpha_pop[j,2]*100)) 
      
print "\nBeta NIOs:\n"
for i in range(0,eigs_beta+1):
    if (i==0):
       print('{:<20s}{:^20s}{:^20s}{:^20s}'.format("NIO Number",listoftitles[0],listoftitles[1],listoftitles[2]))
    else:
       j = i - 1
       print('{:<20d}{:^20.2f}{:^20.2f}{:^20.2f}'.format(i,eigs_beta_pop[j,0],eigs_beta_pop[j,1]*100,eigs_beta_pop[j,2]*100))
print "\n"
print "NIO Analysis done! To visualize the NIOs please open:", filename3

