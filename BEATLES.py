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

# BEATLES: Bundle of Essential and Assistive Tools Library for Electronic Structure

# NBasGrab: reads in a name of .fchk file
# output:  -Number of basis functions
#          -Charge
#          -Multiplicity
#          -Number of Atoms
#          -Cartesian Coordinates
#          -Atomic Symbols
#          -SCF Energy 
#          -Total Energy (needs to be added)

def NBasGrab(filename):
  NBasis = 0  
  NElem = 0
  Charge = 0
  Multiplicity = 0
  NAtoms = 0
  with open(filename, 'r') as origin:
     for line in origin:
        if "Number of basis functions" in line:
           words = line.split()
           for i in words:
              for letter in i:
                 if(letter.isdigit()):
                    NBasis = NBasis*10 + int(letter)
        if "Charge   " in line:
           words = line.split()
           for i in words:
              for letter in i:
                 if(letter.isdigit()):
                    Charge = Charge*10 + int(letter)
        if "Multiplicity" in line:
           words = line.split()
           for i in words:
              for letter in i:
                 if(letter.isdigit()):
                    Multiplicity = Multiplicity*10 + int(letter)
        if "Number of atoms" in line:
           words = line.split()
           for i in words:
              for letter in i:
                 if(letter.isdigit()):
                    NAtoms = NAtoms*10 + int(letter)

        if "SCF Energy" in line:
            words = line.split()
            print "SCF Energy = ", words[3], " Hartree"
            SCFEnergy = float(words[3])
            print "SCF Energy (float) = ", SCFEnergy

#        if "Total Energy" in line:
#            words = line.split()
#            TotalEnergy = float(words[3])
#            print "Total Energy = ", TotalEnergy, " Hartree"

  NElem = NBasis*NBasis
  print "Number of Basis Functions (subroutine) = ", NBasis, "\n"
  print "Charge (subroutine) = ", Charge, "\n"
  return NBasis, NElem, Charge, Multiplicity, NAtoms, SCFEnergy

# GeomGet: reads in the file name, number of atoms
# Output: -One dimensional vector (NAtoms * 3) that includes the cartesian coordinates of each atom
#
 
def GeomGet(filename,NAtoms):
   p = 0
   r = 0
   n = 1
   NElements = NAtoms * 3 
   RawCart = np.zeros(NElements)
   if (NElements%5 == 0):
      n = 0
   RawCartLines = int(NElements/5) + n
   print "Raw Cart lines = ", RawCartLines
   print "Number of Atoms =", NAtoms
   print "Number of coordinates =", NElements
   with open(filename,'r') as origin:
      for i, line in enumerate(origin):
         if "Current cartesian coordinates" in line:
            i = i + 1
            pointer = i
            print "Cartesian Coordinates starts at line :", pointer
            endpointer = pointer + RawCartLines - 1
            print "Cartesian Coordinates ends at line :", endpointer
            for m in range(0,endpointer - pointer +1):
               nextline = origin.next()
               nextline = nextline.split()
               for p in range(p,len(nextline)):
                  RawCart[r] = nextline[p]
                  r = r + 1
               p = 0
   print "Raw Cart (subroutine) = ", RawCart
   return  RawCart

# GetAtoms:  Reads in file name, number of atoms
# output:   -One dimensional vector (NAtoms) that contains the atomic numbers of the atoms 
#

def GetAtoms(filename1,NAtoms):
   p = 0
   r = 0
   n = 1
   AtomicNum = np.zeros(NAtoms)
   if (NAtoms%6 ==0):
     n = 0
   AtomLines = int(NAtoms/6) + n

   with open(filename1,'r') as origin:
      for i, line in enumerate(origin):
         if "Atomic numbers" in line:
            i = i + 1
            pointer = i
            endpointer = pointer + AtomLines -1
            for m in range(0, endpointer - pointer + 1):
                nextline = origin.next()
                nextline = nextline.split()
                for p in range(p,len(nextline)):
                   AtomicNum[r] = nextline[p]
                   r = r + 1
                p = 0
   return AtomicNum

# MatGrab: Reads in filename, NBasis, user-defined switch
# Output: -Alpha MO Coefficients (Done)
#         -Beta MO Coefficients (Done)
#         -Alpha Density Matrix (Done)
#         -Beta Density Matrix (Done)
#         -Alpha MO Energies (Done)
#         -Beta MO Energies (Done)
#
# Switch: 1 = Alpha MO Coefficients
#        -1 = Beta MO Coefficients
#         2 = Alpha and Beta Density Matrices
#         3 = Alpha MO Energies
#        -3 = Beta MO Energies
#
def MatGrab(filename,NBasis,switch):
   if (switch == 1):
      filename1 = filename
      MOElements = NBasis * NBasis
      MOlines = int(MOElements/5) + 1
      if (NBasis%5 == 0):
         MOlines = MOlines - 1
      p = 0
      r = 0
      AOE = 0
      MOrawa = np.zeros(NBasis*NBasis)
      with open(filename1,'r') as origin:
          for i, line  in enumerate(origin):
              if "Alpha Orbital Energies" in line:
                    AOE = i
              if  "Alpha MO coefficients" in line:
                    i=i+1
                    AMO=i
                    print "Alpha MO coefficients starts at line :", i
                    j=i+MOlines-1
                    print "Alpha MO coefficients ends at line :", j
                    for m in range(0,j-i+1):
                       nextline = origin.next()
                       nextline = nextline.split()
                       for p in range(p,len(nextline)):
                          MOrawa[r] = nextline[p]
                          r = r+1
                       p = 0
      print "MO Raw = ", MOrawa
      return MOrawa, AMO
   if (switch == -1):
      filename1 = filename
      MOElements = NBasis * NBasis
      MOlines = int(MOElements/5) + 1
      if (NBasis%5 == 0):
         MOlines = MOlines - 1
      p = 0
      r = 0
      BOE = 0
      MOrawb = np.zeros(NBasis*NBasis)
      with open(filename1,'r') as origin:
         for i, line  in enumerate(origin):
              if "Beta Orbital Energies" in line:
                    BOE = i
              if  "Beta MO coefficients" in line:
                    i=i+1
                    BMO=i
                    print "Beta MO coefficients starts at line :", i
                    j=i+MOlines-1
                    print "Beta MO coefficients ends at line :", j
                    for m in range(0,j-i+1):
                       nextline = origin.next()
                       nextline = nextline.split()
                       for p in range(p,len(nextline)):
                          MOrawb[r] = nextline[p]
                          r = r+1
                       p = 0

         print "MO Raw = ", MOrawb
         return MOrawb, BMO

   if (switch == 2):
      filename1 = filename
      PElements = int(NBasis*(NBasis+1)/2)
      Plines = int(PElements/5) + 1
      TotalPraw = np.zeros(PElements)
      SpinPraw = np.zeros(PElements)
      with open(filename1,'r') as origin:
       for i, line in enumerate(origin):
        if  "Total SCF Density" in line:
              i=i+1
              r = 0
              p = 0
              print "Total SCF Density starts at line :", i
              j=i+Plines-1
              print "Total SCF Density ends at line :", j
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   TotalPraw[r] = nextline[p]
                   r = r+1
                 p = 0
        if  "Spin SCF Density" in line:
              i=i+1
              r = 0
              p = 0
              print "Spin SCF Density starts at line: ", i
              j=i+Plines-1
              print "Spin SCF Density ends at line: ", j
              for m in range(0,j-i+1):
                 nextline = origin.next()
                 nextline = nextline.split()
                 for p in range(p,len(nextline)):
                   SpinPraw[r] = nextline[p]
                   r = r+1
                 p = 0
       PalphaRaw = (np.add(TotalPraw,SpinPraw)) * 0.5
       PbetaRaw = (np.subtract(TotalPraw,SpinPraw)) * 0.5
       Palpha = symmetrize(PalphaRaw)
       Pbeta  = symmetrize(PbetaRaw)       
       return Palpha, Pbeta

   if (switch == 3):
      filename1 = filename
      AlphaMO = np.zeros(NBasis)
      AlphaMOlines = int(NBasis/5) + 1
      if (NBasis % 5 == 0):
         AlphaMOlines = AlphaMOlines - 1
      with open(filename1,'r') as origin:
         for i, line in enumerate(origin):
             if "Alpha Orbital Energies" in line:
                i = i + 1
                r = 0
                p = 0
                print "Alpha MO Energies starts at line: ", i
                j = i + AlphaMOlines - 1
                print "Alpha MO Energies ends at line: ", j
                for m in range(0,j-i+1):
                    nextline = origin.next()
                    nextline = nextline.split()
                    for p in range(p,len(nextline)):
                        AlphaMO[r] = nextline[p]
                        r = r + 1
                    p = 0
      print "Alpha MO energies = ", AlphaMO
      return AlphaMO

   if (switch == -3):
      filename1 = filename
      BetaMO = np.zeros(NBasis)
      BetaMOlines = int(NBasis/5) + 1
      if (NBasis % 5 == 0):
         BetaMOlines = BetaMOlines - 1
      with open(filename1,'r') as origin:
         for i, line in enumerate(origin):
             if "Beta Orbital Energies" in line:
                i = i + 1
                r = 0
                p = 0
                print "Beta MO Energies starts at line: ", i
                j = i + BetaMOlines - 1
                print "Beta MO Energies ends at line: ", j
                for m in range(0,j-i+1):
                    nextline = origin.next()
                    nextline = nextline.split()
                    for p in range(p,len(nextline)):
                        BetaMO[r] = nextline[p]
                        r = r + 1
                    p = 0
      print "Beta MO energies = ", BetaMO
      return BetaMO

# sci_notation:  reads in a number
# output:        prints the number in the desired scientific notation. note that this function has a different output than the one found in nio.py
#
def sci_notation(n):
    a = '%.8f' % n
    return '%.8f' % Decimal(n.real)

# AtomicSymbol:  Reads in atomic number of the element
# Output:       -Atomic Symbol
# 

def AtomicSymbol(AtomicNumber):
    p = AtomicNumber - 1
    PTlist = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','T','Ru','Rh','Pd','Ah','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hb','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Uut','Fl','Uup','Lv','Uus','Uuo']
#    print "There are currently ", len(PTlist), " atoms defined"
    return PTlist[p]

# Symmetrize:  Reads in a packed symmetric column matrix into NBasis x NBasis square matrix 
# Output:     -Matrix(NBasis,NBasis)
#

def symmetrize(a):
  Nbas = int((np.sqrt(8*len(a)+1)-1)/2)
  b = np.zeros((Nbas,Nbas))
  n = 0
  for i in range(0,Nbas):
    for j in range(0,i+1):
      b[i,j]=a[n]
      b[j,i]=a[n]
      n=n+1
  return b

# Column2Square:  Reads in a packed column matrix, number of basis functions.
# Output:        -Matrix(NBasis,NBasis)

def column2square(A,NBasis):
  C = np.zeros((NBasis,NBasis))
  t=0
  for i in range(0,NBasis):
     for j in range(0,NBasis):
       C[j,i]=A[t]
       t=t+1
  return C

# GetOverlap:  Reads in packed symmetric column matrix, number of basis functions.
# Output:     -Overlap Matrix (NBasis,NBasis)

def GetOverlap(A,NBasis):
   C = column2square(A,NBasis)
   CInv = np.linalg.inv(C)
   S = np.dot(np.transpose(CInv),CInv)
   return S  

# PrintSI: Reads in filename, user-defined switch
# Output: -SCF Energy, Charge, Multiplicity, Geometry
#
# Switch:  1 = print to new file (filename1-SI.txt)
#         -1 = print to screen
#

def PrintSI(filename1,switch):
    NBasis, NElementsGrab, Charge, Multiplicity, NAtoms, SCFEnergy  = NBasGrab(filename1)
    AtomicNum = GetAtoms(filename1,NAtoms)
    RawCart = GeomGet(filename1,NAtoms)
    Cart = np.resize(RawCart,(NAtoms,3))
    filename2 = os.path.splitext(filename1)[0] + "-SI.txt"
    if (switch == 1):
       with open(filename2,'w') as f2:
         f2.write("SI info for ")
         f2.write(filename1)
         f2.write("\n\n")
         f2.write("SCF Energy = ")
         f2.write(str(SCFEnergy))
         f2.write(" Hartree")
         f2.write("\n\n")
         f2.write(str(Charge))
         f2.write(" ")
         f2.write(str(Multiplicity))
         f2.write("\n")
         for i in range(0,NAtoms):
             h = i + 1
             z = AtomicNum[i]
             Atom  = AtomicSymbol(int(z))
             f2.write(Atom)
             f2.write(" ")
             for j in range(0,3):
                 if (Cart[i,j] >= 0):
                    f2.write(" ")
                 f2.write(str(sci_notation(Cart[i,j])))
                 f2.write(" ")
             f2.write("\n")
         f2.write(" ")
    if (switch == -1):
      print "SCF Energy = ", SCFEnergy, " Hartree\n"
      print "Charge = ", Charge, "\n"
      print "Multiplicity = ", Multiplicity, "\n"
      print "Cartesian Geometry:\n"
      for i in range(0,NAtoms):
          h = i + 1
          z = AtomicNum[i]
          Atom = AtomicSymbol(int(z))
          print Atom, sci_notation(Cart[i,0]), sci_notation(Cart[i,1]), sci_notation(Cart[i,2])
      print "\n"    
      
# CalcNO: Reads in filename, NBasis
# Output: Natural Orbitals eigenvalues and eigenvectors (both alpha and beta)
# 

def CalcNO(filename,NBasis):
       Palpha, Pbeta = MatGrab(filename,NBasis,2)   
       C, a = MatGrab(filename,NBasis,1)
       S = GetOverlap(C,NBasis)
       Svals, Svecs = np.linalg.eig(S)
       Sval_minhalf = (np.diag(Svals**(0.5)))
       Shalf = np.dot(Svecs,np.dot(Sval_minhalf,np.transpose(Svecs)))
       NOvalsA, NOvecsA = np.linalg.eig(np.dot(Shalf,np.dot(Shalf,Palpha)))
       NOvalsB, NOvecsB = np.linalg.eig(np.dot(Shalf,np.dot(Shalf,Pbeta))) 
       NOvalsA = NOvalsA.real
       NOvalsB = NOvalsB.real 
       NOvecsA = NOvecsA.real
       NOvecsB = NOvecsB.real
       print "Alpha Natural Orbitals Eigenvectors =\n", NOvecsA
       print "Alpha Natural Orbitals Eigenvalues  =\n", NOvalsA
       print "Beta  Natural Orbitals Eigenvectors  =\n", NOvecsB
       print "Beta  Natural Orbitals Eigenvalues   =\n", NOvalsB
       return NOvecsA, NOvecsB, NOvalsA, NOvalsB

# NElec: Reads in filename
# Output: Total number of electrons, Alpha Electrons, Beta Electrons
#

def NElec(filename):
  NElec = 0
  NAlpha = 0
  NBeta = 0
  with open(filename, 'r') as origin:
     for line in origin:
        if "Number of electrons" in line:
           words = line.split()
           for i in words:
              for letter in i:
                 if(letter.isdigit()):
                    NElec = NElec*10 + int(letter)
        if "Number of alpha electrons" in line:
           words = line.split()
           for i in words:
              for letter in i:
                 if(letter.isdigit()):
                    NAlpha = NAlpha*10 + int(letter)
        if "Number of beta electrons" in line:
           words = line.split()
           for i in words:
              for letter in i:
                 if(letter.isdigit()):
                    NBeta = NBeta*10 + int(letter)
  return NElec, NAlpha, NBeta

# OrbTransform: Reads in Alpha Density Matrix, Beta Density Matrix, Overlap Matrix, n
# Output: New Density Matrices: P' = S**(1-n).P.S**(n)
#

def OrbTransform(Pa,Pb,S,n):
       Svals, Svecs = np.linalg.eig(S)
       Sval1 = np.diag(Svals**(n))
       Sval2 = np.diag(Svals**(1-n)) 
       Sdag1 = np.dot(Svecs,np.dot(Sval1,np.transpose(Svecs)))
       Sdag2 = np.dot(Svecs,np.dot(Sval2,np.transpose(Svecs)))
       PdagAlpha = np.dot(Sdag1,np.dot(Pa,Sdag2))
       PdagBeta = np.dot(Sdag1,np.dot(Pb,Sdag2))
       print "OrbTransform Subroutine test:\n"
       print "PdagAlpha = ", PdagAlpha, "\n"
       print "PdagBeta = ", PdagBeta, "\n"
       OvalsA, OvecsA = np.linalg.eig(PdagAlpha)
       OvalsB, OvecsB = np.linalg.eig(PdagBeta) 
       print "OVals A = ", OvalsA, "\n"
       print "OVecs A = ", OvecsA, "\n"
       print "OVals B = ", OvalsB, "\n"
       print "OVecs B = ", OvecsB, "\n"
       return PdagAlpha, PdagBeta, OvecsA, OvecsB, OvalsA, OvalsB

# CartoZmat: Transforms Cartesian coordinates to z-matrix form
# Input: NAtoms, RawCart, AtomicNum
# Output: z-matrix file
#

def DistAB(e1,e2):
    R = 0.0
    for i in range(0,3):
      R = R + (e1[i]-e2[i])**2
    R = math.sqrt(R) 
#    Rab = math.sqrt((xb-xa)**2+(yb-ya)**2+(zb-za)**2)
#    print "Distance = ", Rab
    print " R = ", R
    return R

def AngleABC(e1,e2,e3):
    eab_x = (e1[0] - e2[0]) / DistAB(e1,e2)
    eab_y = (e1[1] - e2[1]) / DistAB(e1,e2)
    eab_z = (e1[2] - e2[2]) / DistAB(e1,e2)

    ebc_x = (e3[0] - e2[0]) / DistAB(e2,e3)
    ebc_y = (e3[1] - e2[1]) / DistAB(e2,e3)
    ebc_z = (e3[2] - e2[2]) / DistAB(e2,e3)

    eab = [eab_x, eab_y, eab_z]
    ebc = [ebc_x, ebc_y, ebc_z]
   
    cos_angle = np.dot(eab,ebc)
    angle = np.arccos(cos_angle) / 3.1415926535 * 180
    print "Cosine of angle =" , cos_angle
    print "Angle = ", angle
    return eab, ebc, angle    

def TorsionABCD(e1,e2,e3,e4):
    eab, ebc, angle1 = AngleABC(e1,e2,e3)
    ebc, ecd, angle2 = AngleABC(e2,e3,e4)
    cos_angle = np.dot(np.cross(eab,ebc),np.cross(ebc,ecd)) / (np.sin(angle1)*np.sin(angle2))
    angle = np.arccos(cos_angle) / 3.1415926535 * 180
    print "Cos angle = ", cos_angle
    print "The Dihedral angle is =", angle

