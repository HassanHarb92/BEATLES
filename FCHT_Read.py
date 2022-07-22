from BEATLES import *


filename = sys.argv[1]
print "Filename = ", filename

NBasis, NElem, Charge, Multiplicity, NAtoms, SCFEnergy = NBasGrab(filename)
print "Number of atoms = ", NAtoms
vibModes = NAtoms*3 - 6
print "Vibrational modes (assuming non linear molecule) =", vibModes
DimDuchinsky = vibModes*vibModes 
print "Duchinsky matrix elements = ", DimDuchinsky


def FCHTDmatGrab(filename,NBasis,switch):
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
#              if "Alpha Orbital Energies" in line:
#                    AOE = i
              if  "FCHT DMat" in line:
                    i=i+1
                    AMO=i
#                    print "Alpha MO coefficients starts at line :", i
                    j=i+MOlines-1
#                    print "Alpha MO coefficients ends at line :", j
                    for m in range(0,j-i+1):
                       nextline = origin.next()
                       nextline = nextline.split()
                       for p in range(p,len(nextline)):
                          MOrawa[r] = nextline[p]
                          r = r+1
                       p = 0
#      print "MO Raw = ", MOrawa
      return MOrawa

RawMatrix = FCHTDmatGrab(filename,vibModes,1)
DuchinskyMatrix = column2square(RawMatrix,vibModes)

DuchinskyMatrix = np.square(DuchinskyMatrix)

print np.trace(DuchinskyMatrix)

absTrace = 0.0
for i in range(0,vibModes):
#   print np.absolute(DuchinskyMatrix[i,i])
   absTrace = absTrace + np.absolute(DuchinskyMatrix[i,i])

#for i in range(0,vibModes):
#    for j in range(0,vibModes):
#        if DuchinskyMatrix[i,j] < 0.00000000000001:
#           print "!!!! ", DuchinskyMatrix[i,j]

print "abs trace = ", absTrace
print "abs trace / vibModes = ", absTrace/vibModes

#### Calculate distance wrt to identity

I = np.identity(vibModes)
print I, np.trace(I)

distance = 0.0

for i in range(0,vibModes):
    for j in range(0,vibModes):
        distance = distance + np.absolute((I[i,j] - DuchinskyMatrix[i,j]))

print "Proposed metric: ", np.around(distance/vibModes,3)
print np.linalg.norm(I)
print np.around(np.linalg.norm(DuchinskyMatrix)/vibModes,3)

