from BEATLES import *


filename = sys.argv[1]
NElectrons, NAlpha, NBeta = NElec(filename)
NBasis, NElem, Charge, Multiplicity, NAtoms, SCFEnergy = NBasGrab(filename)

NBas, NBasUse = ifNBasUse(filename)
print "NBasis = ", NBasis, NBas
print "NBasis use = ", NBasUse

NBasis = NBasUse
AlphaMO = MatGrab(filename,NBasis,3)
print "Alpha LUMO energy = ", np.around(AlphaMO[NAlpha]*27.21,2), " eV"
print "Alpha HOMO energy = ", np.around(AlphaMO[NAlpha-1]*27.21,2), " eV"
print "Alpha HOMO-1 energy = ", np.around(AlphaMO[NAlpha-2]*27.21,2), " eV"
