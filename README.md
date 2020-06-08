# BEATLES
**BEATLES**: **B**undle of **E**ssential and **A**ssistive **T**ools **L**ibrary for **E**lectronic **S**tructure

Python library that interfaces with Gaussian formatted checkpoint files and matrix files.

For details about the available functions, please check the descriptions in the BEATLES.py file.

Several applications to the BEATLES library are also available here:

1. Natural Ionization Orbitals analysis **nio_beatles.py** 
   usage: python nio_beatles.py ground_state.fchk detached_state.fchk
   
2. Corresponding Orbitals analysis **corresponding_orbtials_v1.py**
   usage: python corresponding_orbtials_v1.py molecule.fchk
   
3. Supporting Information Generator: **SI_Generator.py**
   usage: python SI_Generator.py outputfile.txt file1.fchk file2.fchk ....
   
4. Get Basis sets from www.basisssetexchange.org: **get_basis**
   usage: python get_basis.py basisset format element1 element2 ...
   
5. Natural Orbitals Transformation: **NaturalOrbitals.py**
   usage: python NaturalOrbitals.py molecule.fchk
   
6. Mix and Match MOs: **MixandMatch.py**
   usage: python MixandMatch.py molecule1.fchk molecule2.fchk
   
7. Distances Calculator: **Distances.py**
   usage: python Distances.py molecule.fchk
   
**Fun Part**: Everytime you write a code that uses BEATLES.py, make sure to end it with PrintLyrics() to print a random quote from the Beatles when the code runs.
