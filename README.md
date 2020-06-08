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
   
4. Get Basis sets from www.basissetexchange.org: **get_basis**
   usage: python get_basis.py basisset format element1 element2 ...
   
5. Natural Orbitals Transformation: **NaturalOrbitals.py**
   usage: python NaturalOrbitals.py molecule.fchk
   
6. Mix and Match MOs: **MixandMatch.py**
   usage: python MixandMatch.py molecule1.fchk molecule2.fchk
   
7. Distances Calculator: **Distances.py**
   usage: python Distances.py molecule.fchk
   
**Fun Part**: Everytime you write a code that uses BEATLES.py, make sure to end it with PrintLyrics() to print a random quote from the Beatles when the code runs.


**References**
1. J. Chem. Phys. 144, 204117 (2016); DOI: 10.1063/1.4951738 
2. J. Phys. Chem. Solids 65, 781 (2004). DOI: 10.1016/j.jpcs.2003.11.015
3. Proc. R. Soc. London, Ser. A 263, 483-493 (1961). Online at https://www.jstor.org/stable/2414327
4. J. Chem. Phys. 47, 1936 (1967); DOI: 10.1063/1.1712221
5. Phys. Rev. 101, 1730; DOI: 10.1103/PhysRev.101.1730
6. J. Chem. Inf. Model. 2019, 59(11), 4814-4820, doi:10.1021/acs.jcim.9b00725.

