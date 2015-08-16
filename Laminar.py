from Joukowski import make_joukowski

# Q is the degree of the polynomial used to represent elements. For Finite Volume/Difference codes, this should be Q=1 for linear elements.
# Finite Element codes are encouraged to use super-parametric elements with Q=4, or the highest available
Q = 1

#The range of refinement levels to generate
refmin = 3
refmax = 3

#Set to True for triangle grids, and False for qauds
TriFlag=False

#Used to specify the file format to dump
FileFormat="p2d"

for ref in xrange(refmin,refmax+1):
    make_joukowski(ref, Q, TriFlag, FileFormat, reynolds=1000,
                   filename_base="Joukowski_Laminar")
    
