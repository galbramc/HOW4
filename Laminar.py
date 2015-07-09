from Joukowski import make_airfoil

# Q is the degree of the polynomial used to represent elements. For Finite Volume/Difference codes, this should be Q=1 for linear elements.
# Finite Element codes are encouraged to use super-parametric elements with Q=4, or the highest available
Q = 4

#The range of refinement levels to generate
refmin = 0
refmax = 5

#Set to True for triangle grids, and False for qauds
TriFlag=False

for ref in xrange(refmin,refmax+1):
    make_airfoil(100, ref, Q, TriFlag, FileFormat="grm", nchordwise=8, nxwake=8, nnormal=16,
                 rnormal=4, rnormalfar=4, rxwakecenter=3.65, reynolds=1000,
                 filename_base="Joukowski_Laminar")
    
