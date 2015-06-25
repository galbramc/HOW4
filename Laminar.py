from Joukowski import make_airfoil

Q = 4

dy_te = 1. / 1000**0.5 / 2**4
print dy_te

for ref in xrange(0,3):
    make_airfoil(100, ref, Q, False, "msh", nchordwise=8, nxwake=8, nnormal=14,
                 rnormal=4, rnormalfar=4, rxwakecenter=3.65, reynolds=1000,
                 filename_base="Joukowski_Laminar")
    print("Done with level " + str(ref));
