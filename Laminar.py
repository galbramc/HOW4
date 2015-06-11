from Joukowski import make_airfoil

Q = 1
for ref in xrange(0,4):
    make_airfoil(100, ref, Q, False, "grm", nchordwise=8, nxwake=8, nnormal=14,
                 rnormal=4, rnormalfar=4, rxwakecenter=3.65, reynolds=1000,
                 filename_base="Joukowski")
    print("Done with level " + str(ref));
