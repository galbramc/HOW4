import numpy as npy


#-----------------------------------------------------------
# writes a curve geometry files
def writeCurve(X, Y, nWK):
    
    # Airfoil
    Afx = X[nWK+1:-nWK,0]
    Afy = Y[nWK+1:-nWK,0]
    nAf = len(Afx)

    # Inflow
    Inx = X[:,-1]
    Iny = Y[:,-1]
    nIn = len(Inx)

    # Outflow Upper
    Outx = X[-1,-2:0:-1]
    Outy = Y[-1,-2:0:-1]

    # Outflow Lower
    Outx = npy.append(Outx, X[0,:-1] )
    Outy = npy.append(Outy, Y[0,:-1])
    nOut = len(Outx)

    filename = 'joukwoski_af.csv'
    Af = open(filename, 'w')
    print 'Writing ', filename

    filename = 'joukwoski_in.csv'
    In = open(filename, 'w')
    print 'Writing ', filename

    filename = 'joukwoski_out.csv'
    Out = open(filename, 'w')
    print 'Writing ', filename

    #----------#
    # Vertices #
    #----------#
    floatformat = "{:3.16e}"
    
    for i in xrange(nAf):
        Af.write(floatformat.format(Afx[i]) + ' ' + floatformat.format(Afy[i]) + '\n')
    Af.close()

    for i in xrange(nIn):
        In.write(floatformat.format(Inx[i]) + ' ' + floatformat.format(Iny[i]) + '\n')
    In.close()

    for i in xrange(nOut):
        Out.write(floatformat.format(Outx[i]) + ' ' + floatformat.format(Outy[i]) + '\n')
    Out.close()
