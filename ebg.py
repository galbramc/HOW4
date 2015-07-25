import numpy as npy


#-----------------------------------------------------------
# writes an ebg geometry file
def writeEBG(filename, X, Y, nWK):
    f = open(filename, 'w')
    print 'Writing ', filename
    
    f.write('nEmbeddedBoundaryFaceGroups\n')
    f.write('4\n')

    # Airfoil
    Vx = X[nWK+1:-nWK,0]
    Vy = Y[nWK+1:-nWK,0]
    nAf = len(Vx)

    # Inflow
    Vx = npy.append(Vx, X[:,-1])
    Vy = npy.append(Vy, Y[:,-1])
    nIn = len(Vx)

    # Outflow Lower
    Vx = npy.append(Vx, X[0,:-1])
    Vy = npy.append(Vy, Y[0,:-1])
    nOutL = len(Vx)

    # Outflow Upper
    Vx = npy.append(Vx, X[-1,1:-1])
    Vy = npy.append(Vy, Y[-1,1:-1])
    nOutU = len(Vx)

    nnode = len(Vx)
    f.write('Vertices\n')
    f.write(str(nnode) + '\n')

    #----------#
    # Vertices #
    #----------#
    floatformat = "{:3.16e}"
    
    for i in xrange(nnode):
        f.write(floatformat.format(Vx[i]) + ' ' + floatformat.format(Vy[i]) + '\n')
      
    #-------#
    # Faces #
    #-------#
    f.write('Faces\n')
    f.write(str(nAf) + '\n')
    for i in xrange(nAf-1):
        f.write(str(i+1) + ' ' + str(i+2) + '\n')
    f.write(str(nAf) + ' 1\n') #Close the loop

    f.write(str(nIn-nAf-1) + '\n')
    for i in xrange(nAf, nIn-1):
        f.write(str(i+1) + ' ' + str(i+2) + '\n')

    f.write(str(nOutL-nIn) + '\n')
    for i in xrange(nIn, nOutL-1):
        f.write(str(i+1) + ' ' + str(i+2) + '\n')
    f.write(str(nOutL) + ' ' + str(nAf+1) + '\n') 

    f.write(str(nOutU-nOutL+1) + '\n')
    f.write(str(nIn+1) + ' ' + str(nOutL+1) + '\n') 
    for i in xrange(nOutL, nOutU-1):
        f.write(str(i+1) + ' ' + str(i+2) + '\n')
    f.write(str(nOutU) + ' ' + str(nIn) + '\n') 

    f.close()