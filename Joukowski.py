from __future__ import division
import numpy as npy
from plot3d import writePlot2D, writeLaballiur, writePlot3D, writePlot3Dxz, writeOVERFLOW
from grm import writeGRM
from vtk import writeVTK
from fec import writeFEC
from gmsh import writeGMSH
from ebg import writeEBG
from geo import writeGEO

from Joukowski_Challenge import make_joukowski_challenge
from Joukowski_Classic import make_joukowski_classic

import matplotlib.pyplot as plt

def make_joukowski(ref, Q, TriFlag, Distribution, FileFormat, reynolds, filename_base):

    if Distribution == "Classic":
        XC, YC = make_joukowski_classic(ref, Q, reynolds)
    elif Distribution == "Challenge":
        XC, YC = make_joukowski_challenge(ref, Q, reynolds)
    else:
        raise ValueError("Distribution should be Classic or Challenge")

    nLE = 16*Q*2**ref+1
    nWK = 8*Q*2**ref+1
    nWB = XC.shape[0]
    nr  = XC.shape[1]

    fac = 2 if TriFlag else 1
    print 'Cell size ' + str( int((nWB-1)/Q) ) + 'x' + str( int((nr-1)/Q) ) + ' with '  + str( fac*int((nWB-1)/Q)*int((nr-1)/Q) ) + ' Elements'
    
    if FileFormat == 'p2d':
        writePlot2D(filename_base + '_ref'+str(ref)+ '_Q'+str(Q)+'.p2d.x', XC, YC)
    if FileFormat == 'labl':
        assert Q == 1
        writeLaballiur(filename_base + '_ref'+str(ref)+ '_Q'+str(Q)+'.labl', XC, YC, nWK)
    if FileFormat == 'p3dxy':
        writeNMF(filename_base + '_ref'+str(ref)+ '_Q'+str(Q)+'.nmf', XC, nLE, nWK, nWB, nr, 'z')
        writePlot3D(filename_base + '_ref'+str(ref)+ '_Q'+str(Q)+'.p3d', XC, YC)
    if FileFormat == 'p3dxz':
        writeNMF(filename_base + '_ref'+str(ref)+ '_Q'+str(Q)+'.nmf', XC, nLE, nWK, nWB, nr, 'y')
        writePlot3Dxz(filename_base + '_ref'+str(ref)+ '_Q'+str(Q)+'.p3d', XC, YC)
    if FileFormat == 'in':
        writeOVERFLOW('grid.in.'+str(ref), XC, YC)
    if FileFormat == 'hypgen':
        writePlot2D('joukowski_c.crv', XC[:,0:1], YC[:,0:1])
    if FileFormat == 'ebg':
        writeEBG('joukowski.ebg', XC, YC, nWK)
    if FileFormat == 'geo':
        writeGEO('joukowski.geo', XC, YC, nWK)
         
    
    #--------------------#
    # Vertices, unrolled #
    #--------------------#
    V = npy.zeros((nWB*nr,2),float)
    V[:,0] = XC.T.reshape(nWB*nr)
    V[:,1] = YC.T.reshape(nWB*nr)
    
    #pyl.plot(XC.reshape(nWB*nr),YC.reshape(nWB*nr),'o')
    #pyl.show()

    #pyl.plot(V[:,0],V[:,1],'o')
    #pyl.show()

    #---------------------------------------------#
    # node number matrices for writing out blocks #
    #---------------------------------------------#

    NC = npy.arange(nWB*nr).reshape( (nr, nWB) ).T+1
    V = npy.delete(V,NC[nWB-nWK:nWB,0]-1,0)
    NC[nWB-nWK:nWB,0] = NC[nWK-1::-1,0]
    NC[:,1:] = NC[:,1:]-nWK
    
    #---------------#
    # form elements #
    #---------------#
    E = block_elem(NC, Q);

    #---------------#
    # write file    #
    #---------------#

    if FileFormat == 'grm':
        writeGRM(filename_base, ref, Q, TriFlag, E, V, nLE, NC, nWK, nWB, nr);
    if FileFormat == 'fec':
        writeVTK(filename_base, ref, Q, E, V);
        writeFEC(filename_base, ref, Q, E, V, nLE, NC, nWK, nWB, nr);
    if FileFormat == 'msh':
        writeGMSH(filename_base, ref, Q, TriFlag, E, V, nLE, NC, nWK, nWB, nr);

    print("Done with refinement " + str(ref))
    
    
#-----------------------------------
def block_elem(N, Q):
    nx, ny = N.shape;
    #if (Q != 1) and ((mod(nx,Q) != 1) or (mod(ny,Q) != 1)): print('ERROR 2'); return;
    mx = int((nx-1)/Q);
    my = int((ny-1)/Q);
    E = npy.zeros( (mx*my,(Q+1)*(Q+1)),int);
    i = 0;
    for imy in xrange(my):
        for imx in xrange(mx):
            ix = Q*(imx+1)-(Q-1)-1;
            iy = Q*(imy+1)-(Q-1)-1;
            k = 0;
            for ky in xrange(Q+1):
                for kx in xrange(Q+1):
                    E[i,k] = N[ix+kx,iy+ky]
                    k = k+1;

            i = i + 1;
      
    return E

#-----------------------------------------------------------
# writes BC information for FUN3D
def writeNMF(fname, X, nLE, nWK, nWB, nr, sym):

    ni, nj = X.shape; nk = 2

    # NMF file expects all comments and extra empty lines

    f = open(fname, 'w')
    f.write('# ===== Neutral Map File generated by Python  =====\n')
    f.write('# =================================================\n')
    f.write('# Block# IDIM JDIM KDIM\n')
    f.write('# -------------------------------------------------\n')
    f.write('1\n\n')
    f.write('1 ' + str(ni) + ' ' + str(nj) + ' ' + str(nk) + '\n\n')
    f.write('# =================================================\n')
    f.write('# Type         B1 F1 S1 E1 S2 E2 B2 F2 S1 E1 S2 E2 Swap\n')
    f.write('# ---------------------------------------------------\n')
    f.write("'symmetry_" + sym + "'    1 1   1 " + str(ni) + " 1 " +str(nj) + "\n")
    f.write("'symmetry_" + sym + "'    1 2   1 " + str(ni) + " 1 " +str(nj) + "\n")
    f.write("'one-to-one'    1 5   1 " + str(nk) + " 1 " +str(nWK) + "  1 5  1 " + str(nk) + " " + str(ni) + " " + str(nWK+nLE-1) +" False\n")
    f.write("'viscous_solid' 1 5   1 " + str(nk) + " " + str(nWK) + " " + str(nWK+nLE-1) + "\n")
    f.write("'farfield_riem' 1 6   1 " + str(nk) + " 1 " +str(ni) + "\n");
    f.write("'farfield_riem' 1 3   1 " + str(nj) + " 1 " +str(nk) + "\n");
    f.write("'farfield_riem' 1 4   1 " + str(nj) + " 1 " +str(nk) + "\n");



if __name__ == "__main__":
    
    ref = 2
    Q = 1
    reynolds = 1000
    plt.clf()
    plt.axis('equal')
    
    XC, YC = make_joukowski_challenge(ref, Q, reynolds)
    
    plt.pcolor(XC, YC, 0*XC, edgecolor='k', cmap='Greens')
    
    XC, YC = make_joukowski_classic(ref, Q, reynolds)
    
    plt.pcolor(XC, YC, 0*XC, edgecolor='r', cmap='Greens')

    plt.show()
    plt.draw()
    