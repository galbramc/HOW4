from __future__ import division
import numpy as npy
from numpy import pi, sin, cos, tan, log10
from scipy.interpolate import interp1d
from scipy import integrate
from scipy import optimize
from plot3d import writePlot2D, writePlot3D

import pylab as pyl

def Distance(i, dx_te, ratio):
    return dx_te*(ratio**i - 1) / (ratio-1)

def GradDist(i, dx_te, ratio):
    return dx_te*((i-1)*ratio**i - i*ratio**(i-1) + 1) / (ratio-1)**2

def FindStretching(n, h_min, Hc):
    # Find the ratio of successive cell sizes to get a total length of Hc from
    # n cells, starting with the first cell size = h_min.
    
    # This guess is exact for very large size ratios
    guess = (Hc / h_min)**(1/(n-1))
    finished = False
    while not finished:
        func = Distance(n, h_min, guess) - Hc
        grad = GradDist(n, h_min, guess)
        delta = - func / grad
        guess += delta
        #print guess, ' ', delta
        if (abs(delta) < 1.e-12):
            finished = True

    return guess;
    

def make_airfoil(Dfarfield, ref, Q, TriFlag, farang=0.0, nchordwise=20,
                 nxwake=9, rxwakecenter=3.0, rxwakefary=0.35, nnormal=14,
                 rnormal=2.8, rnormalfar=3.0, TEfac=1.0, Ufac=1.0,
                 wakeangle=0.0, reynolds=1.e6, filename_base="Joukowski"):
    # function make_airfoil(xyfile, UseExact, Dfarfield, ref, Q, TriFlag, farang, varargin)
    #
    # Makes a quad or tri .gri mesh for an airfoil using the points 
    # supplied in the file xyfile.  This file must have two numbers
    # per line, each representing an (x,y) coordinate of a point.
    # The points should start at the trailing edge and loop clockwise.  
    # The trailing-edge is assumed closed (no gap), and the trailing-
    # edge point should not be repeated.  The number of points in
    # the file should be sufficient to represent the geometry well, but
    # it need not be a multiple of Q as the points will be re-splined.
    # An optional hard-coded analytical geometry function can be used to
    # nudge points to the true geometry (if using the spline is not enough).
    # The spacing of points on the geometry is done via a quasi-curvature
    # based method -- the optional mesh size input controls this.
    # The generated mesh is of the "C" type (see graphic make_airfoil.png).
    #
    # INPUTS:
    #   Dfarfield : approximate distance to the farfield, in chords (e.g. 50)
    #   ref       : refinement number (useful for convergence studies)
    #   Q         : geometry order (e.g. 1,2,3,4,...)
    #   TriFlag   : 0 = quad, 1 = tri
    #   farang    : angle from horizontal of farfield inflow at min/max y
    #               (useful to keep inflow a true inflow for nonzero alpha)
    #   varargin  : mesh size/spacing structure (optional, else default one
    #               defined below will be used).
    #
    # OUTPUTS:
    #   A file named: airfoil.gri
    
    #--------------------------------------------#
    # pull off meshsize structure from variable  # 
    # input argument; or define it here          #
    #--------------------------------------------#
    
    # nchordwise = 20    # number of elements along one side of the airfoil geometry
    # nxwake = 9         # x-wake on centerline
    # rxwakecenter = 3.0 # x-wake stretching on centerline
    # rxwakefary = 0.35  # x-wake stretching far from airfoil (+/-y)
    # nnormal = 14       # points normal to airfoil surface, 14=viscous, 12=inviscid
    # rnormal = 2.8      # normal-direction stretching, x close to airfoil
    # rnormalfar = 3.0   # normal-direction stretching, x far from airfoil
    
    # TEfac = 1.0        # factor controlling bunching at TE (high = bunched)
    # Ufac  = 1.0        # factor controlling uniformity of chordwise spacing (high = uniform)
    # wakeangle=0.0      # angle of wake leaving the airfoil
    
    #--------------------#
    # load/spline points #
    #--------------------#
    X = Joukowski(nchordwise*2**ref,Q)
    # print X;
    
    c = max(X[:,0]) - min(X[:,0])          # chord length
    Hc = Dfarfield*c                       # farfield distance
    
    xte = X[0,:];                          # TE point
    XLE = X[npy.append(range(len(X)),0),:] # rest of airfoil
    nLE = len(XLE)
    
    #-------------------------------------#
    # put points down along farfield, FLE #
    #-------------------------------------#
    dx = XLE[1:,:]-XLE[:-1,:]
    ds = (dx[:,0]**2 + dx[:,1]**2)**0.5+0.1
    s  = npy.zeros(nLE)
    for i in xrange(1,nLE):
        s[i] = s[i-1]+ds[i-1]
    x0     = tan(farang)*Hc
    radius = (x0**2 + Hc**2)**0.5
    t0     = s/max(s)*(pi-2*farang) + 3*pi/2 + farang
    FLE    = npy.zeros([nLE,2])
    FLE[:,0] = x0 - radius*cos(t0)
    FLE[:,1] =      radius*sin(t0)
    
    #----------------------#
    # x-wake on centerline #
    #----------------------#
    nr0 = nxwake*2**ref # 9=inviscid, 11=viscous
    a   = 0.1
    b   = rxwakecenter  # 2.9 for NACA, inviscid
    dx_te = X[0,0] - X[Q,0];
    #print "TE locations\n"
    #print dx_te, X[0,0], X[Q,0]
    #print "TE locations done\n"
    re  = (npy.logspace(a,b,nr0+1) - 10**a)/(10**b-10**a)
    #print re;

    ratio = FindStretching(nr0, dx_te, Hc)
    for i in xrange(0, nr0+1):
        re[i] = Distance(i, dx_te, ratio)/Hc
    #print re;

    rw  = spaceq(re, ref, Q)
    
    #----------------------------------#
    # C-grid: put points on wake first #
    #----------------------------------#
    
    XWK = npy.flipud(npy.array([rw*Hc+xte[0], npy.zeros(len(rw))]).transpose())
    XWK[:,1] = (XWK[:,0]-xte[0])*tan(wakeangle)
    XWK2 = npy.flipud(XWK)

    nWK = len(XWK)
    
    #----------------------------------------#
    # x-wake spacing far from airfoil (+/-y) #
    #----------------------------------------#
    a  = 0.1
    b  = rxwakefary
    re = (npy.logspace(a,b,nr0+1) - 10**a)/(10**b-10**a)
    rbot = npy.flipud(spaceq(re, ref, Q)*(Hc+xte[0]))
    
    FWK1 = npy.array([rbot,              XWK[:,1]  - Hc - rbot*x0/Hc]).transpose()
    FWK2 = npy.array([npy.flipud(rbot), XWK2[:,1] + Hc + npy.flipud(rbot)*x0/Hc]).transpose()
    
    #-------------------#
    # Wake and boundary #
    #-------------------#
    XWB = npy.append(XWK,  XLE[1:-1,:], axis = 0)
    XWB = npy.append(XWB,  XWK2,        axis = 0)
    FWB = npy.append(FWK1, FLE[1:-1,:], axis = 0)
    FWB = npy.append(FWB,  FWK2,        axis = 0)
    
    nWB = len(XWB)
    
    #------------------#
    # points on C grid #
    #------------------#
    nr0 = nnormal*2**ref

    # The old spacing; log-linear
    a  = 0.1
    b  = rnormal
    re = (npy.logspace(a,b,nr0+1) - 10**a)/(10**b-10**a)
    # print re

    # The new spacing; exponential
    if (reynolds > 5e5):
        # Turbulent.  y+=5 for the first cell at the TE on the coarse mesh
        coarse_yplus = 5
        dy_te = 5.82 * (coarse_yplus / reynolds**0.9) / 2**ref
        wake_power = 0.8
    else:
        # Laminar.  Put two cells across the BL at the TE on the coarse mesh
        dy_te = 1. / reynolds**0.5 / 2**ref
        wake_power = 0.5

    nr = 1 + (len(re)-1)*Q
    XC = npy.zeros([nWB, nr])
    YC = npy.array(XC)
    a  = 0.1
    b  = rnormalfar
    re = (npy.logspace(a,b,nr0+1) - 10**a)/(10**b-10**a)
    r1 = spaceq(re, ref, Q)
    #print re
    
    for i in xrange(nWB):
        iplus = min(nWB-1, i+1)
        iminus = max(0, i-1)
        ds = ((XWB[iplus,0] - XWB[iminus,0])**2 +
              (XWB[iplus,1] - XWB[iminus,1])**2)**0.5/ (iplus - iminus)
        dy = min(dy_te, ds*5*Q)
        dy = dy_te * max(XWB[i,0],1)**wake_power
        # print XWB[iplus,0], XWB[iminus,0], ds, dy, iplus, iminus
        ratio = FindStretching(nr0, dy, Hc)
        for i2 in xrange(0, nr0+1):
            #print i2, Distance(i2, dy, ratio)/Hc
            re[i2] = Distance(i2, dy, ratio)/Hc
        r0 = spaceq(re, ref, Q)
        #print re
    
        r = r0
        #if i < nWK-1 or i > nWB-nWK-1:
        #    xx = (XWB[i,0]-XWK[-1,0])/max(XWB[:,0])
        #    r = r0 * (1-xx) + r1 * xx
        XC[i,:] = XWB[i,0] + r*(FWB[i,0]-XWB[i,0])
        YC[i,:] = XWB[i,1] + r*(FWB[i,1]-XWB[i,1])
    
    #pyl.plot(XC,YC,'o')
    #pyl.show()
   
    writePlot2D(filename_base + '_ref'+str(ref)+ '_Q'+str(Q)+'.p2d', XC, YC)
    writePlot3D(filename_base + '_ref'+str(ref)+ '_Q'+str(Q)+'.p3d', XC, YC)
 
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

    #writeGRM(filename_base, ref, Q, E, V, nLE, NC, nWK, nWB, nr);
    #writeVTK(filename_base, ref, Q, E, V);
    #writeFEC(filename_base, ref, Q, E, V, nLE, NC, nWK, nWB, nr);

    return
    
    #-----------------------------------
def block_elem(N, Q):
    nx, ny = N.shape;
    #if (Q != 1) and ((mod(nx,Q) != 1) or (mod(ny,Q) != 1)): print('ERROR 2'); return;
    mx = int((nx-1)/Q);
    my = int((ny-1)/Q);
    E = npy.zeros( (mx*my,(Q+1)*(Q+1)),int);
    dy = ny;
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

def writeGRM(filename_base, ref, Q, E, V, nLE, NC, nWK, nWB, nr):
    #=========================#
    # Write out the grid file #
    #=========================#

    f = open(filename_base + '_ref'+str(ref)+ '_Q'+str(Q)+'.grm', 'w')
     
    nelem = E.shape[0];
    #if (TriFlag): nelem = nelem*2;
    print 'Elements : ', nelem
    nnode = V.shape[0];
    f.write('2 ' + str(nnode) + ' 1 3\n') #dim nNodes negrp nbfgrp
     
    #----------#
    # Vertices #
    #----------#
    floatformat = "{:3.16e}"
    
    for i in xrange(nnode):
        f.write(floatformat.format(V[i,0]) + ' ' + floatformat.format(V[i,1]) + '\n')
      
    #----------------#
    # Boundary faces #
    #----------------#
        
    # Airfoil
    nb = int((nLE-1)/Q) #(nLE-1)/2+(nTE-1)/2;
    f.write(str(nb) + '\n');
    f.write('PXE_Shape_Edge\n')
    for i in xrange(int((nLE-1)/Q)):
        f.write(str(NC[nWK-1+Q*i,0]) + ' ' + str(NC[nWK-1+Q*(i+1),0]) + '\n')
    
      
    # Farfield inflow
    f.write(str(int((nWB-1)/Q)) + '\n')
    f.write('PXE_Shape_Edge\n')
    for i in xrange(int((nWB-1)/Q)):
      f.write(str(NC[Q*i,nr-1]) + ' ' + str(NC[Q*(i+1),nr-1]) + '\n')
      
    # Farfield Outflow
    nb = int(2*(nr-1)/Q);
    f.write(str(nb) + '\n')
    f.write('PXE_Shape_Edge\n')
    for i in xrange(int((nr-1)/Q)):
        f.write(str(NC[0,Q*i]) + ' ' + str(NC[0,Q*(i+1)]) + '\n')

    for i in xrange(int((nr-1)/Q)):
        f.write(str(NC[nWB-1,Q*i]) + ' ' + str(NC[nWB-1,Q*(i+1)]) + '\n')
    
      
      
    #----------#
    # Elements #
    #----------#
      
#     if (TriFlag)
#       fprintf(fid, '%d %d TriLagrange\n', nelem, Q);
#       j=1; for ic=(Q+1):-1:1, for ir=1:ic,N1a(j)=(ir-1)*(Q+1)+ic; j=j+1; end; end;
#       j=1; for ic=1:(Q+1), for ir=(Q+1):-1:ic,N2a(j)=(ir-1)*(Q+1)+ic; j=j+1; end; end;
#       j=1; for ir=1:(Q+1), for ic=1:(Q+2-ir),N1b(j)=(ir-1)*(Q+1)+ic; j=j+1; end; end;
#       j=1; for ic=(Q+1):-1:1, for ir=(Q+2-ic):Q+1,N2b(j)=(ir-1)*(Q+1)+ic; j=j+1; end; end;
#       for k=1:nelem/2,
#         xy = V(E(k,6), :); % a node inside the elem
#         if (xy(2) > 0.)
#           N1 = N1a; 
#           N2 = N2a; 
#         else
#           N1 = N1b;
#           N2 = N2b;
#         end
#         fprintf(fid, '%d %d %d %d %d %d %d %d %d %d\n', E(k,N1));
#         fprintf(fid, '%d %d %d %d %d %d %d %d %d %d\n', E(k,N2));
#       end
#     else
    f.write(str(nelem)+'\n')
    f.write(str(Q)+'\n')
    f.write('PXE_Shape_Quad\n')
    f.write('UniformNodeDistribution\n')
    for e in xrange(nelem):
        for k in xrange((Q+1)*(Q+1)):
            f.write(str(E[e,k])+' ')
        f.write('\n')
      
    f.close()
    return

def writeVTK(filename_base, ref, Q, E, V):
    #=========================#
    # Write out the grid file #
    #=========================#

    f = open(filename_base + '_ref'+str(ref)+ '_Q'+str(Q)+'.vtk', 'w')

    f.write('# vtk DataFile Version 2\n');
    f.write(filename_base + ', level ' + str(ref) + ' order ' + str(Q) + '\n');
    f.write('ASCII\n\n');
    f.write('DATASET UNSTRUCTURED_GRID\n');
    
    nelem = E.shape[0];
    nnode = V.shape[0];

    #----------#
    # Vertices #
    #----------#
    f.write('POINTS ' + str(nnode) + ' float\n');
    floatformat = "{:3.16e}"
    for i in xrange(nnode):
        f.write(floatformat.format(V[i,0]) + ' ' +
                floatformat.format(V[i,1]) + ' 0\n')
      
    #----------#
    # Elements #
    #----------#
    # Write as linear quads
    
    f.write('CELLS '+ str(nelem)+ ' ' + str(5*nelem) + '\n');
    for e in xrange(nelem):
        f.write('4   ');
        f.write(str(E[e,1-1]-1)+' ');
        f.write(str(E[e,Q+1-1]-1)+' ');
        f.write(str(E[e,(Q+1)*(Q+1)-1]-1)+' ');
        f.write(str(E[e,Q*(Q+1)+1-1]-1)+' ');
        f.write('\n');

    #------------------------------------------------#
    # Element types: Linear tri = 5; linear quad = 9 #
    #------------------------------------------------#
    f.write('CELL_TYPES ' + str(nelem) + '\n');
    for e in xrange(nelem):
        f.write('9\n');
                
    f.close()
    return

def writeFEC(filename_base, ref, Q, E, V, nLE, NC, nWK, nWB, nr):
    #=========================#
    # Write out the grid file #
    # UBC curved FE format    #
    #                         #
    # Also writes the boundary#
    # shape in GRUMMP .bdry   #
    # format                  #
    #=========================#

    if (Q != 3):
        print("Error: can't write cubic FE data with Q = " + str(Q) + '\n');
        return;
    
    f = open(filename_base + '_ref'+str(ref)+'.fec', 'w')

    nelem = E.shape[0];
    nnode = V.shape[0];

    f.write(str(nnode) + ' ' + str(nelem) + ' 4\n');
    
    #----------#
    # Vertices #
    #----------#
    floatformat = "{:3.16e}"
    for i in xrange(nnode):
        f.write(floatformat.format(V[i,0]) + ' ' +
                floatformat.format(V[i,1]) + '\n')
      
    #----------#
    # Elements #
    #----------#
    # Write as cubic quads
    
    for e in xrange(nelem):
        f.write('12   ');
        # First write the four corners
        f.write(str(E[e,0]-1)+' ');
        f.write(str(E[e,3]-1)+' ');
        f.write(str(E[e,15]-1)+' ');
        f.write(str(E[e,12]-1)+' ');
        # Now write the edge nodes, CCW from vert 0
        f.write(str(E[e,1]-1)+' ');
        f.write(str(E[e,2]-1)+' ');
        f.write(str(E[e,7]-1)+' ');
        f.write(str(E[e,11]-1)+' ');
        f.write(str(E[e,14]-1)+' ');
        f.write(str(E[e,13]-1)+' ');
        f.write(str(E[e,8]-1)+' ');
        f.write(str(E[e,4]-1)+' ');
        f.write('\n');

    f.close()

    f = open(filename_base + '_ref'+str(ref)+'.bdry', 'w')
    #----------------#
    # Boundary faces #
    #----------------#

    # Number of points, number of curves
    f.write(str(nLE + nWB + 2*nr - 4) + ' 3\n');
    # Airfoil coords
    f.write('# airfoil coordinates ' + str(nLE-1) + ' points\n');
    for i in xrange(0,nLE-1):
        index = NC[nWK-2+i,0]
        f.write(floatformat.format(V[index,0]) + ' ' +
                floatformat.format(V[index,1]) + '\n')
#        f.write(str(i) + ' ' + floatformat.format(V[index,0]) + ' ' +
#                floatformat.format(V[index,1]) + '\n')

      
    # Farfield inflow
    f.write('#inflow ' + str(nWB) + ' points\n')
    for i in xrange(nWB):
        index = NC[i, nr-1] - 1
        f.write(floatformat.format(V[index,0]) + ' ' +
                floatformat.format(V[index,1]) + '\n')
#        f.write(str(i  +nLE-1) + ' ' + floatformat.format(V[index,0]) + ' ' +
#                floatformat.format(V[index,1]) + '\n')

    # Farfield Outflow
    nb = int(2*(nr-1));
    f.write('#outflow ' + str(nb) + ' points \n')
    for i in xrange(1,nr-1):
        index = NC[0,nr - i - 1] - 1
        f.write(floatformat.format(V[index,0]) + ' ' +
                floatformat.format(V[index,1]) + '\n')
#        f.write(str(i + nLE-1 + nWB-1) + ' ' + floatformat.format(V[index,0]) + ' ' +
#                floatformat.format(V[index,1]) + '\n')
    f.write('#outflow part 2\n')
    for i in xrange(nr-1):
        index = NC[nWB-1,i] - 1
        f.write(floatformat.format(V[index,0]) + ' ' +
                floatformat.format(V[index,1]) + '\n')
#        f.write(str(i + nLE-1 + nWB-1 + nr) + ' ' + floatformat.format(V[index,0]) + ' ' +
#                floatformat.format(V[index,1]) + '\n')

    # Airfoil surface
    f.write('spline r 1 b 1 ' + str(nLE) + ' ');
    for i in xrange(0, nLE-1):
        f.write(str(i) + ' ')
    f.write(' 0\n')
    
    # Inflow part of domain
    f.write('spline b 5 r 1 ' + str(nWB) + ' ')
    for i in xrange(nWB):
        f.write(str(i+nLE-1) + ' ')
    f.write('\n')
      
    # Outflow end of the domain
    f.write('spline r 1 b 5 ' + str(nb+1) + ' ' + str(nLE-1) + ' ')
    for i in xrange(nb-1):
        f.write(str(i + nLE + nWB - 1) + ' ')
    f.write(str(nLE-1 + nWB-1) + '\n')
    f.close()

    return

#-----------------------------------
def Joukowski_xy(s,a):
    den  = 1 + 2*a*(1 + a)*(1 + cos(pi*s)) ;
    xnum = (1 + a*(1 + 2*a)*(1 + cos(pi*s)))*(sin(0.5*pi*s))**2 ;
    ynum = 0.5*a*(1 + 2*a)*(1 + cos(pi*s))*sin(pi*s) ;
    x = xnum/den ;
    y = ynum/den ;
    
    return x, y

#-----------------------------------
def Joukowski_dxy_ds(s,a):
    den  = 1 + 2*a*(1 + a)*(1 + cos(pi*s)) ;
    xnum = (1 + a*(1 + 2*a)*(1 + cos(pi*s)))*(sin(0.5*pi*s))**2 ;
    ynum = 0.5*a*(1 + 2*a)*(1 + cos(pi*s))*sin(pi*s) ;

    den_ds  = -2*a*(1 + a)*pi*sin(pi*s) ;
    xnum_ds = pi*cos((pi*s)/2.)*(1 + a*(1 + 2*a)*(1 + cos(pi*s)))*sin((pi*s)/2.) - a*(1 + 2*a)*pi*sin((pi*s)/2.)**2*sin(pi*s) ;
    ynum_ds = (a*(1 + 2*a)*pi*cos(pi*s)*(1 + cos(pi*s)))/2. - (a*(1 + 2*a)*pi*sin(pi*s)**2)/2. ;

    dxds = xnum_ds/den - xnum*den_ds/den**2 ;
    dyds = ynum_ds/den - ynum*den_ds/den**2 ;
    
    return dxds, dyds

#-----------------------------------
def Joukowski(nn, Q):
    # hardcoded analytical function
    
    X = npy.zeros([2*nn*Q,2])
    a = 0.1

    # The Joukowski airfoil is already defined in a cosine parametric space,
    # so linspace is correct here, not cos(linspace).
    #s = 1-npy.linspace(0,1,nn+1)
    s = 1-0.5*(1-npy.cos(pi*npy.linspace(0,1,nn+1)))
    #print nn, s
    sL = spaceqarc(s, a, Q)
    #print sL;
    sU = sL[::-1]

    xL, yL = Joukowski_xy(sL,a)
    xU, yU = Joukowski_xy(sU,a)
    yL = -yL
    #print xL;
    
    X[:,0] = npy.append(xL,xU[1:-1])
    X[:,1] = npy.append(yL,yU[1:-1])

    return X

#===============================================================================
def spaceq(re, ref, Q):
    nsub = Q
    nre = len(re) - 1
    nr  = nsub*nre
    r = npy.zeros(nr+1)
    for i in xrange(nre):
        for j in xrange(nsub):
            f = j/nsub
            r[i*nsub+j] = re[i]*(1.0-f) + re[i+1]*f
    r[nr] = re[nre]
    
    return r

#===============================================================================
def spaceqarc(se, a, Q):
    
    def arc(s):
        dxds, dyds = Joukowski_dxy_ds(s,a)
        return npy.sqrt( dxds**2 + dyds**2 )

    nsub = Q
    ns = len(se) - 1
    nr  = nsub*ns
    s = npy.zeros(nr+1)
    for i in xrange(ns):
        
        arclength = integrate.quad( arc, se[i], se[i+1] )[0]
        
        s[i*nsub] = se[i]
        for j in xrange(1,nsub):
            f = j/float(nsub)
            
            s[i*nsub+j] = optimize.bisect(lambda t:integrate.quad( arc, se[i], t )[0]-arclength*f, se[i] + 1e-8*arclength, se[i+1]-1e-8*arclength)

    s[nr] = se[ns]
    
    return s



if __name__ == '__main__':
    Q = 3
    for ref in xrange(0,3):
        make_airfoil(100, ref, Q, False, nchordwise=8, nxwake=8, nnormal=14,
                     rnormal=4, rnormalfar=4, rxwakecenter=3.65, reynolds=1.e6,
                     filename_base="Joukowski")
        print("Done with level " + str(ref));
    #import pylab as pyl
    #X = Joukowski(500, 1)
    #pyl.plot(X[:,0],X[:,1],'-o')
    #pyl.axis( [0,1,-0.5,0.5] )
    #pyl.show()
