from __future__ import division

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
