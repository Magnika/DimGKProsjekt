# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 16:38:14 2018

@author: bjohau
"""
import numpy as np
#import triangles_LF as tri
import triangles_with_TODO as tri
#import quads_LF_bh as quad
import quads_with_TODO as quad
import fem_utilities as fem_util
import fem_models

# Element Type
numElementNodes = 9  # Valid numbers 3, 4, 6, 9

# Number of nodes: Should be odd numbers in order to handle 9 node quad and 6 node triangle
numNodesX = 21
numNodesY = 13

# Cantilever with dimensions H x L x thickness
H         =  2.0
L         = 10.0
thickness =  0.1

model = fem_models.CantileverModel(L, H, numElementNodes, numNodesX, numNodesY)

# Distributed load in x and y, load pr unit area
eq = np.array([0.,1.0e3])
eq = np.array([0.,0.])
#End load, Given as resultant

endLoadXY = np.array([0.0,3.0e6])
#endLoadXY = np.array([3.0e6,0])
#endLoadXY = np.array([4.2e9,0.0]) # Should give unit disp at Poisson = 0

# Material properties and thickness

ep = [1,thickness]
E  = 2.1e11
nu = 0.3
Dmat = np.array([
        [ 1.0,  nu,  0.],
        [  nu, 1.0,  0.],
        [  0.,  0., (1.0-nu)/2.0]]) * E/(1.0-nu**2)

#model.vtu_write_connected_mesh("GeometryOnly_ConnectedMesh.vtu")
#model.vtu_write_stl_style_mesh("GeometryOnly_STL_style.vtu")

# Assemble stiffness matrix

ndofs = len(model.coords) * 2
K = np.zeros((ndofs,ndofs))
R = np.zeros((ndofs,1))

#Set the load at the right hand edge
for i in range(numNodesY):
    R[-(i*2+2),0] = endLoadXY[0] / numNodesY
    R[-(i*2+1),0] = endLoadXY[1] / numNodesY

numElements = len(model.elnods)
for iel in range(numElements):
    ex_el = []
    ey_el = []

    nElNodes = len(model.elnods[iel])
    for inod in range(nElNodes):
        ex_el.append(model.coords[model.elnods[iel][inod]][0])
        ey_el.append(model.coords[model.elnods[iel][inod]][1])

    if nElNodes == 3:
        K_el, f_el = tri.tri3_Kmatrix(ex_el,ey_el,Dmat,thickness,eq)
    elif nElNodes == 6:
        K_el, f_el = tri.tri6_Kmatrix(ex_el,ey_el,Dmat,thickness,eq)
    elif nElNodes == 4:
        K_el, f_el = quad.quad4_Kmatrix(ex_el,ey_el,Dmat,thickness,eq)
    elif nElNodes == 9:
        K_el, f_el = quad.quad9_Kmatrix(ex_el,ey_el,Dmat,thickness,eq)

    fem_util.assem(model.eldofs[iel], K, K_el, R, f_el)

r, R0 = fem_util.solveq(K, R, model.bcdofs)

# Compute element corner stresses
elementCornerStresses = []
for iel in range(numElements):
    ex_el = []
    ey_el = []

    nElNodes = len(model.elnods[iel])
    for inod in range(nElNodes):
        ex_el.append(model.coords[model.elnods[iel][inod]][0])
        ey_el.append(model.coords[model.elnods[iel][inod]][1])

    nElDofs = len(model.eldofs[iel])
    elDisp = np.zeros(nElDofs)
    for idof in range(nElDofs):
        elDisp[idof] = r[model.eldofs[iel][idof],0]

    if nElNodes == 3:
        cornerStresses = tri.tri3_cornerstresses(ex_el,ey_el,Dmat,thickness,elDisp)
    elif nElNodes == 6:
        #cornerStresses= tri.tri6_cornerstresses(ex_el,ey_el,Dmat,thickness,elDisp)
        cornerStresses = [[1,4,0],[2,5,0],[3,6,0]]
    elif nElNodes == 4:
        #cornerStresses = quad.quad4_cornerstresses(ex_el,ey_el,Dmat,thickness,elDisp)
        cornerStresses = [[1,4,0],[2,5,0],[3,6,0],[0,0,0]]
    elif nElNodes == 9:
        cornerStresses = [[1,4,0],[2,5,0],[3,6,0],[0,0,0]]
        #cornerStresses = quad.quad9_cornerstresses(ex_el,ey_el,Dmat,thickness,elDisp)

    elementCornerStresses.append(cornerStresses)

nodMiddle = numNodesY//2 +1  # Mid nod on right edge
xC = r[-(nodMiddle*2)  ,0] # 2 dofs per node, so this is the middle dof on end
yC = r[-(nodMiddle*2)+1,0] # 2 dofs per node, so this is the middle dof on end
print("Displacement center node right end,  x:{:12.3e}   y:{:12.3e}".format(xC, yC))

# Sum uf reaction forces
R0Sum = np.zeros(2,'f')
for i in range(0,(numNodesY*2),2):
    R0Sum[0] += R0[i  ,0]
    R0Sum[1] += R0[i+1,0]
print("Total reaction force in x:{:12.3e} y:{:12.3e})".format(R0Sum[0],R0Sum[1]))

# Draw the displacements and stresses
model.vtu_write_stl_style_mesh("Results.vtu",dispVector=r,elementCornerStresses=elementCornerStresses)

