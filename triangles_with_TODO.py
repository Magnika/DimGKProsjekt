# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 08:15:51 2018

@author: bjohau
"""
import numpy as np

def tri3_area(ex, ey):
    """
    Compute the area of a triangle element

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :return area of the triangle
    """

    tmp = np.matrix([[1, ex[0], ey[0]],
                     [1, ex[1], ey[1]],
                     [1, ex[2], ey[2]]])

    A2 = np.linalg.det(tmp)  # Double of triangle area
    A = A2 / 2.0
    return A


def tri3_Bmatrix(ex, ey):
    """
    Compute the strain displacement matrix for a 3 node triangular membrane element.

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :return  [3 x 6] strain displacement matrix
    """

    A = tri6_area(ex, ey)
    A2 = 2.0 * A

    cyclic_ijk = [0, 1, 2, 0, 1]  # Cyclic permutation of the nodes i,j,k

    zi_px = np.zeros(3)  # Partial derivative with respect to x
    zi_py = np.zeros(3)  # Partial derivative with respect to y

    for i in range(3):
        j = cyclic_ijk[i + 1]
        k = cyclic_ijk[i + 2]
        zi_px[i] = (ey[j] - ey[k]) / A2
        zi_py[i] = (ex[k] - ex[j]) / A2

    B = np.array([
        [zi_px[0], 0, zi_px[1], 0, zi_px[2], 0],
        [0, zi_py[0], 0, zi_py[1], 0, zi_py[2]],
        [zi_py[0], zi_px[0], zi_py[1], zi_px[1], zi_py[2], zi_px[2]]])

    return B


def tri3_Kmatrix(ex, ey, D, th, eq=None):
    """
    Compute the stiffness matrix for a two dimensional beam element.

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :param list eq: distributed loads, local directions [bx, by]
    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: consistent load vector [6 x 1] (if eq!=None)
    """

    A = tri6_area(ex, ey)
    A2 = 2.0 * A

    cyclic_ijk = [0, 1, 2, 0, 1]  # Cyclic permutation of the nodes i,j,k

    zi_px = np.zeros(3)  # Partial derivative with respect to x
    zi_py = np.zeros(3)  # Partial derivative with respect to y

    for i in range(3):
        j = cyclic_ijk[i + 1]
        k = cyclic_ijk[i + 2]
        zi_px[i] = (ey[j] - ey[k]) / A2
        zi_py[i] = (ex[k] - ex[j]) / A2

    B = tri3_Bmatrix(ex, ey)

    Ke = (B.T @ D @ B) * A * th

    if eq is None:
        return Ke
    else:
        fx = A * th * eq[0] / 3.0
        fy = A * th * eq[1] / 3.0
        fe = np.array([[fx], [fy], [fx], [fy], [fx], [fy]])
        return Ke, fe


def tri3_cornerstresses(ex, ey, D, th, elDispVec):
    """
    Compute the corner stresses for all 3 corner nodes

    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    :param list D : 2D constitutive matrix
    :param list th: element thickness
    :param list eq: distributed loads, local directions [bx, by]
    :return mat Ke: element stiffness matrix [6 x 6]
    :return mat fe: consistent load vector [6 x 1] (if eq!=None)
    """

    B = tri3_Bmatrix(ex, ey)

    strain = B @ elDispVec
    stress = D @ strain

    cornerStresses = []
    for inod in range(3):
        cornerStresses.append([stress[0], stress[1], stress[2]])

    return cornerStresses
    
def zeta_partials_x_and_y(ex,ey):
    """
    Compute partials of area coordinates with respect to x and y.
    
    :param list ex: element x coordinates [x1, x2, x3]
    :param list ey: element y coordinates [y1, y2, y3]
    """
    
    tmp = np.array([[1,ex[0],ey[0]],
                    [1,ex[1],ey[1]],
                    [1,ex[2],ey[2]]])
    
    A2 = np.linalg.det(tmp)  # Double of triangle area
       
    cyclic_ijk = [0,1,2,0,1]      # Cyclic permutation of the nodes i,j,k
    
    zeta_px = np.zeros(3)           # Partial derivative with respect to x
    zeta_py = np.zeros(3)           # Partial derivative with respect to y

    
    #Dette blir bare samme som for 3 noder
    for i in range(3):
        j = cyclic_ijk[i + 1]
        k = cyclic_ijk[i + 2]
        zeta_px[i] = (ey[j] - ey[k]) / A2
        zeta_py[i] = (ex[k] - ex[j]) / A2

    return zeta_px, zeta_py

# Functions for 6 node triangle
    
def tri6_area(ex,ey):
        
    tmp = np.array([[1,ex[0],ey[0]],
                    [1,ex[1],ey[1]],
                    [1,ex[2],ey[2]]])
    
    A = np.linalg.det(tmp) / 2
    
    return A


def tri6_shape_functions(zeta):
    
    cyclic_ijk = [0,1,2,0,1]      # Cyclic permutation of the nodes i,j,k

    N6 = np.zeros(6)

    #Zeta er triangular cordinates fra node 1, 2 og 3

    for i in range(3):
        j = cyclic_ijk[i + 1]
        N6[i] = zeta[i]*(2*zeta[i] - 1)
        N6[i + 3] = 4*zeta[i]*zeta[j]


    return N6


def tri6_shape_function_partials_x_and_y(zeta,ex,ey):
    
    zeta_px, zeta_py = zeta_partials_x_and_y(ex,ey)
    
    N6_px = np.zeros(6)
    N6_py = np.zeros(6)
    
    cyclic_ijk = [0,1,2,0,1]      # Cyclic permutation of the nodes i,j,k

    for i in range(3):          #Bruker dN/dx = dN/dz * dz/dx
        j = cyclic_ijk[i +1]

        N6_px[i] = (4*zeta[i] - 1)*zeta_px[i] 
        N6_py[i] = (4*zeta[i] - 1)*zeta_py[i]

        N6_px[i+3] = 4*(zeta[j]*zeta_px[i] + zeta[i]*zeta_px[j])
        N6_py[i+3] = 4*(zeta[j]*zeta_py[i] + zeta[i]*zeta_py[j])

    return N6_px, N6_py


def tri6_Bmatrix(zeta,ex,ey):
    
    nx,ny = tri6_shape_function_partials_x_and_y(zeta, ex, ey)

    Bmatrix = np.zeros((3,12))
    
    for j in range(6): #Fyller inn B matrisen med dN/dx og dN/dy
        Bmatrix[0][j*2] = nx[j]
        Bmatrix[1][j*2+1] = ny[j]
        Bmatrix[2][j*2] = ny[j]
        Bmatrix[2][j*2+1] = nx[j]


    return Bmatrix


def tri6_Kmatrix(ex,ey,D,th,eq=None):
    
    zetaInt = np.array([[0.5,0.5,0.0],
                        [0.0,0.5,0.5],
                        [0.5,0.0,0.5]])
    
    wInt = np.array([1.0/3.0,1.0/3.0,1.0/3.0]) # Vektene

    A    = tri6_area(ex,ey)

    Ke = np.zeros((12,12))

    Ke = np.zeros((12, 12))

    for i in range(3): #Bruker zetaInt til å få midtpunktene på kantene til trekanten og regner ut B med disse 3.
        zeta = zetaInt[i]
        w = wInt[i]
        B = tri6_Bmatrix(zeta, ex, ey)
        Ke += (B.T @ D @ B) * w * A * th


    if eq is None:
        return Ke
    else:
        fe = np.zeros((12,1))  

        for i in range(3):      #fe = integrale over N*q dV
            N = tri6_shape_functions(zetaInt[i])
            w = wInt[i]
            for j in range(6):
                fe[2*j] = A*th*w*N[j]*eq[0]
                fe[2*j + 1] = A*th*w*N[j]*eq[1]

        return Ke, fe
