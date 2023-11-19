# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 16:38:14 2018

@author: bjohau
"""

import numpy as np

import triangles_with_TODO as tri

# ----- Topology -------------------------------------------------
ex = np.array([0.,1.,0.]) #element x coordinates [x1, x2, x3]
ey = np.array([0.,0.,1.]) #element y coordinates [y1, y2, y3]

th = 0.1        #Element tykkhet
ep = [1,th]

E  = 2.1e11
nu = 0.3 # Poisson's ratio

D = np.array([          #E matrise i kompendie. s=E*s Hvor s er spenningsvektor og e er tøyningsvektor
        [ 1.0,  nu,  0.],
        [  nu, 1.0,  0.],
        [  0.,  0., (1.0-nu)/2.0]]) * E/(1.0-nu**2)

eq = [1.0, 3.0] #distributed loads, local directions [bx, by]


Ke, fe = tri.tri3_Kmatrix(ex,ey,D,th,eq)

print('Stiffness matrix:\n', Ke)
print('Consistent forces:\n', fe)

rigX = np.array([[1.0],[0.0],[1.0],[0.0],[1.0],[0.0]])
rigY = np.array([[0.0],[1.0],[0.0],[1.0],[0.0],[1.0]])
rigR = np.array([[-ey[0]],[ex[0]],[-ey[1]],[ex[1]],[-ey[2]],[ex[2]]])

fx = Ke @ rigX
fy = Ke @ rigY
fr = Ke @ rigR

print('Force from rigX translation:\n',fx)
print('Force from rigY translation:\n',fy)
print('Force from rigR rotation:\n',fr)
