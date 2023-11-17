import numpy as np
import meshio

class FemModel:
    def __init__(self):
        self.coords = []
        self.elnods = []
        self.eldofs = []
        self.bcdofs = []

    def vtu_write_connected_mesh(self,fileName):
        # Write geomtry output to Paraview as a .vtu file
        points = []
        for ip in range(len(self.coords)):
            points.append([self.coords[ip][0], self.coords[ip][1], 0])

        cells = []
        for ie in range(len(self.elnods)):
            numElementNodes = len(self.elnods[ie])
            if numElementNodes == 3 or numElementNodes == 6:
                cells.append(("triangle", [self.elnods[ie][:3]]))
            elif numElementNodes == 4 or numElementNodes == 9:
                cells.append(("quad", [self.elnods[ie][:4]]))

        mesh = meshio.Mesh(points, cells)
        mesh.write("geo_connected_nodes.vtu")


    def vtu_write_stl_style_mesh(self, fileName, dispVector=None, elementCornerStresses=None):
        # Write geomtry output to Paraview as a .vtu file
        # Both 3 and 6 node triangles are written as 3 node triangles
        # Both 4 and 9 node quads are written as 4 node quads
        points = []
        cells = []
        dispArr = []
        stressArr = []

        for iel in range(len(self.elnods)):
            numElementNodes = len(self.elnods[iel])
            if numElementNodes == 3 or numElementNodes == 6:
                numCorners = 3
            elif numElementNodes == 4 or numElementNodes == 9:
                numCorners = 4

            # Add the points
            for ip in range(numCorners):
                inode = self.elnods[iel][ip]
                points.append([self.coords[inode][0], self.coords[inode][1], 0])

                if dispVector is not None:
                    dispArr.append([dispVector[inode*2,0],dispVector[inode*2+1,0],0])

                if elementCornerStresses is not None:
                    stressArr.append(elementCornerStresses[iel][ip])

            ip = len(points) - numCorners
            if numCorners == 3:
                cells.append(("triangle", [[ip,ip+1,ip+2]]))
            elif numCorners == 4:
                cells.append(("quad", [[ip,ip+1,ip+2,ip+3]]))

        point_data = {}
        if dispVector is not None:
            point_data["displacement"] = dispArr
        if elementCornerStresses is not None:
            point_data["stress"] = stressArr

        mesh = meshio.Mesh(points, cells, point_data=point_data)
        mesh.write(fileName)


class CantileverModel(FemModel):
    def __init__(self,length, heigth, numElementNodes, numNodesX, numNodesY):
        FemModel.__init__(self)
        # number of patches that will fit a 9 node element
        numPatchX = (numNodesX-1) // 2
        numPatchX = 1 if numPatchX < 1 else numPatchX
        numPatchY = (numNodesY-1) // 2
        numPatchY = 1 if numPatchY < 1 else numPatchY

        numNodesX = numPatchX*2 + 1
        numNodesY = numPatchY*2 + 1

        if numElementNodes == 6 or numElementNodes == 9:
            numElementsX = (numNodesX-1) // 2
            numElementsY = (numNodesY-1) // 2
        else:
            numElementsX = numNodesX -1
            numElementsY = numNodesY -1

        # Cantilever with dimensions H x L x thickness
        H         = heigth
        L         = length
        thickness =  0.1

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

        numNodes    = numNodesX * numNodesY
        numElements = numElementsX * numElementsY
        if numElementNodes in [3,6]:
            numElements *= 2

        L_elx = L / (numNodesX-1)
        L_ely = H / (numNodesY-1)

        # Set the node coordinates and node dofs

        for i in range(numNodesX):
            for j in range(numNodesY):
                self.coords.append([L_elx * i, L_ely * j])

        # Set the element connectivites and element dofs
        for ip in range(numPatchX):
            ii = ip*2
            for jp in range(numPatchY):
                jj = jp*2
                # 0 based node numbers, 9 nodes of a 3x3 patch
                nod9 = np.array([
                    (ii  )*numNodesY + (jj  ),
                    (ii+1)*numNodesY + (jj  ),
                    (ii+2)*numNodesY + (jj  ),
                    (ii  )*numNodesY + (jj+1),
                    (ii+1)*numNodesY + (jj+1),
                    (ii+2)*numNodesY + (jj+1),
                    (ii  )*numNodesY + (jj+2),
                    (ii+1)*numNodesY + (jj+2),
                    (ii+2)*numNodesY + (jj+2)],'i')

                if numElementNodes == 3:
                    for i in range(2):
                        for j in range(2):
                            self.elnods.append([nod9[3*i+j],nod9[3*i+j+1],nod9[3*(i+1)+j+1]])
                            self.elnods.append([nod9[3*(i+1)+j+1],nod9[3*(i+1)+j],nod9[3*i+j]])
                elif numElementNodes == 6:
                    self.elnods.append([nod9[0],nod9[2],nod9[8],nod9[1],nod9[5],nod9[4]])
                    self.elnods.append([nod9[8],nod9[6],nod9[0],nod9[7],nod9[3],nod9[4]])
                elif numElementNodes == 4:
                    for i in range(2):
                        for j in range(2):
                            self.elnods.append([nod9[3*i+j],nod9[3*i+j+1],nod9[3*(i+1)+j+1],nod9[3*(i+1)+j]])
                elif numElementNodes == 9:
                    self.elnods.append([nod9[0],nod9[2],nod9[8],nod9[6],
                                     nod9[1],nod9[5],nod9[7],nod9[3],
                                     nod9[4]])


        for iel in range(len(self.elnods)):
            dofs = []
            for inod in range(len(self.elnods[iel])):
                dofs.append(self.elnods[iel][inod] * 2) # The x dofs
                dofs.append(self.elnods[iel][inod] * 2 + 1)  # The y dofs

            self.eldofs.append(dofs)


        # Set fixed boundary condition on left side, i.e. nodes 0-nNody
        bc = np.array(np.zeros(numNodesY*2),'i')
        idof = 1
        for i in range(numNodesY):
            self.bcdofs.append(i*2  )  # 0-based dof i.e. fixing x
            self.bcdofs.append(i*2+1)  # 0-based dof i.e. fixing y
