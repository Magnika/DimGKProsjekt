import numpy as np
import sys

def gauss_points(iRule):
    """
    Returns gauss coordinates and weight given integration number

    Parameters:

        iRule = number of integration points

    Returns:

        gp : row-vector containing gauss coordinates
        gw : row-vector containing gauss weight for integration point

    """
    gauss_position = [[ 0.000000000],
                      [-0.577350269,  0.577350269],
                      [-0.774596669,  0.000000000,  0.774596669],
                      [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116],
                      [-0.9061798459, -0.5384693101, 0.0000000000, 0.5384693101, 0.9061798459]]
    gauss_weight   = [[2.000000000],
                      [1.000000000,   1.000000000],
                      [0.555555556,   0.888888889,  0.555555556],
                      [0.3478548451,  0.6521451549, 0.6521451549, 0.3478548451],
                      [0.2369268850,  0.4786286705, 0.5688888889, 0.4786286705, 0.2369268850]]


    if iRule < 1 and iRule > 5:
        sys.exit("Invalid number of integration points.")

    idx = iRule - 1
    return gauss_position[idx], gauss_weight[idx]


def quad4_shapefuncs(xsi, eta):
    """
    Calculates shape functions evaluated at xi, eta
    """
    # ----- Shape functions -----
    # TODO: fill inn values of the  shape functions
    N = np.array([(1+xsi)*(1+eta)/4,(1-xsi)*(1+eta)/4,(1-xsi)*(1-eta)/4,(1+xsi)*(1-eta)/4])  #Defined by zero-lines
    
    return N

def quad4_shapefuncs_grad_xsi(xsi, eta):
    """
    Calculates derivatives of shape functions wrt. xsi
    """
    # ----- Derivatives of shape functions with respect to xsi -----
    # TODO: fill inn values of the  shape functions gradients with respect to xsi

    Ndxi = np.array([(1+eta)/4,-1*(1+eta)/4,-1*(1-eta)/4,(1-eta)/4])   #Calculated by hand
    
    return Ndxi

def quad4_shapefuncs_grad_eta(xsi, eta):
    """
    Calculates derivatives of shape functions wrt. eta
    """
    # ----- Derivatives of shape functions with respect to eta -----
    # TODO: fill inn values of the  shape functions gradients with respect to xsi
    Ndeta = np.array([(1+xsi)/4,(1-xsi)/4,-1*(1-xsi)/4,-1*(1+xsi)/4])   #Calculated by hand
    
    return Ndeta

def quad4_Kmatrix(ex, ey, D, thickness, eq=None):
    """
    Calculates the stiffness matrix for a 8 node isoparametric element in plane stress

    Parameters:

        ex  = [x1 ... x4]           Element coordinates. Row matrix
        ey  = [y1 ... y4]
        D   =           Constitutive matrix
        thickness:      Element thickness
        eq = [bx; by]       bx:     body force in x direction
                            by:     body force in y direction

    Returns:

        Ke : element stiffness matrix (8 x 8)
        fe : equivalent nodal forces (4 x 1)

    """
    t = thickness

    if eq is None:
        f = np.zeros((2,1))  # Create zero matrix for load if load is zero
    else:
        f = np.array([eq]).T  # Convert load to 2x1 matrix

    Ke = np.zeros((8,8))        # Create zero matrix for stiffness matrix
    fe = np.zeros((8,1))        # Create zero matrix for distributed load

    numGaussPoints = 2  # Number of integration points
    gp, gw = gauss_points(numGaussPoints)  # Get integration points and -weight

    for iGauss in range(numGaussPoints):  # Solves for K and fe at all integration points
        for jGauss in range(numGaussPoints):

            xsi = gp[iGauss]
            eta = gp[jGauss]

            Ndxsi = quad4_shapefuncs_grad_xsi(xsi, eta)
            Ndeta = quad4_shapefuncs_grad_eta(xsi, eta)
            N1    = quad4_shapefuncs(xsi, eta)  # Collect shape functions evaluated at xi and eta

            # Matrix H and G defined according to page 52 of Waløens notes
            H = np.transpose([ex, ey])    # Collect global x- and y coordinates in one matrix
            G = np.array([Ndxsi, Ndeta])  # Collect gradients of shape function evaluated at xi and eta

            #TODO: Calculate Jacobian, inverse Jacobian and determinant of the Jacobian
            #J = np.eye(2) #TODO: Correct this
            J = G @ H
            
            invJ = np.linalg.inv(J)  # Inverse of Jacobian
            detJ = np.linalg.det(J)  # Determinant of Jacobian

            dN = invJ @ G  # Derivatives of shape functions with respect to x and y
            dNdx = dN[0]
            dNdy = dN[1]

            # Strain displacement matrix calculated at position xsi, eta

            #TODO: Fill out correct values for strain displacement matrix at current xsi and eta
            B = np.array([[dNdx[0],0,dNdx[1],0,dNdx[2],0,dNdx[3],0],
                        [0,dNdy[0],0,dNdy[1],0,dNdy[2],0,dNdy[3]],
                        [dNdy[0],dNdx[0],dNdy[1],dNdx[1],dNdy[2],dNdx[2],dNdy[3],dNdx[3]]])   #From Waløen page 51


            #TODO: Fill out correct values for displacement interpolation xsi and eta
            N2 = np.array([[N1[0],0,N1[1],0,N1[2],0,N1[3],0],
                        [0,N1[0],0,N1[1],0,N1[2],0,N1[3]]])        #From Waløen page 58

            # Evaluates integrand at current integration points and adds to final solution
            Ke += (B.T) @ D @ B * detJ * t * gw[iGauss] * gw[jGauss] 
            fe += (N2.T) @ f    * detJ * t * gw[iGauss] * gw[jGauss]

    return Ke, fe  # Returns stiffness matrix and nodal force vector

def quad9_Kmatrix(x_koordinater, y_koordinater, D_matrise, tykkelse, ekstern_last=None):
    tyk = tykkelse
    StivhetsMatrise = np.zeros((18, 18))
    LastVektor = np.zeros((18, 1))

    if ekstern_last is None:
        last = np.zeros((2, 1))
    else:
        last = np.array([ekstern_last]).T

    gaussPunkter = 5
    gaussPos, gaussVekt = gauss_points(gaussPunkter)

    for i in range(gaussPunkter):
        for j in range(gaussPunkter):
            xsi = gaussPos[i]
            eta = gaussPos[j]

            dN_dXsi = grad_xsi_ni_node(xsi, eta)
            dN_dEta = grad_eta_ni_node(xsi, eta)
            N = formfunksjoner_ni_node(xsi, eta)

            H_matrise = np.transpose([x_koordinater, y_koordinater])
            G_matrise = np.array([dN_dXsi, dN_dEta])

            Jacobian = G_matrise @ H_matrise
            invJacobian = np.linalg.inv(Jacobian)
            detJacobian = np.linalg.det(Jacobian)

            dN = invJacobian @ G_matrise
            dN_dx = dN[0]
            dN_dy = dN[1]

            B = np.array([utvid_0e(dN_dx), utvid_0s(dN_dy), utvid_0e(dN_dy) + utvid_0s(dN_dx)])
            N_utvidet = np.array([utvid_0e(N), utvid_0s(N)])

            StivhetsMatrise += B.T @ D_matrise @ B * detJacobian * tyk * gaussVekt[i] * gaussVekt[j]
            LastVektor += N_utvidet.T @ last * detJacobian * tyk * gaussVekt[i] * gaussVekt[j]

    return StivhetsMatrise, LastVektor if ekstern_last is not None else StivhetsMatrise

def formfunksjoner_ni_node(xsi, eta):
    return np.array([1/4*(xsi**2+xsi)*(eta**2+eta),1/4*(xsi**2-xsi)*(eta**2+eta),1/4*(xsi**2-xsi)*(eta**2-eta),
            1/4*(xsi**2+xsi)*(eta**2-eta),1/2*(1-xsi**2)*(eta**2+eta),1/2*(xsi**2-xsi)*(1-eta**2),
            1/2*(1-xsi**2)*(eta**2-eta),1/2*(xsi**2+xsi)*(1-eta**2),(1-xsi**2)*(1-eta**2)]) 

def grad_xsi_ni_node(xsi, eta):
    return np.array([eta*(1+eta)*(1+2*xsi)/4,eta*(1+eta)*(2*xsi-1)/4,(eta**2-eta)*(2*xsi-1)/4,
                    (eta**2-eta)*(2*xsi+1)/4,-xsi*eta*(eta+1),
                    xsi+(eta**2-1)/2-xsi*eta**2, -xsi*(eta**2-eta),
                    xsi+(1-eta**2)/2-xsi*eta**2,-2*xsi*(1-eta**2)])
def grad_eta_ni_node(xsi, eta):
    return np.array([xsi*(1+xsi)*(1+2*eta)/4,-xsi*(1-xsi)*(1+2*eta)/4,(xsi**2-xsi)*(2*eta-1)/4,
                    (xsi**2+xsi)*(2*eta-1)/4,eta+(1-xsi**2)/2-xsi**2*eta,
                    -eta*(xsi**2-xsi), eta+(xsi**2-1)/2-xsi**2*eta,-eta*(xsi**2+xsi),
                    -2*eta*(1-xsi**2)])

def utvid_0e(liste):
    
    z = np.zeros(len(liste)*2)
    j = 0
    for i in range(0,len(z),2):
        z[i] = liste[j]
        j += 1
    
    return z

def utvid_0s(liste):
    z = np.zeros(len(liste)*2)
    j = 0
    for i in range(0,len(z),2):
        z[i+1] = liste[j]
        j += 1
    
    return z



def beregn_elastisitetsmatrise(E, nu):
    return np.matrix([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1 - nu) / 2]
    ]) * E / (1 - nu ** 2)

def beregn_stivhet_4_node_quad(ex, ey, D, thickness):
    return quad4_Kmatrix(ex, ey, D, thickness)

def beregn_stivhet_9_node_quad(ex, ey, D, thickness):
    return quad9_Kmatrix(ex, ey, D, thickness)

tykkelse = 0.1
E_modul = 2.1e11
poissons_ratio = 0.3
D = beregn_elastisitetsmatrise(E_modul, poissons_ratio)

# 4 Node Quad
ex_q4 = np.array([1, -1, -1, 1])
ey_q4 = np.array([1, 1, -1, -1])
K4, feq4 = beregn_stivhet_4_node_quad(ex_q4, ey_q4, D, tykkelse)
print('Stivhet for 4 node quad:\n', K4, '\n')

# 9 Node Quad
ex_q9 = np.array([1, -1, -1, 1, 0, -1, 0, 1, 0])
ey_q9 = np.array([1, 1, -1, -1, 1, 0, -1, 0, 0])
K9, feq9 = beregn_stivhet_9_node_quad(ex_q9, ey_q9, D, tykkelse)
print('Stivhet for 9 node quad:\n', K9, '\n')