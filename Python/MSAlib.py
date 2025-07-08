import numpy as np

def Molarity2VolumeDensity(molarity):

    # molarity en moles/litro
    
    Av = 6.02214129 * 1e23		# Constante de Avogadro mol^(-1)
    return molarity * Av * 1e3  # Resultado en m^(-3)


# Construcción de los parámetros ξ en https://doi.org/10.1143/JPSJ.27.1415 (# Ec. 3.17)
def BuildXi(volumeDensity, diameter):

    # Diameter en Angstroms
    # VolumeDensity en m^(-3)
    
    xi = np.empty((0))
    
    eta = (np.pi / 6.0) * volumeDensity
    
    for i in range(4):
        xi = np.append(xi, np.sum( eta * np.power(diameter * 1e-10, i) ) )
        
    return xi

# Construcción de la función D(Γ) (llamada FGamma) y ∂D(Γ)/∂Γ (llamada dFGamma)
# Referencia -> https://doi.org/10.1080/00268977700101011 (Ec. 3)
def GammaFunction(gamma, xi, z, diameter, rho, e0, BETA, E):

    #valence = z          adimensional       
    #diameter = R         en angstroms
    #volumeDensity = rho  en m^(-3)
    
    
    #*******CAST ANGSTROMS TO METERS******
    
    R = np.copy(diameter)
    R *= 1e-10
    
    
    #***********Gamma Functions***********
    
    cTilde = (np.pi/2.0) / (1.0 - xi[3])
    
    zita = 1.0 + gamma * R
    
    chiNumerator = -cTilde * np.sum( (rho * R * z)*np.power(zita, -1.0) )
    
    chiDenominator = 1.0 + cTilde * np.sum( (rho * R*R*R) * np.power(zita, -1.0) )
    
    chi = chiNumerator / chiDenominator
    
    X = (z + R*R*chi) / zita
    
    DGamma = np.sum( rho * X * X )
    
    omega = (BETA * e0*e0 ) / (4.0 * np.pi * E)
    
    FGamma = gamma*gamma - np.pi * omega *DGamma
    
    #***********Gamma Functions Derivarives***********
    
    dchiNumerator =  cTilde * np.sum( (rho * R*R * z) * np.power(zita, -2.0) )
    
    dchiDenominator = -cTilde * np.sum( (rho * np.power(R, 4.0)) * np.power(zita, -2.0) )
    
    dchi = (dchiNumerator*chiDenominator - chiNumerator*dchiDenominator) / np.power(chiDenominator, 2.0)
    
    dX = (R*dchi*zita - z - R*R*chi) * R / np.power(zita, 2.0)
    
    dDGamma = np.sum( rho * 2.0 * X * dX )
    
    dFGamma = 2.0*gamma - np.pi * omega * dDGamma
    
    
    return FGamma, dFGamma


# Calcular el valor de Γ mediante Newton-Raphson
def GammaNewtonRaphson(xi, maxError, valence, diameter, volumeDensity, e0, BETA, E):
    
    gammaError = 1.0
    gammaSeed = 1.0
    
    # print("\n----------------------------")
    # print("CALCULANDO EL VALOR DE GAMMA")
    # print("----------------------------")
    
    while(gammaError > maxError):
        
        Fgamma, dFgamma = GammaFunction(gammaSeed, xi, valence, diameter, volumeDensity, e0, BETA, E)
        
        gamma = gammaSeed - (Fgamma / dFgamma)
        
        gammaError = np.abs((gammaSeed - gamma) / gamma)
        
        gammaSeed = gamma
    
    
    return gamma, gammaError


# Construcción de los coeficientes de esfera dura y la parte eléctrica.
# Se pueden consultar las referencias https://doi.org/10.1143/JPSJ.27.1415 y https://doi.org/10.1080/00268977700101011
# o la ecuación A34 de la tesis https://repositorioinstitucional.uaslp.mx/xmlui/handle/i/8296
def BuildCijMatrices(i, j, xi, gamma, z, rho, diameter, e0, E, BETA):
    
    #valence = z
    #diameter = R
    #volumeDensity = rho
    
    #**********CAST ANGSTROMS TO METERS**********
    R = np.copy(diameter)
    R *= 1e-10
    
    
    #**********BASIC FUNCTIONS**********
    
    Rij = (R[j] + R[i])/2.0
    
    lamij = np.abs(R[j] - R[i])/2.0
    
    
    omega = (BETA * e0*e0 ) / (4.0 * np.pi * E)
    
    eta = (np.pi / 6.0) * rho
    
    cTilde = (np.pi/2.0) / (1.0 - xi[3])
    
    zita = 1.0 + gamma * R
    
    chiNumerator = -cTilde * np.sum( (rho * R * z)*np.power(zita, -1.0) )
    
    chiDenominator = 1.0 + cTilde * np.sum( (rho * R*R*R) * np.power(zita, -1.0) )
                         
    chi = chiNumerator / chiDenominator
    
    X = (z + R*R*chi) / zita
    
    N = (R*chi - gamma*z) / zita
    
    
    #**********Cij FUNCTIONS**********
    
    if (R[j] > R[i]):
        beta0 = 2.0*omega * (z[i]*N[j] - X[i]*R[i]*chi + R[i]*R[i]*R[i]*chi*chi/3.0)
    else:
        beta0 = 2.0*omega * (z[j]*N[i] - X[j]*R[j]*chi + R[j]*R[j]*R[j]*chi*chi/3.0)
    
    
    alpha0 = omega*lamij*lamij * ((X[i]+X[j])*chi - Rij*Rij*chi*chi + N[i]*N[j])
    
    
    alpha1 = omega * ( (chi*chi/3.0) * (np.power(R[i], 3.0) + np.power(R[j],3.0)) - \
                    ( X[i] - X[j] ) * ( N[i] - N[j] ) - \
                    ( X[i]*X[i] + X[j]*X[j] ) * gamma - \
                    2.0 * Rij * N[i] * N[j] )
    
    
    alpha2 = omega * ( (X[i] + X[j])*chi + N[i]*N[j] - \
                       (chi*chi/2.0) * (R[i]*R[i] + R[j]*R[j]) )
    
    
    alpha3 = omega * chi * chi / 3.0
    
    
    ni = 1.0 / (1.0-xi[3])
    
    a = 3.0 * xi[2] * ni*ni
    
    b = 3.0 * ni*ni * ( xi[1]  + (3.0 * xi[2]*xi[2] * ni ) )
    
    c = 3.0 * ni*ni * ( xi[0] + (6.0 * xi[1]*xi[2] * ni) + \
                        (9.0 * np.power(xi[2], 3.0) * ni*ni) ) 
    
    
    q = ni + (a*R) + (b*R*R) + (c*R*R*R/3.0)

    v = (-a/2.0) - (b*R) - (c*R*R/2.0)
    
    
    #**********OTHER CONSTANTS**********
    
    Uij = z[i]*z[j] * e0*e0 / (4.0 * np.pi * E)
    
    pij = - (lamij*lamij/2.0) * (a + 2.0*b*Rij + c*Rij*Rij)
    
    qij = (q[i] + q[j]) / 2.0
    
    vij = (v[i] + v[j]) / 2.0
    
    w = (1.0/2.0) * np.sum(eta * q)
    
    
    #**********CREATION OF MATRICES**********
    cFConst = np.zeros((3,4))
    
    
    #**********FILLING THE MATRICES**********
    #filling for 0 <= s <= lamij

    if (R[j] > R[i]):
        cFConst[0, 1] = beta0 - q[i]
    else:
        cFConst[0, 1] = beta0 - q[j]
    
    #filling for lamij < s <= Rij
    cFConst[1, 0] = alpha0 - pij
    cFConst[1, 1] = alpha1 - qij
    cFConst[1, 2] = alpha2 - vij
    cFConst[1, 3] = alpha3 - w

    #filling for s > Rij
    cFConst[2, 0] = - BETA * Uij
        
    
    #**********CAST METERS TO ANGSTROMS**********
    cFConst[:,0] *= 1e10
    cFConst[:,2] *= 1e-10
    cFConst[:,3] *= 1e-30
    
    return cFConst

# Construcción de la función de correlación directa c(r) para MSA
def BuildCrMSA(r, z, rho, diameter, e0, E, BETA):

    nSpecies = len(diameter)

    Cij = np.zeros((nSpecies, nSpecies, len(r)))

    ######################################################################################

    maxError = 1e-8    # Umbral para el cálculo de Γ mediante Newton-Raphson

    # Construcción de los parámetros ξ
    xi = BuildXi(rho, diameter)
    
    ######################################################
    
    # Calculate gamma value 
    gamma, gammaError = GammaNewtonRaphson(xi, maxError, z, diameter, rho, e0, BETA, E)
    # print("gamma = ", gamma, "m^(-1)")
    # print("gammaError = ", gammaError, "\n")

    ######################################################################################

    for i in range(nSpecies):
        for j in range(nSpecies):

            Rij = (diameter[j] + diameter[i])/2.0
            lamij = np.abs(diameter[j] - diameter[i])/2.0

            cFConst = BuildCijMatrices(i, j, xi, gamma, z, rho, diameter, e0, E, BETA)

            for k in range(len(r)):

                if (r[k] < lamij):    #  0 <= r <= lamij

                    Cij[i,j,k] = cFConst[0, 1]

                elif (r[k] < Rij):    #  lamij < s <= Rij

                    Cij[i,j,k] = cFConst[1, 1] + \
                                 cFConst[1, 2] * r[k] + \
                                 cFConst[1, 3] * np.power(r[k], 3)                   

                    if (np.abs(lamij) > 0.0):
                    
                        Cij[i,j,k] += cFConst[1, 0] / r[k]

                else:                 #  r >= Rij
                    if (np.abs(Rij) > 0.0):
                        Cij[i,j,k] = cFConst[2, 0] / r[k]
    
    return Cij
    

# Construcción de la función de correlación directa C(k) para MSA
def BuildCkMSA(kVec, z, rho, diameter, e0, E, BETA, rmax):

    nSpecies = len(diameter)

    Cij = np.zeros((nSpecies, nSpecies, len(kVec)))

    ######################################################################################

    maxError = 1e-8    # Umbral para el cálculo de Γ mediante Newton-Raphson

    # Construcción de los parámetros ξ
    xi = BuildXi(rho, diameter)
    
    ######################################################
    
    # Calculate gamma value 
    gamma, gammaError = GammaNewtonRaphson(xi, maxError, z, diameter, rho, e0, BETA, E)
    # print("gamma = ", gamma, "m^(-1)")
    # print("gammaError = ", gammaError, "\n")

    ######################################################################################

    for i in range(nSpecies):
        for j in range(nSpecies):

            Rij = (diameter[j] + diameter[i])/2.0
            lamij = np.abs(diameter[j] - diameter[i])/2.0

            cFConst = BuildCijMatrices(i, j, xi, gamma, z, rho, diameter, e0, E, BETA)

            for k in range(len(kVec)):

                if (kVec[k] < 0.002):

                    Cij[i,j,k] = cFConst[1, 1] * FT_TaylorSeries_at_k0(kVec[k], 0, lamij, Rij) + \
                                 cFConst[1, 2] * FT_TaylorSeries_at_k0(kVec[k], 1, lamij, Rij) + \
                                 cFConst[1, 3] * FT_TaylorSeries_at_k0(kVec[k], 3, lamij, Rij) + \
                                 cFConst[0, 1] * FT_TaylorSeries_at_k0(kVec[k], 0, 0, lamij)                                 
                    
                    if (np.abs(lamij) > 1e-10):
                        Cij[i,j,k] += cFConst[1, 0] * FT_TaylorSeries_at_k0(kVec[k], -1, lamij, Rij)

                    if (np.abs(Rij) > 0.0):
                        Cij[i,j,k] += cFConst[2, 0] * FT_TaylorSeries_at_k0(kVec[k], -1, Rij, rmax)

                else:

                    Cij[i,j,k] = cFConst[1, 1] * FT(kVec[k], 0, lamij, Rij) + \
                                 cFConst[1, 2] * FT(kVec[k], 1, lamij, Rij) + \
                                 cFConst[1, 3] * FT(kVec[k], 3, lamij, Rij) + \
                                 cFConst[0, 1] * FT(kVec[k], 0, 0, lamij)
                    
                    if (np.abs(lamij) > 1e-10):
                        Cij[i,j,k] += cFConst[1, 0] * FT(kVec[k], -1, lamij, Rij)

                    if (np.abs(Rij) > 0.0):
                        Cij[i,j,k] += cFConst[2, 0] * FT(kVec[k], -1, Rij, rmax)


    return (4*np.pi) * Cij


def BuildPhik(kVec, z, rho, diameter, e0, E, BETA, rmax):

    nSpecies = len(diameter)

    Phi_ij = np.zeros((nSpecies, nSpecies, len(kVec)))

    for i in range(nSpecies):
        for j in range(nSpecies):

            Uij = z[i]*z[j] * e0*e0 / (4.0 * np.pi * E)
            phi_ij = BETA * Uij * 1e10

            for k in range(len(kVec)):
                    
                if (kVec[k] < 0.002):

                    Phi_ij[i,j,k] = phi_ij * FT_TaylorSeries_at_k0(kVec[k], -1, 0, rmax)
                
                else:

                    Phi_ij[i,j,k] = phi_ij * FT(kVec[k], -1, 0, rmax)


    return (4*np.pi) * Phi_ij


def BuildQk(kVec, z, rho, diameter, e0, E, BETA, rmax):

    nSpecies = len(diameter)

    Qij = np.zeros((nSpecies, nSpecies, len(kVec)))

    kappa = np.sqrt(BETA * e0*e0 * np.sum(rho*z*z) / E)
    kappa *= 1e-10    # cast m^(-1) to A^(-1)

    if (np.abs(kappa) < 1e-12):
        return Qij


    for i in range(nSpecies):
        for j in range(nSpecies):

            Uij = z[i]*z[j] * e0*e0 / (4.0 * np.pi * E)
            qij = - BETA * Uij * 1e10

            for k in range(len(kVec)):
                    
                if (kVec[k] < 0.002):

                    Qij[i,j,k] = qij * FTexp_TaylorSeries_at_k0(kVec[k], kappa, 0, rmax)
                
                else:

                    Qij[i,j,k] = qij * FTexp(kVec[k], kappa, 0, rmax)

    return (4*np.pi) * Qij

# Construcción del factor de estructura S(k) empleando la Ec. 10 de https://doi.org/10.1063/1.441426
def BuildSkMSA(kVec, Cij, z, rho, diameter, e0, E, BETA, rmax):

    # Cij = BuildCkMSA(kVec, z, rho, diameter, e0, E, BETA, rmax)
    phi = BuildPhik(kVec, z, rho, diameter, e0, E, BETA, rmax)
    Q   = BuildQk(kVec, z, rho, diameter, e0, E, BETA, rmax)

    F = Cij + phi

    nSpecies = len(rho)
    N = len(Cij[0,0,:])

    I = np.identity(len(rho))

    rhoMat = np.zeros((nSpecies, nSpecies))
    for i in range(nSpecies):
        for j in range(nSpecies):
            rhoMat[i, j] = np.sqrt(rho[i]*rho[j]) * 1e-30


    Sij_k = np.zeros(Cij.shape)

    for k in range(N):
        
        Fij = rhoMat * F[:,:,k]
        Qij = rhoMat * Q[:,:,k]

        A = np.matmul(I+Qij, Fij)
        B = np.matmul(A, Qij+Fij)
        C = np.matmul(Qij, Fij) + B
        D = np.linalg.inv(I - A)
        G = Qij + Fij + I
        
        Sij_k[:,:,k] = G + np.matmul(D, C)

    return Sij_k


def BuildHkMSA(Sij, rho):

    nSpecies = len(rho)
    N = len(Sij[0,0,:])

    I = np.identity(len(rho))

    rhoMat = np.zeros((nSpecies, nSpecies))
    for i in range(nSpecies):
        for j in range(nSpecies):
            rhoMat[i, j] = np.sqrt(rho[i]*rho[j]) * 1e-30

    Hij_k = np.zeros(Sij.shape)

    for k in range(N):
        
        Hij_k[:,:,k] = (1.0/rhoMat) * (Sij[:,:,k] - I)

    return Hij_k

# Construcción del factor de estructura S(k) empleando la relación S(k) = [I - rho*C(k)]^(-1)
def BuildSk(nSpecies, rhoVec, Cij_k):

    N = len(Cij_k[0,0,:])

    I = np.identity(nSpecies)

    Sij_k = np.zeros(Cij_k.shape)

    rhoMat = np.zeros((nSpecies, nSpecies))
    for i in range(nSpecies):
        for j in range(nSpecies):
            rhoMat[i, j] = np.sqrt(rhoVec[i]*rhoVec[j]) * 1e-30

    for i in range(N):
        Cij_submat = Cij_k[:,:,i]

        Sij_submat = np.linalg.inv(I - rhoMat*Cij_submat)

        Sij_k[:,:,i] = Sij_submat

    return Sij_k



def Ck_PYHS(kVec, phi, sigma):

    result = np.zeros(kVec.shape)

    eta = phi

    alpha = np.power(1+2*eta, 2) / np.power(1-eta, 4)
    beta = -6*eta * np.power(1+0.5*eta, 2) / np.power(1-eta, 4)
    gamma = 0.5*eta*alpha

    # print(f"==================")
    # print(f"q = {alpha}")
    # print(f"v = {beta}")
    # print(f"w = {gamma}")

    for ind, k in enumerate(kVec):

        ks = k*sigma

        AA = alpha * np.power(ks,3)*(np.sin(ks) - ks*np.cos(ks))
        BB = beta * np.power(ks,2)*(2*ks*np.sin(ks) - (np.power(ks,2)-2)*np.cos(ks) - 2)
        CC = gamma * ((4*np.power(ks,3) - 24*ks)*np.sin(ks) - (np.power(ks,4) - 12*np.power(ks,2) + 24)*np.cos(ks) + 24)

        result[ind] = AA + BB + CC
        result[ind] *= np.power(sigma, 3)
        
        result[ind] *= -1/np.power(ks,6)

    return result


def BuildCompressibility(r, Cr, rho):

    nSpecies = len(rho)

    rhoTotal = np.sum(rho)

    compressibility = 0.0

    for i in range(nSpecies):
        for j in range(nSpecies):

            rhoFactor = 1e-30 * rho[i]*rho[j] / rhoTotal

            temp = 0.0

            for k in range(1,len(r)):

                temp += (Cr[i, j, k] + Cr[i, j, k-1]) * (r[k]-r[k-1]) / 2.0

            compressibility += rhoFactor * temp
    
    return 1 - 4*np.pi * compressibility



def BuildNumericHr(k, r, Hk, rhoVec):

    nSpecies = len(rhoVec)

    Hr = np.zeros((nSpecies, nSpecies, len(r)))

    rhoMat = np.zeros((nSpecies, nSpecies))
    for i in range(nSpecies):
        for j in range(nSpecies):
            rhoMat[i, j] = np.sqrt(rhoVec[i]*rhoVec[j]) * 1e-30

    for i in range(nSpecies):
        for j in range(nSpecies):
            for o in range(len(r)):

                integrand = k * Hk[i,j,:] * np.sin(r[o]*k) / r[o]

                Hr[i, j, o] = 4.0 * np.pi * trapz(k, integrand) / (np.power(2*np.pi,3))

    return Hr



# ################################################################################
# FÓRMULAS ANALÍTICAS Y MÉTODOS NUMÉRICOS
# ################################################################################

def FT(k, power, A, B):

    # FT es la transformada de fourier de la forma
    # ∫_{A}^{B} [sin(k r) r^(n+1) / k] dr
    # n es el parámetro power que se da como input en FT

    if (power == -1):
        result = ( (-np.cos(k*B)) - \
                   (-np.cos(k*A)) ) * np.power(k, -2)

    elif (power == 0):
        result = ( (np.sin(k*B) - k*B*np.cos(k*B)) - \
                   (np.sin(k*A) - k*A*np.cos(k*A)) ) * np.power(k, -3)

    elif (power == 1):
        result = ( ((2 - np.power(k*B, 2))*np.cos(k*B) + 2*k*B*np.sin(k*B)) - \
                   ((2 - np.power(k*A, 2))*np.cos(k*A) + 2*k*A*np.sin(k*A)) )  * np.power(k, -4)

    elif (power == 3):

        result = ( (4*k*B * (np.power(k*B, 2) - 6) * np.sin(k*B) - (np.power(k*B, 4) - 12*np.power(k*B, 2) + 24)*np.cos(k*B)) - 
                   (4*k*A * (np.power(k*A, 2) - 6) * np.sin(k*A) - (np.power(k*A, 4) - 12*np.power(k*A, 2) + 24)*np.cos(k*A)) ) * np.power(k, -6)

    return result


def FTexp(k, kappa, A, B):

    # FT es la transformada de fourier de la forma
    # ∫_{A}^{B} [sin(k r) exp(-kappa r) / k] dr

    result = - np.exp(-kappa*B) * (kappa*np.sin(k*B) + k*np.cos(k*B)) / (k*(kappa*kappa + k*k)) + \
               np.exp(-kappa*A) * (kappa*np.sin(k*A) + k*np.cos(k*A)) / (k*(kappa*kappa + k*k))

    return result


def FT_TaylorSeries_at_k0(k, power, A, B):

    # Expansión, alrededor de k=0, en series de Taylor de las transformadas de Fourier
    # ∫_{A}^{B} [sin(k r) r^(n+1) / k] dr
    # n es el parámetro power que se da como input en FT

    if (power == -1):
    
        result = (B**2 - A**2)/2.0 - (k**2)*(B**4 - A**4)/24.0 + (k**4)*(B**6 - A**6)/720.0

    elif (power == 0):

        result = (B**3 - A**3)/3.0 - (k**2)*(B**5 - A**5)/30.0 + (k**4)*(B**7 - A**7)/840.0

    elif (power == 1):

        result = (B**4 - A**4)/4.0 - (k**2)*(B**6 - A**6)/36.0 + (k**4)*(B**8 - A**8)/960.0

    elif (power == 3):

        result = (B**6 - A**6)/6.0 - (k**2)*(B**8 - A**8)/48.0 + (k**4)*(B**10 - A**10)/1200.0
        
    return result


def FTexp_TaylorSeries_at_k0(k, kappa, A, B):

    # Expansión, alrededor de k=0, en series de Taylor de las transformadas de Fourier
    # ∫_{A}^{B} [sin(k r) exp(-kappa r) / k] dr

    KB  = kappa*B
    KB2 = np.power(kappa*B, 2)
    KB3 = np.power(kappa*B, 3)
    KB4 = np.power(kappa*B, 4)
    KB5 = np.power(kappa*B, 5)
 
    KA  = kappa*A
    KA2 = np.power(kappa*A, 2)
    KA3 = np.power(kappa*A, 3)
    KA4 = np.power(kappa*A, 4)
    KA5 = np.power(kappa*A, 5)
 
    AB = -np.exp(-KB) * (KB + 1) / np.power(kappa,2)
    AA = -np.exp(-KA) * (KA + 1) / np.power(kappa,2)

    BB = k*k * np.exp(-KB) * (KB3 + 3*KB2 + 6*KB + 6) / (6*np.power(kappa,4))
    BA = k*k * np.exp(-KA) * (KA3 + 3*KA2 + 6*KA + 6) / (6*np.power(kappa,4))

    CB = k*k*k*k * np.exp(-KB) * (KB5 + 5*KB4 + 20*KB3 + 60*KB2 + 120*KB + 120) / (120*np.power(kappa,6))
    CA = k*k*k*k * np.exp(-KA) * (KA5 + 5*KA4 + 20*KA3 + 60*KA2 + 120*KA + 120) / (120*np.power(kappa,6))

    result = (AB-AA) + (BB-BA) - (CB-CA)

    return result

# Integral ∫_{A}^{B} f(r) dr por el método del trapecio
# r debe estar ordenado de menor a mayor
def trapz(r, fr):

    # Check if r and fr have the same length

    if (len(r) != len(fr)):
        print("r y f(r) no tienen la misma longitud")
        return 0

    integral = 0.0

    for i in range(1, len(r)):
        integral += (fr[i] + fr[i-1])*(r[i] - r[i-1]) / 2.0

    return integral


# ################################################################################
# Fórmulas para esfera dura (HS)
# ################################################################################

# ESTA FUNCIÓN SE PUEDE DESCARTAR O GUARDAR PARA FUTUROS CÁLCULOS 
# DONDE SOLO SE REQUIERA ESFERA DURA
# Construcción de los coeficientes de esfera dura solamente
def BuildCijMatricesHS(i, j, xi, rho, diameter):

    #valence = z
    #diameter = R
    #volumeDensity = rho
    
    #**********CAST ANGSTROMS TO METERS**********
    R = np.copy(diameter)
    R *= 1e-10
    
    #**********BASIC FUNCTIONS**********
    
    Rij = (R[j] + R[i])/2.0
    
    lamij = np.abs(R[j] - R[i])/2.0

    eta = (np.pi / 6.0) * rho

    #**********Cij FUNCTIONS**********
    
    ni = 1.0 / (1.0-xi[3])

    a = 3.0 * xi[2] * ni*ni
    
    b = 3.0 * ni*ni * ( xi[1]  + (3.0 * xi[2]*xi[2] * ni ) )
    
    c = 3.0 * ni*ni * ( xi[0] + (6.0 * xi[1]*xi[2] * ni) + \
                        (9.0 * np.power(xi[2], 3.0) * ni*ni) ) 
    
    
    q = ni + (a*R) + (b*R*R) + (c*R*R*R/3.0)

    v = (-a/2.0) - (b*R) - (c*R*R/2.0)
    
    
    #**********OTHER CONSTANTS**********
    
    pij = - (lamij*lamij/2.0) * (a + 2.0*b*Rij + c*Rij*Rij)
    
    qij = (q[i] + q[j]) / 2.0
    
    vij = (v[i] + v[j]) / 2.0
    
    w = (1.0/2.0) * np.sum(eta * q)
    
    
    #**********CREATION OF MATRICES**********
    cFConst = np.zeros((3,4))
    
    
    #**********FILLING THE MATRICES**********
    #filling for 0 <= s <= lamij

    if (R[j] > R[i]):
        cFConst[0, 1] = - q[i]
    else:
        cFConst[0, 1] = - q[j]
    
    #filling for lamij < s <= Rij
    cFConst[1, 0] = - pij
    cFConst[1, 1] = - qij
    cFConst[1, 2] = - vij
    cFConst[1, 3] = - w
    
    
    #**********CAST METERS TO ANGSTROMS**********
    cFConst[:,0] *= 1e10
    cFConst[:,2] *= 1e-10
    cFConst[:,3] *= 1e-30
    
    return cFConst

def BuildCrHS(r, rho, diameter):

    nSpecies = len(diameter)

    xi = BuildXi(rho, diameter)

    Cij = np.zeros((nSpecies, nSpecies, len(r)))

    for i in range(nSpecies):
        for j in range(nSpecies):

            Rij = (diameter[j] + diameter[i])/2.0
            lamij = np.abs(diameter[j] - diameter[i])/2.0

            cFConst = BuildCijMatricesHS(i, j, xi, rho, diameter)

            for k in range(len(r)):

                if (r[k] < lamij):

                    Cij[i,j,k] = cFConst[0, 1]

                elif (r[k] < Rij):

                    Cij[i,j,k] = cFConst[1, 1] + \
                                cFConst[1, 2] * r[k] + \
                                cFConst[1, 3] * np.power(r[k], 3)

                    if (np.abs(lamij) > 1e-10):

                        Cij[i,j,k] += cFConst[1, 0] / r[k]
                                     
                else:
                    Cij[i,j,k] = 0
    
    return Cij

def BuildCkHS(kVec, rho, diameter):

    nSpecies = len(diameter)

    xi = BuildXi(rho, diameter)

    Cij = np.zeros((nSpecies, nSpecies, len(kVec)))

    for i in range(nSpecies):
        for j in range(nSpecies):

            Rij = (diameter[j] + diameter[i])/2.0
            lamij = np.abs(diameter[j] - diameter[i])/2.0

            cFConst = BuildCijMatricesHS(i, j, xi, rho, diameter)

            for k in range(len(kVec)):

                if (kVec[k] < 0.002):

                    Cij[i,j,k] = cFConst[1, 1] * FT_TaylorSeries_at_k0(kVec[k], 0, lamij, Rij) + \
                                 cFConst[1, 2] * FT_TaylorSeries_at_k0(kVec[k], 1, lamij, Rij) + \
                                 cFConst[1, 3] * FT_TaylorSeries_at_k0(kVec[k], 3, lamij, Rij) + \
                                 cFConst[0, 1] * FT_TaylorSeries_at_k0(kVec[k], 0, 0, lamij)
                    
                    if (np.abs(lamij) > 1e-10):

                        Cij[i,j,k] += cFConst[1, 0] * FT_TaylorSeries_at_k0(kVec[k], -1, lamij, Rij)

                else:

                    Cij[i,j,k] = cFConst[1, 1] * FT(kVec[k], 0, lamij, Rij) + \
                                 cFConst[1, 2] * FT(kVec[k], 1, lamij, Rij) + \
                                 cFConst[1, 3] * FT(kVec[k], 3, lamij, Rij) + \
                                 cFConst[0, 1] * FT(kVec[k], 0, 0, lamij)
                    
                    if (np.abs(lamij) > 1e-10):

                        Cij[i,j,k] += cFConst[1, 0] * FT(kVec[k], -1, lamij, Rij)

    return (4*np.pi) * Cij