import numpy as np
import matplotlib.pyplot as plt

# CONSTANTS DEFINITION
T_ = 288.15
Mair = 0.028964
R = 8.314/Mair
gamma_ = 1.4
P_ = 101315

def MachNumber(g=gamma_, T=T_)->float:
    return np.sqrt( g * R * T )

def P_ratio(M1:float, g=gamma_)->float:
    return (1 + (2*g)/(g+1) * (M1**2 - 1) )

def U_ratio(M1:float, g=gamma_)->float:
    return ( (2 + (g-1)*M1**2)/((g+1)*M1**2)  )

def RHO_ratio(M1:float, g=gamma_)->float:
    return 1/U_ratio(M1,g)

def T_ratio(M1:float, g=gamma_)->float:
    return P_ratio(M1, g) * U_ratio(M1, g)

def M_ratio(M1:float, g=gamma_)->float:
    M22 = (1 + (g-1)/2 * M1**2) / (g*M1**2 - (g-1)/2)
    return np.sqrt(M22)/M1

def M2U(M:float, g=gamma_, T=T_)->float:
    return M * MachNumber(g,T)

def U2M(U:float, g=gamma_, T=T_)->float:
    return U / MachNumber(g,T)

def T0_ratio(M,g=gamma_):
    return 1 + (g-1)/2 * M**2

def P0_ratio(M,g=gamma_):
    return T0_ratio(M,g) ** (g/(g-1))

def RHO0_ratio(M,g=gamma_):
    return T0_ratio(M,g) ** (1/(g-1))

def getTheta(beta, M, g=gamma_):
    beta = np.deg2rad(beta)
    theta = np.arctan( 2 * 1/np.tan(beta) * (M**2 * np.sin(beta)**2 - 1 ) / ( M**2 * (g + np.cos(2*beta)) + 2 ) )
    return np.rad2deg(theta)

def getBeta(theta, M, bool_weak=True, g=gamma_):
    theta = np.deg2rad(theta)
    if not bool_weak:
        B_init = np.deg2rad(80)
    else:
        B_init = np.deg2rad(40)

    a = (1 + (g-1)/2 * M**2) * np.tan(theta) 
    b = (M**2 - 1)
    c = (1 + (g+1)/2 * M**2) * np.tan(theta)
    alg = algoNewton(a, b, c, np.tan(B_init))
    res = alg[0]
    err = alg[1]
    
    while err > 1e-6:
        alg = algoNewton(a, b, c, np.tan(res))
        res = alg[0]
        err = alg[1]

    return np.rad2deg(res)

def algoNewton(a, b, c, x):
    res = (2*a*x**3 - b*x**2 - 1)/( 3*a*x**2 - 2*b*x + c )
    theta1 = np.arctan(res)
    theta0 = np.arctan(x)
    return [theta1, np.abs(theta1 - theta0)]


def PrandtlMeyerFunction(M,g=gamma_):
    A = np.sqrt((g+1)/(g-1))
    B = np.arctan( np.sqrt((g-1)*(M**2-1) / (g+1)) )
    C = np.arctan( np.sqrt(M**2 - 1) )
    return A * B - C

def InvPMF(target, g=gamma_, tol=1e-6):
    M_low, M_high = 1.0, 50.0
    while M_high - M_low > tol:
        M_mid = 0.5*(M_low + M_high)
        if PrandtlMeyerFunction(M_mid, g) < target:
            M_low = M_mid
        else:
            M_high = M_mid
    return 0.5*(M_low + M_high)

def T_ratio_PM(M1,M2,g=gamma_):
    return T0_ratio(M1,g)/T0_ratio(M2,g)

def P_ratio_PM(M1,M2,g=gamma_):
    return T_ratio_PM(M1,M2,g) ** (g/(g-1))

def RHO_ratio_PM(M1,M2,g=gamma_):
    return T_ratio_PM(M1,M2,g) ** (1/(g-1))

def isPerfectGas(state):
    return state.P/(state.rho * R * state.T)
