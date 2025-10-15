import numpy as np

def ShockSimilarity(K,gamma=1.4):
    return 2 * ( (gamma + 1)/4 + np.sqrt( (gamma + 1)/4)**2 + 1/K**2 )

def PrandtlMeyerSimilarity(K,gamma=1.4):
    return 2/(gamma * K) * ( (1 - (gamma-1)/2 * K)**((2*gamma)/(gamma-1)) - 1 ) 
