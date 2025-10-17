import numpy as np

def ShockSimilarity(K, gamma):
    term = ((gamma + 1) / 4.0)
    return 2.0 * (term + np.sqrt(term**2 + 1.0 / K**2))

def PrandtlMeyerSimilarity(K, gamma):
    factor = 2.0 / (gamma * K**2)
    base = 1.0 - ((gamma - 1.0) / 2.0) * K
    return factor * (base**((2.0 * gamma) / (gamma - 1.0)) - 1.0)
