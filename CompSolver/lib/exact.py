import numpy as np

def CpExact(M0, p0, p1, gamma):
    return 2/(gamma*M0**2) * ((p1/p0) - 1)