import CompSolver as cs
import numpy as np

s0 = cs.State(M=10)

PLATE = cs.FlatPlate('fp1',-18,s0)
PLATE.Solve()

print(PLATE)
