import CompSolver as cs

s0 = cs.State(M=10)
print(s0.CpExpSimilarity(18))
s1 = s0.PrandtlMeyer(theta=18)
print(s1.Cp(s0))
