import CompSolver as cs

s = cs.State(M=2)
ss = s.PrandtlMeyer(23.38,bool_print = True)

print(s)
print(ss)

