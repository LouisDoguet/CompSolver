import CompSolver as cs

inf_state = cs.State(M=10)
flat_plate = cs.FlatPlate(18, inf_state)
flat_plate.Solve()
print(flat_plate)