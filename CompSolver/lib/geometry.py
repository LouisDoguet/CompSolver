import numpy as np
import CompSolver.lib.newtonian as nw
import CompSolver.lib.similarity as sim
import CompSolver.lib.state as state

class Geometry:
    def __init__(self,name):
        self.name = name


class FlatPlate(Geometry):
    def __init__(self, name, theta, state0:state.State):
        super().__init__(name)

        self.theta = theta
        self.previous_state = state0

        if theta > 0:
            self.state = state0.AngleShock(theta=theta)
        else:
            self.state = state0.PrandtlMeyer(theta=theta)

        self.Cp = {}
        self.Cl = {}
        self.Cd = {}
        self.LD = {}

    def CpSimilarity(self):
        K = self.previous_state.M * self.theta
        if self.theta > 0 and self.previous_state.M > 1:
            self.Cp['Similarity'] = sim.ShockSimilarity(K,self.state.gamma)*np.deg2rad(self.theta)**2
        elif self.theta <= 0 and self.previous_state.M > 1:
            self.Cp['Similarity'] = sim.PrandtlMeyerSimilarity(K,self.state.gamma)*np.deg2rad(self.theta)**2
        
    def CpNewtonian(self):
        if self.theta > 0 and self.previous_state.M > 1:
            self.Cp['Newtonian'] = nw.CpNewtonian(self.theta)
        elif self.theta <= 0 and self.previous_state.M > 1:
            self.Cp['Newtonian'] = 0
    
    def CpExact(self):
        self.Cp['Exact'] = 2/(self.state.gamma*self.previous_state.M**2) * ((self.state.P/self.previous_state.P) - 1)

    def Solve(self):
        self.CpExact()
        self.CpSimilarity()
        self.CpNewtonian()