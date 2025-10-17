import numpy as np
import CompSolver.lib.newtonian as nw
import CompSolver.lib.similarity as sim
import CompSolver.lib.exact as ex
import CompSolver.lib.state as state


class Geometry:
    def __init__(self):
        pass


class FlatPlate(Geometry):
    def __init__(self, theta, state0:state.State):
        super().__init__()

        self.theta = theta
        self.previous_state = state0

        self.shock_state = state0.AngleShock(theta=np.abs(theta))
        self.expansion_state = state0.PrandtlMeyer(theta=np.abs(theta))

        self.Cp = {}
        self.Cl = {}
        self.Cd = {}
        self.LD = {}

    def CpSimilarity(self):

        K = self.previous_state.M * np.abs(np.deg2rad(self.theta))

        Cp_shock = sim.ShockSimilarity(K,self.shock_state.gamma) * np.deg2rad(self.theta)**2
        Cp_expansion = sim.PrandtlMeyerSimilarity(K,self.expansion_state.gamma) * np.deg2rad(self.theta)**2

        if self.theta < 0:
            Cpp = Cp_shock
            Cpm = -Cp_expansion
        else:
            Cpp = Cp_expansion
            Cpm = -Cp_shock

        self.Cp['Similarity'] = [Cpp, Cpm]

        
    def CpNewtonian(self):

        Cp_shock = nw.CpNewtonian(self.theta)
        Cp_expansion = 0

        if self.theta < 0:
            Cpp = Cp_shock
            Cpm = -Cp_expansion
        else:
            Cpp = Cp_expansion
            Cpm = -Cp_shock

        self.Cp['Newtonian'] = [Cpp, Cpm]
    
    def CpExact(self):

        Cp_shock = np.abs( ex.CpExact(self.previous_state.M, self.previous_state.P, self.shock_state.P, self.shock_state.gamma) )
        Cp_expansion = np.abs( ex.CpExact(self.previous_state.M, self.previous_state.P, self.expansion_state.P, self.shock_state.gamma) )

        if self.theta < 0:
            Cpp = Cp_shock
            Cpm = -Cp_expansion
        else:
            Cpp = Cp_expansion
            Cpm = -Cp_shock

        self.Cp['Exact'] = [Cpp, Cpm]


    def setCd(self):
        self.Cd['Exact'] = (self.Cp['Exact'][0] - self.Cp['Exact'][1]) * np.sin(np.deg2rad(self.theta))
        self.Cd['Similarity'] = (self.Cp['Similarity'][0] - self.Cp['Similarity'][1]) * np.sin(np.deg2rad(self.theta))
        self.Cd['Newtonian'] = (self.Cp['Newtonian'][0] - self.Cp['Newtonian'][1]) * np.sin(np.deg2rad(self.theta))

    def setCl(self):
        self.Cl['Exact'] = (self.Cp['Exact'][0] - self.Cp['Exact'][1]) * np.cos(np.deg2rad(self.theta))
        self.Cl['Similarity'] = (self.Cp['Similarity'][0] - self.Cp['Similarity'][1]) * np.cos(np.deg2rad(self.theta))
        self.Cl['Newtonian'] = (self.Cp['Newtonian'][0] - self.Cp['Newtonian'][1]) * np.cos(np.deg2rad(self.theta))

    def setLD(self):
        self.LD['Exact'] = np.abs(self.Cl['Exact'] / self.Cd['Exact'])
        self.LD['Similarity'] = np.abs(self.Cl['Similarity'] / self.Cd['Similarity'])
        self.LD['Newtonian'] = np.abs(self.Cl['Newtonian'] / self.Cd['Newtonian'])

    def Solve(self):
        self.CpExact()
        self.CpSimilarity()
        self.CpNewtonian()
        self.setCd()
        self.setCl()
        self.setLD()

    def __str__(self):
        s =  "      ======================================\n"
        s += "     ¦    Exact   ¦  Newtonian ¦ Similarity ¦\n"
        s += "     ¦     --     ¦     --     ¦     --     ¦\n"
        s += f" Cp+ ¦ {self.Cp['Exact'][0]:>10.5f} ¦ {self.Cp['Newtonian'][0]:>10.5f} ¦ {self.Cp['Similarity'][0]:>10.5f} ¦\n"
        s += f" Cp- ¦ {self.Cp['Exact'][1]:>10.5f} ¦ {self.Cp['Newtonian'][1]:>10.5f} ¦ {self.Cp['Similarity'][1]:>10.5f} ¦\n"
        s += "     ¦     --     ¦     --     ¦     --     ¦\n"
        s += f" Cd  ¦ {self.Cd['Exact']:>10.5f} ¦ {self.Cd['Newtonian']:>10.5f} ¦ {self.Cd['Similarity']:>10.5f} ¦\n"
        s += f" Cl  ¦ {self.Cl['Exact']:>10.5f} ¦ {self.Cl['Newtonian']:>10.5f} ¦ {self.Cl['Similarity']:>10.5f} ¦\n"
        s += "     ¦     --     ¦     --     ¦     --     ¦\n"
        s += f" LD  ¦ {self.LD['Exact']:>10.5f} ¦ {self.LD['Newtonian']:>10.5f} ¦ {self.LD['Similarity']:>10.5f} ¦\n"
        s += "      ======================================\n"
        return s
