import numpy as np
import matplotlib.pyplot as plt

import CompSolver.lib.thermo as th
import CompSolver.lib.utils as ut
import CompSolver.lib.newtonian as nw

ID = 0

class State:
    """
    Classe storing state information and generating modified states (see methods)
    """
    
    # ID generator
    _ID = 0

    def __init__(self,
                 U=None,
                 M=None,
                 rho=None,
                 P=None, 
                 T=None,
                 gamma=None,
                 ID=None,
                 bool_total=False):
        """
        Generate a State object.
        
        Select only one velocity (the other is computed from the aerothermal data)

        If a parameter is not set by the user, Patm,T288,gamma1.4 is assumed

        :param U: Velocity (m/s)
        :param M: Velocity (Mach)
        :param rho: Density (kg/m3)
        :param P: Pressure (Pa)
        :param T: Temperature (K)
        :param gamma: Cp/Cv (default 1.4)
        :param ID: Manual setup of the ID
        :param bool_total: [DEV] stops recursivity

        """


        # ID DEFINITION (automatic increase)
        self.ID = ID
        if not ID:
            self.ID = type(self)._ID
            type(self)._ID += 1
       
        # SET GAMMA
        self.gamma = gamma
        if not gamma:
            self.gamma = th.gamma_
        
        # SET TEMPERATURE
        self.T = T
        if not T:
            self.T = th.T_

        # SET VELOCITIES
        self.a = th.MachNumber(self.gamma,self.T)
        self.U = U
        self.M = M

        if not (U or M):
            raise ValueError("Init Error")

        if not U:
            self.M = M
            self.U = th.M2U(M)

        if not M:
            self.U = U
            self.M = th.U2M(U)
        
        # SET THERMO VARIABLES
        
        self.P = P
        if not P:
            self.P = th.P_
        
        self.rho = rho
        if not rho:
            self.rho = self.P / (th.R * self.T)

        if not bool_total:
            M0 = 1
            T0r = th.T0_ratio(self.M,self.gamma)
            P0r = th.P0_ratio(self.M,self.gamma)
            Rho0r = th.RHO0_ratio(self.M,self.gamma)
            self.TotalState = State(
                    M=M0,
                    rho = Rho0r * self.rho,
                    P = P0r * self.P,
                    T = T0r * self.T,
                    ID = -1,
                    bool_total = True,
                    )



    def __str__(self)->str:
        if self.ID != -1:
            str = f"#### STATE {self.ID} : #### \n"
            str+= f" - M = {self.M} \n"
            str+= f" - U = {self.U} m/s \n"
            str+= f" - P = {self.P} Pa \n"
            str+= f" - rho = {self.rho} kg/m3 \n"
            str+= f" - T = {self.T} K \n"
            str+= f"{self.TotalState}"
            str+= f" - P0/P = {self.TotalState.P / self.P} \n"
            str+= f" - T0/T = {self.TotalState.T / self.T} \n"
            str+= "-------------------- \n"
            str+= f" - Perfect : {th.isPerfectGas(self)}\n"
        else:
            str =  "~~~ Total State : ~~~\n"
            str+= f" - P0 = {self.P} Pa\n"
            str+= f" - T0 = {self.T} K\n"
            str+= f" - Rho0 = {self.rho} kg/m3\n"
        return str

    def NormalShock(self,bool_print=False):
        if self.M < 1:
            raise ValueError("M<1 : No shock")
        Ur = th.U_ratio(self.M,self.gamma)
        Mr = th.M_ratio(self.M,self.gamma)
        Rhor = th.RHO_ratio(self.M,self.gamma)
        Pr = th.P_ratio(self.M,self.gamma)
        Tr = th.T_ratio(self.M,self.gamma)
        outputState = State(
                    U = Ur * self.U,
                    M = Mr * self.M,
                    rho = Rhor * self.rho,
                    P = Pr * self.P,
                    T = Tr * self.T,
                    gamma = self.gamma
                )
        if bool_print:
            str = "-- Shock ratios : -- \n"
            str+= "-------------------- \n"
            str+=f" - U2/U1 = {Ur} \n"
            str+=f" - M2/M1 = {Mr} \n"
            str+= "-------------------- \n"
            str+=f" - Rho2/Rho1 = {Rhor} \n"
            str+=f" - P2/P1 = {Pr} \n"
            str+=f" - T2/T1 = {Tr} \n"
            str+= "-------------------- \n"
            str+=f" - P02/P01 = {outputState.TotalState.P / self.TotalState.P} \n"
            str+=f" - T02/T01 = {outputState.TotalState.T / self.TotalState.T} \n"
            str+=f" - Rho02/Rho01 = {outputState.TotalState.rho / self.TotalState.rho} \n"
            str+= "-------------------- \n"
            print(str)

        return outputState   

    def CpShockSimilarity(self, angle):
        K = self.M * np.deg2rad(angle)
        return nw.ShockSimilarity(K,self.gamma)*np.deg2rad(angle)**2

    def CpExpSimilarity(self, angle):
        K = self.M * np.deg2rad(angle)
        return nw.PrandtlMeyerSimilarity(K,self.gamma)*np.deg2rad(angle)**2

    def Cp(self, state0):
        return 2/(self.gamma*state0.M**2) * ((self.P/state0.P) - 1)
    
    def AngleShock(self, beta=None, theta=None, bool_weak=True, bool_print=False):

        if not theta:
            theta = th.getTheta(beta, self.M)
            print(f"TBM -> Theta : {theta} deg \n" +
                  f"       Beta  : {beta} deg \n" +
                  f"       Only 1 solution (weak/strong bool ignored) \n")
        if not beta:
            beta = th.getBeta(theta, self.M, bool_weak)
            print(f"       Theta : {theta} deg \n" +
                  f"TBM -> Beta : {beta} deg \n" +
                  f"       Weak shock : {bool_weak} \n")

        b_rad = np.deg2rad(beta)
        t_rad = np.deg2rad(theta)

        NState = State(
                M = self.M * np.sin(b_rad),
                rho = self.rho,
                T = self.T,
                P = self.P,
                gamma=self.gamma,
                ID=-1
                ).NormalShock(bool_print = bool_print)

        OutputState = State(
                M = NState.M / np.sin(b_rad - t_rad),
                rho = NState.rho,
                T = NState.T,
                P = NState.P,
                gamma=self.gamma,
                ID=self.ID+1)

        return OutputState

    def PrandtlMeyer(self,theta,bool_print=False):
        theta = np.deg2rad(theta)
        nuM2 = theta + th.PrandtlMeyerFunction(self.M)
        M2 = th.InvPMF(nuM2)
        M1 = self.M
        Pr = th.P_ratio_PM(M1,M2,self.gamma)
        Tr = th.T_ratio_PM(M1,M2,self.gamma)
        Rhor = th.RHO_ratio_PM(M1,M2,self.gamma)
        outputState = State(
                M = M2,
                rho = Rhor * self.rho,
                P = Pr * self.P,
                T = Tr * self.T,
                gamma=self.gamma
                )
        if bool_print:
            str = "-- Prandtl-Meyer ratios : -- \n"
            str+= "---------------------------- \n"
            str+=f" - U2/U1 = {outputState.U / self.U} \n"
            str+=f" - M2/M1 = {M2/M1} \n"
            str+= "---------------------------- \n"
            str+=f" - Rho2/Rho1 = {Rhor} \n"
            str+=f" - P2/P1 = {Pr} \n"
            str+=f" - T2/T1 = {Tr} \n"
            str+= "---------------------------- \n"
            str+=f" - P02/P01 = {outputState.TotalState.P / self.TotalState.P} \n"
            str+=f" - T02/T01 = {outputState.TotalState.T / self.TotalState.T} \n"
            str+=f" - Rho02/Rho01 = {outputState.TotalState.rho / self.TotalState.rho} \n"
            str+= "---------------------------- \n"
            str+=f" - Mu1 = {np.rad2deg(np.arcsin(1/M1))} \n"
            str+=f" - Mu2 = {np.rad2deg(np.arcsin(1/M2))} \n"
            str+= "---------------------------- \n"
            print(str)

        return outputState
