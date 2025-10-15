import CompSolver.lib.state as cs

def solveNormalShock():
    """
    Function to solve normal shock relations given initial state and either Mach number or velocity.
    """

    print("\n\n####### NORMAL SHOCK RESOLUTION #######")
    print("####### First state definition :")
    definition = input("## Use the 15C/0m/Patm/g1.4 conditions ? (Y/N) ")
    p = None
    t = None
    r = None
    g = None
    if not definition.lower() == "y":
        try: 
            p = float(input("## Pressure (Pa)   : "))
        except:
            pass
        try:
            t = float(input("## Temperature (K) : "))
        except:
            pass
        try:
            r = float(input("## Density (kg/m3) : "))
        except:
            pass
        try:
            g = float(input("## Gamma           : "))
        except:
            pass
    
    vel_input = input("## Velocity type : M or U ? (M/U) ")
    if vel_input.lower() == "m":
        m = float(input("## Velocity (Mach) : "))
        u = None
    else:
        m = None
        u = float(input("## Velocity (m/s)  : "))
    print("####### RESULTS : \n")
    state0 = cs.State( M=m, U=u, P=p, T=t, rho=r, gamma=g )
    state1 = state0.NormalShock(bool_print=True)
    print(state0)
    print(state1)
