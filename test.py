import CompSolver as cs
import matplotlib.pyplot as plt
import numpy as np

M = np.linspace(5,35,50)
theta = np.linspace(1,89,89)

CD = {
    'exact':[],
    'newtonian':[],
    'sim':[]
}

CL = {
    'exact':[],
    'newtonian':[],
    'sim':[]
}

LD = {
    'exact':[],
    'newtonian':[],
    'sim':[]
}


for t in theta:
    s = cs.State(M=10)
    fp = cs.FlatPlate(theta=t,state0=s)
    fp.Solve()

    CD['exact'].append(fp.Cd['Exact'])
    CD['newtonian'].append(fp.Cd['Newtonian'])
    CD['sim'].append(fp.Cd['Similarity'])
    
    CL['exact'].append(fp.Cl['Exact'])
    CL['newtonian'].append(fp.Cl['Newtonian'])
    CL['sim'].append(fp.Cl['Similarity'])

    LD['exact'].append(fp.LD['Exact'])
    LD['newtonian'].append(fp.LD['Newtonian'])
    LD['sim'].append(fp.LD['Similarity'])


plt.plot(M,CL['exact'])
plt.show()
