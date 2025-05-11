import pandas as pd

data = pd.read_csv('../optimization_data.csv')
iter = pd.DataFrame(data['iter']).to_numpy()
obj = pd.DataFrame(data[' Scalar Mixing']).to_numpy()
KKT = pd.DataFrame(data[' KKTmax']).to_numpy()
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import matplotlib as mpl

fig, ax = plt.subplots(2,1, figsize=(6, 7), dpi = 200)
ax[0].plot(iter, obj, '-k')
ax[0].set_ylabel("Objective function")
ax[0].set_xlim(iter[0], iter[-1])
ax[1].plot(iter[1:-1], KKT[1:-1], '-k')
ax[1].set_xlabel("iteration")
ax[1].set_ylabel("KKT")
ax[1].set_xlim(iter[0], iter[-1])

fig.savefig('convergence_history.png')
