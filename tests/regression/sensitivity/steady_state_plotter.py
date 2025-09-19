import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('steady_state_data.csv')

reset_indices = df.index[df['time'].diff() < 0].tolist()

segment_bounds = [0] + reset_indices + [len(df)]

# Variables to plot and colours
vars_to_plot = {
    'u': 'tab:blue',
    'v': 'tab:orange',
    'w': 'tab:green',
    'p': 'tab:red',
    't': 'tab:purple'
}

plt.figure(figsize=(10, 6))

# Plot each segment
for i in range(len(segment_bounds) - 1):
    segment = df.iloc[segment_bounds[i]:segment_bounds[i + 1]]
    for var, color in vars_to_plot.items():
        if i == 0:  # only label the first segment
            plt.semilogy(segment['time'], segment[var], color=color, label=var)
        else:
            plt.semilogy(segment['time'], segment[var], color=color)

plt.xlabel('Time')
plt.ylabel('Values (log scale)')
plt.title('u, v, w, p, t vs Time (Semilog Y)')
plt.legend(ncol=3)
plt.grid(True, which="both", linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()

