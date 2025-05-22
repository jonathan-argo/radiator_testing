import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.animation import FFMpegWriter

# Load and downsample
df = pd.read_csv("state_sol.csv", usecols=range(6))


# Define which columns are angles and lengths
angle_cols = ['theta1', 'theta2', 'theta3', 'theta4', 'theta5']
df = df.iloc[::10].reset_index(drop=True)
length_cols = [0.385, 0.72, 0.72, 0.72, 0.72]
num_vectors = len(angle_cols)

# Set up the figure
fig, ax = plt.subplots(constrained_layout=True)
fig.set_size_inches(4, 6)

ax.set_title("Radiator Deployment Animation", pad=10)

lines = [ax.plot([], [], 'o-', markersize=4)[0] for _ in range(num_vectors)]
time_text = ax.text(0.96, 0.98, '', transform=ax.transAxes,
                    ha='right', va='top', fontsize=10, bbox=dict(facecolor='white', alpha=0.8))

# Compute bounds
ax.set_xlim(-0.9144, 0.9144)
ax.set_ylim(0, 3.6576)
ax.set_aspect('equal')

# Update function
def update(frame):
    x_prev, y_prev = 0, 0
    for i in range(num_vectors):
        theta = df.loc[frame, angle_cols[i]]
        length = length_cols[i]
        x_new = x_prev + length * np.cos(theta)
        y_new = y_prev + length * np.sin(theta)
        lines[i].set_data([x_prev, x_new], [y_prev, y_new])
        x_prev, y_prev = x_new, y_new

    time_text.set_text(f"t = {df.loc[frame, 'time']:.2f} s")
    return lines + [time_text]


# Animate
ani = FuncAnimation(fig, update, frames=len(df), interval=100, blit=True)
ani.save("radiator_deployment.gif", writer=PillowWriter(fps=10))

# writer = FFMpegWriter(fps=10)
# ani.save("radiator_deployment.mp4", writer=writer)