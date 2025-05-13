import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Import updated deployment angle arrays
from deployment_simulation_v2 import theta_root, theta_1, theta_2, theta_3, theta_4

# Adjust angles for animation visualization
angles_root = theta_root + np.pi / 2
angles_1 = (theta_1 / 2) + np.pi / 2
angles_2 = (theta_2 / 2) + np.pi / 2
angles_3 = (theta_3 / 2) + np.pi / 2
angles_4 = (theta_4 / 2) + np.pi / 2

# Panel lengths
L_root = 0.385  # m
L_panel = 0.720  # m

# Frame settings
fps = 10
speedup_factor = 20
num_frames = len(angles_root)
dt = 0.01

# Function to compute joint positions
def compute_positions(theta_r, theta_1, theta_2, theta_3, theta_4):
    x0, y0 = 0, 0
    x1 = x0 + L_root * np.sin(theta_r)
    y1 = y0 - L_root * np.cos(theta_r)
    
    x2 = x1 + L_panel * np.sin(theta_1)
    y2 = y1 + L_panel * np.cos(theta_1)
    
    x3 = x2 + L_panel * np.sin(theta_2)
    y3 = y2 - L_panel * np.cos(theta_2)
    
    x4 = x3 + L_panel * np.sin(theta_3)
    y4 = y3 + L_panel * np.cos(theta_3)
    
    x5 = x4 + L_panel * np.sin(theta_4)
    y5 = y4 - L_panel * np.cos(theta_4)

    return [(x0, y0), (x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5)]

# Plot stowed and deployed configurations
fig1, ax = plt.subplots(1, 2, figsize=(12, 6))

# --- Stowed Configuration ---
ax[0].set_title("Stowed Configuration")
stowed = compute_positions(np.pi, np.pi, np.pi, np.pi, np.pi)
st_x, st_y = zip(*stowed)

# Draw panels
ax[0].plot([st_x[0], st_x[1]], [st_y[0], st_y[1]], color='red', lw=2)
ax[0].plot([st_x[1], st_x[2]], [st_y[1], st_y[2]], color='blue', lw=2)
ax[0].plot([st_x[2], st_x[3]], [st_y[2], st_y[3]], color='green', lw=2)
ax[0].plot([st_x[3], st_x[4]], [st_y[3], st_y[4]], color='purple', lw=2)
ax[0].plot([st_x[4], st_x[5]], [st_y[4], st_y[5]], color='orange', lw=2)

# Add black dots at hinge joints (not the last tip)
ax[0].scatter(st_x[:5], st_y[:5], color='black', zorder=5)
ax[0].set_xlim(-1, 1)
ax[0].set_ylim(-1.5, 1.5)
ax[0].set_aspect('equal')

# --- Deployed Configuration ---
ax[1].set_title("Deployed Configuration")
deployed = compute_positions(np.pi/2, np.pi/2, np.pi/2, np.pi/2, np.pi/2)
dep_x, dep_y = zip(*deployed)

# Draw panels
ax[1].plot([dep_x[0], dep_x[1]], [dep_y[0], dep_y[1]], color='red', lw=2)
ax[1].plot([dep_x[1], dep_x[2]], [dep_y[1], dep_y[2]], color='blue', lw=2)
ax[1].plot([dep_x[2], dep_x[3]], [dep_y[2], dep_y[3]], color='green', lw=2)
ax[1].plot([dep_x[3], dep_x[4]], [dep_y[3], dep_y[4]], color='purple', lw=2)
ax[1].plot([dep_x[4], dep_x[5]], [dep_y[4], dep_y[5]], color='orange', lw=2)

# Add black dots at hinge joints
ax[1].scatter(dep_x[:5], dep_y[:5], color='black', zorder=5)
ax[1].set_xlim(-1, 5)
ax[1].set_ylim(-1.5, 1.5)
ax[1].set_aspect('equal')

plt.legend(["Root", "Panel 1", "Panel 2", "Panel 3", "Panel 4"], loc='lower right')
plt.tight_layout()
plt.show()



# Animation setup
fig2, ax2 = plt.subplots(figsize=(8, 6))
ax2.set_title("Radiator Deployment Animation")
ax2.set_xlim(-1, 5)
ax2.set_ylim(-1.5, 1.5)
ax2.set_aspect('equal')

# Initialize plot lines for panels
line_root, = ax2.plot([], [], 'r-', lw=2)     # Root panel - red
line_1, = ax2.plot([], [], 'b-', lw=2)        # Panel 1 - blue
line_2, = ax2.plot([], [], 'g-', lw=2)        # Panel 2 - green
line_3, = ax2.plot([], [], 'm-', lw=2)        # Panel 3 - purple (magenta)
line_4, = ax2.plot([], [], '-', color='orange', lw=2)  # Panel 4 - orange

# Add black hinge joint markers (will be updated every frame)
hinge_dots = ax2.scatter([], [], color='black', zorder=5)

# Add time text in bottom right
time_text = ax2.text(0.95, 0.05, '', transform=ax2.transAxes,
                     ha='right', va='bottom', fontsize=12,
                     bbox=dict(facecolor='white', alpha=0.7, edgecolor='black'))

# Real deployment time
deployment_time = dt * len(theta_root)

def update(frame):
    positions = compute_positions(
        angles_root[frame], angles_1[frame], angles_2[frame],
        angles_3[frame], angles_4[frame]
    )
    x, y = zip(*positions)
    
    # Update panel lines
    line_root.set_data([x[0], x[1]], [y[0], y[1]])
    line_1.set_data([x[1], x[2]], [y[1], y[2]])
    line_2.set_data([x[2], x[3]], [y[2], y[3]])
    line_3.set_data([x[3], x[4]], [y[3], y[4]])
    line_4.set_data([x[4], x[5]], [y[4], y[5]])

    # Update hinge dots
    hinge_dots.set_offsets(np.column_stack((x[:5], y[:5])))

    # Update time text (real-world time)
    sim_time = frame * dt
    time_text.set_text(f"t = {sim_time:.2f} s")

    return line_root, line_1, line_2, line_3, line_4, hinge_dots, time_text

# Set how many times faster you want the playback to be
speedup_factor = 5

# Drop frames accordingly
frame_indices = np.arange(0, num_frames, speedup_factor)

ani = animation.FuncAnimation(
    fig2, update, frames=frame_indices,
    interval=1000 / 30, blit=True  # 30 fps playback
)

ani.save("z_fold_deployment_ideal.gif", writer="pillow", fps=30)


