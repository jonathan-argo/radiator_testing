import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Import your real deployment angles
from deployment_simulation_v2 import (
    theta_root_real, theta_1_real, theta_2_real, theta_3_real, theta_4_real
)

# Convert to visualization angles
angles_root = theta_root_real + np.pi / 2
angles_1 = (theta_1_real / 2) + np.pi / 2
angles_2 = (theta_2_real / 2) + np.pi / 2
angles_3 = (theta_3_real / 2) + np.pi / 2
angles_4 = (theta_4_real / 2) + np.pi / 2

# Geometry
L_root = 0.385
L_panel = 0.720
dt = 0.01
fps = 10
speedup_factor = 20

# Determine animation length
def get_deploy_frame(theta):
    idx = np.where(theta <= 0)[0]
    return idx[0] if len(idx) > 0 else len(theta)

num_frames = max(map(get_deploy_frame, [
    theta_root_real, theta_1_real, theta_2_real, theta_3_real, theta_4_real
])) + 1

deployment_time = dt * num_frames

def compute_positions(th_r, th1, th2, th3, th4):
    x0, y0 = 0, 0
    x1, y1 = x0 + L_root * np.sin(th_r), y0 - L_root * np.cos(th_r)
    x2, y2 = x1 + L_panel * np.sin(th1), y1 + L_panel * np.cos(th1)
    x3, y3 = x2 + L_panel * np.sin(th2), y2 - L_panel * np.cos(th2)
    x4, y4 = x3 + L_panel * np.sin(th3), y3 + L_panel * np.cos(th3)
    x5, y5 = x4 + L_panel * np.sin(th4), y4 - L_panel * np.cos(th4)
    return [(x0, y0), (x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5)]

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_title("Radiator Deployment with Tracked Midpoints")
ax.set_xlim(-1, 5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal')

# Panel lines
line_root, = ax.plot([], [], 'r-', lw=2)
line_1, = ax.plot([], [], 'b-', lw=2)
line_2, = ax.plot([], [], 'g-', lw=2)
line_3, = ax.plot([], [], 'm-', lw=2)
line_4, = ax.plot([], [], '-', color='orange', lw=2)

# Hinge joints (large black)
hinge_dots = ax.scatter([], [], color='black', zorder=5, s=20)

# Midpoint dots (small, dark)
mid_dots = {
    "root": ax.plot([], [], 'o', color='#800000', markersize=4)[0],   # dark red
    "p1": ax.plot([], [], 'o', color='#000080', markersize=4)[0],     # dark blue
    "p2": ax.plot([], [], 'o', color='#006400', markersize=4)[0],     # dark green
    "p3": ax.plot([], [], 'o', color='#4B0082', markersize=4)[0],     # dark purple
    "p4": ax.plot([], [], 'o', color='#8B4513', markersize=4)[0],     # dark orange
}

# Trace lines (omit root)
trace_lines = {
    "p1": ax.plot([], [], '-', color='#ADD8E6', lw=1)[0],
    "p2": ax.plot([], [], '-', color='#90EE90', lw=1)[0],
    "p3": ax.plot([], [], '-', color='#D8BFD8', lw=1)[0],
    "p4": ax.plot([], [], '-', color='#FFE4B5', lw=1)[0],
}
trace_coords = {key: ([], []) for key in trace_lines}

# Time text
time_text = ax.text(0.95, 0.05, '', transform=ax.transAxes,
                    ha='right', va='bottom', fontsize=12,
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='black'))

def update(frame):
    pts = compute_positions(
        angles_root[frame], angles_1[frame], angles_2[frame],
        angles_3[frame], angles_4[frame]
    )
    x, y = zip(*pts)

    # Panel lines
    line_root.set_data([x[0], x[1]], [y[0], y[1]])
    line_1.set_data([x[1], x[2]], [y[1], y[2]])
    line_2.set_data([x[2], x[3]], [y[2], y[3]])
    line_3.set_data([x[3], x[4]], [y[3], y[4]])
    line_4.set_data([x[4], x[5]], [y[4], y[5]])

    # Hinges
    hinge_dots.set_offsets(np.column_stack((x[:5], y[:5])))

    # Midpoints + trace buildup
    mids = {
        "root": ((x[0] + x[1]) / 2, (y[0] + y[1]) / 2),
        "p1": ((x[1] + x[2]) / 2, (y[1] + y[2]) / 2),
        "p2": ((x[2] + x[3]) / 2, (y[2] + y[3]) / 2),
        "p3": ((x[3] + x[4]) / 2, (y[3] + y[4]) / 2),
        "p4": ((x[4] + x[5]) / 2, (y[4] + y[5]) / 2),
    }

    for key, (xm, ym) in mids.items():
        mid_dots[key].set_data([xm], [ym])

        if key in trace_lines:
            if frame == 0:
                trace_coords[key][0].clear()
                trace_coords[key][1].clear()
                trace_lines[key].set_data([], [])
            else:
                trace_coords[key][0].append(xm)
                trace_coords[key][1].append(ym)
                trace_lines[key].set_data(trace_coords[key][0], trace_coords[key][1])

    time_text.set_text(f"t = {frame * dt:.2f} s")

    return (
        line_root, line_1, line_2, line_3, line_4,
        hinge_dots, time_text,
        *mid_dots.values(), *trace_lines.values()
    )


# Set how many times faster you want the playback to be
speedup_factor = 5

# Drop frames accordingly
frame_indices = np.arange(0, num_frames, speedup_factor)

ani = animation.FuncAnimation(
    fig, update, frames=frame_indices,
    interval=1000 / 30, blit=True  # 30 fps playback
)

ani.save("z_fold_deployment_tracking.gif", writer="pillow", fps=30)