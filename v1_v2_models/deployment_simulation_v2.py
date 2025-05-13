import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------
# Toggle plot visibility (set to True to show)
# ---------------------------------------
plot_flags = {
    "plot_1_constant_torque_vs_time": False,
    "plot_2_max_torque_vs_time": False,
    "plot_3_overlay_constant_vs_spring": True,
    "plot_4_torque_vs_time_fixedT": True,
    "plot_5_deployment_angle_const_torque": False,
    "plot_6_deployment_angle_spring": False,
    "plot_7_overlay_angle_spring_vs_const": False,
    "plot_8_user_defined_spring": True,
}

# Updated panel properties
m_root = 0.2252  # kg
m_panel = 0.4185  # kg (for each full panel)
w_root = 0.385  # Root width (m)
w_panel = 0.720  # Panel width (m)

# Angular inertia of each panel (as a rod rotating about one end)
I_root = (1/3) * m_root * w_root**2
I_panel = (1/3) * m_panel * w_panel**2

# Initial angles
theta_initial_root = np.deg2rad(90)   # root starts at 90 degrees
theta_initial_other = np.deg2rad(180) # all others start folded

# Compute effective moments of inertia for each hinge (panels it must deploy)
I_root_eff = I_root + 4 * I_panel
I_1_eff = 4 * I_panel
I_2_eff = 3 * I_panel
I_3_eff = 2 * I_panel
I_4_eff = 1 * I_panel

# Target deployment time (seconds)
deployment_time = 20  # <- change this to whatever you want

# Use angular displacement values
theta_0_root = theta_initial_root        # 90° → 0
theta_0_other = theta_initial_other      # 180° → 0

# Time stepping
dt = 0.01
t_max = 80
time = np.arange(0, t_max, dt)

# Initialize angles and angular velocities
theta_root = np.full_like(time, theta_initial_root)
theta_1 = np.full_like(time, theta_initial_other)
theta_2 = np.full_like(time, theta_initial_other)
theta_3 = np.full_like(time, theta_initial_other)
theta_4 = np.full_like(time, theta_initial_other)


# ---------------------------------------
# Plot 1: Constant Torque vs Deployment Time
# ---------------------------------------
deployment_times = np.linspace(5, 60, 300)  # Range of deployment durations

# Compute required torque curves
T_root_required_curve = I_root_eff * (2 * theta_0_root) / deployment_times**2
T_1_required_curve = I_1_eff * (2 * theta_0_other) / deployment_times**2
T_2_required_curve = I_2_eff * (2 * theta_0_other) / deployment_times**2
T_3_required_curve = I_3_eff * (2 * theta_0_other) / deployment_times**2
T_4_required_curve = I_4_eff * (2 * theta_0_other) / deployment_times**2

# Plotting
if plot_flags["plot_1_constant_torque_vs_time"]:
    plt.figure(figsize=(10,6))
    plt.plot(deployment_times, T_root_required_curve, label="Root Hinge", color="red")
    plt.plot(deployment_times, T_1_required_curve, label="Hinge 1-2", color="blue")
    plt.plot(deployment_times, T_2_required_curve, label="Hinge 2-3", color="green")
    plt.plot(deployment_times, T_3_required_curve, label="Hinge 3-4", color="purple")
    plt.plot(deployment_times, T_4_required_curve, label="Hinge 4-5", color="orange")

    plt.xlabel("Deployment Time (s)")
    plt.ylabel("Required Torque (N·m)")
    plt.title("Required Hinge Torque vs Total Deployment Time (V1 Model)", fontweight='bold')
    plt.grid()
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.show()

# ---------------------------------------
# Plot 2: Max Torque (Torsion Spring) vs Deployment Time
# ---------------------------------------
# Derived stiffness and max torque
def get_k_and_tau_max(I, theta_0, T):
    k = (np.pi**2 * I) / (4 * T**2)
    tau_max = k * theta_0
    return k, tau_max

k_root, tau_max_root = get_k_and_tau_max(I_root_eff, theta_0_root, deployment_times)
k_1, tau_max_1 = get_k_and_tau_max(I_1_eff, theta_0_other, deployment_times)
k_2, tau_max_2 = get_k_and_tau_max(I_2_eff, theta_0_other, deployment_times)
k_3, tau_max_3 = get_k_and_tau_max(I_3_eff, theta_0_other, deployment_times)
k_4, tau_max_4 = get_k_and_tau_max(I_4_eff, theta_0_other, deployment_times)

# Plot: Max Torque vs Deployment Time (derived from ODE)
if plot_flags["plot_2_max_torque_vs_time"]:
    plt.figure(figsize=(10,6))
    plt.plot(deployment_times, tau_max_root, label="Root Hinge", color="red")
    plt.plot(deployment_times, tau_max_1, label="Hinge 1-2", color="blue")
    plt.plot(deployment_times, tau_max_2, label="Hinge 2-3", color="green")
    plt.plot(deployment_times, tau_max_3, label="Hinge 3-4", color="purple")
    plt.plot(deployment_times, tau_max_4, label="Hinge 4-5", color="orange")

    plt.xlabel("Deployment Time (s)")
    plt.ylabel("Max Torque (N·m)")
    plt.title("Max Torque for Torsion Springs vs Deployment Time", fontweight='bold')
    plt.grid()
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


# ---------------------------------------
# PLOT 3
# ---------------------------------------
if plot_flags["plot_3_overlay_constant_vs_spring"]:
    plt.figure(figsize=(10,6))

    # Torsion spring max torque curves (from ODE-based dynamics)
    plt.plot(deployment_times, tau_max_root, label="Root Hinge (Torsion Spring)", color="red")
    plt.plot(deployment_times, tau_max_1, label="Hinge 1-2 (Torsion Spring)", color="blue")
    plt.plot(deployment_times, tau_max_2, label="Hinge 2-3 (Torsion Spring)", color="green")
    plt.plot(deployment_times, tau_max_3, label="Hinge 3-4 (Torsion Spring)", color="purple")
    plt.plot(deployment_times, tau_max_4, label="Hinge 4-5 (Torsion Spring)", color="orange")

    # Constant torque reference curves (black dashed lines)
    plt.plot(deployment_times, T_root_required_curve, '--', color="black", alpha=0.7, label="Root Hinge (Constant Torque)")
    plt.plot(deployment_times, T_1_required_curve, '--', color="black", alpha=0.7)
    plt.plot(deployment_times, T_2_required_curve, '--', color="black", alpha=0.7)
    plt.plot(deployment_times, T_3_required_curve, '--', color="black", alpha=0.7)
    plt.plot(deployment_times, T_4_required_curve, '--', color="black", alpha=0.7)

    plt.xlabel("Deployment Time (s)")
    plt.ylabel("Torque (N·m)")
    plt.title("Torsion Spring vs Constant Torque Requirements", fontweight='bold')
    plt.grid()
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


# ---------------------------------------
# Plot 4: Torque vs Time During Deployment (at fixed deployment_time)
# ---------------------------------------
t_deploy = np.arange(0, deployment_time, dt)

# Torque as a function of time: τ(t) = τ_max * (1 - t/T)
T_root_time = T_root_required_curve[deployment_times == min(deployment_times, key=lambda x: abs(x - deployment_time))] * (1 - t_deploy / deployment_time)
T_1_time = T_1_required_curve[deployment_times == min(deployment_times, key=lambda x: abs(x - deployment_time))] * (1 - t_deploy / deployment_time)
T_2_time = T_2_required_curve[deployment_times == min(deployment_times, key=lambda x: abs(x - deployment_time))] * (1 - t_deploy / deployment_time)
T_3_time = T_3_required_curve[deployment_times == min(deployment_times, key=lambda x: abs(x - deployment_time))] * (1 - t_deploy / deployment_time)
T_4_time = T_4_required_curve[deployment_times == min(deployment_times, key=lambda x: abs(x - deployment_time))] * (1 - t_deploy / deployment_time)

# Plotting
if plot_flags["plot_4_torque_vs_time_fixedT"]:
    plt.figure(figsize=(10,6))
    plt.plot(t_deploy, T_root_time, label="Root Hinge", color="red")
    plt.plot(t_deploy, T_1_time, label="Hinge 1-2", color="blue")
    plt.plot(t_deploy, T_2_time, label="Hinge 2-3", color="green")
    plt.plot(t_deploy, T_3_time, label="Hinge 3-4", color="purple")
    plt.plot(t_deploy, T_4_time, label="Hinge 4-5", color="orange")

    plt.xlabel("Time (s)")
    plt.ylabel("Torque Output (N·m)")
    plt.title(f"Torsion Spring Torque During Deployment (T = {deployment_time}s)", fontweight='bold')
    plt.grid()
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.show()

# Plot 5: Deployment Angle vs Time (Constant Torque)
# θ(t) = θ₀ · (1 - (t / T))²
theta_const_root = theta_0_root * (1 - (t_deploy / deployment_time)**2)
theta_const_1 = theta_0_other * (1 - (t_deploy / deployment_time)**2)
theta_const_2 = theta_0_other * (1 - (t_deploy / deployment_time)**2)
theta_const_3 = theta_0_other * (1 - (t_deploy / deployment_time)**2)
theta_const_4 = theta_0_other * (1 - (t_deploy / deployment_time)**2)

if plot_flags["plot_5_deployment_angle_const_torque"]:
    plt.figure(figsize=(10,6))
    plt.plot(t_deploy, np.rad2deg(theta_const_root), label="Root Hinge", color="red")
    plt.plot(t_deploy, np.rad2deg(theta_const_1), label="Hinge 1-2", color="blue")
    plt.plot(t_deploy, np.rad2deg(theta_const_2), label="Hinge 2-3", color="green")
    plt.plot(t_deploy, np.rad2deg(theta_const_3), label="Hinge 3-4", color="purple")
    plt.plot(t_deploy, np.rad2deg(theta_const_4), label="Hinge 4-5", color="orange")

    plt.xlabel("Time (s)")
    plt.ylabel("Deployment Angle (deg)")
    plt.title(f"Plot 5: Deployment Angle vs Time (Constant Torque, T = {deployment_time}s)", fontweight='bold')
    plt.grid()
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.show()


# ---------------------------------------
# Plot 6: Deployment Angle vs Time (Torsion Spring – Hooke's Law, with torque & k in legend)
# ---------------------------------------

# Cosine-based angle profile for each hinge
theta_spring_root = theta_0_root * np.cos(np.pi * t_deploy / (2 * deployment_time))
theta_spring_1 = theta_0_other * np.cos(np.pi * t_deploy / (2 * deployment_time))
theta_spring_2 = theta_0_other * np.cos(np.pi * t_deploy / (2 * deployment_time))
theta_spring_3 = theta_0_other * np.cos(np.pi * t_deploy / (2 * deployment_time))
theta_spring_4 = theta_0_other * np.cos(np.pi * t_deploy / (2 * deployment_time))

# Plot with torque and stiffness in legend (rounded for readability)
if plot_flags["plot_6_deployment_angle_spring"]:
    plt.figure(figsize=(10,6))

    plt.plot(t_deploy, np.rad2deg(theta_spring_root), 
            label=f"Root Hinge\nk={k_root[0]:.4f} Nm/rad, τmax={tau_max_root[0]:.4f} Nm", color="red")

    plt.plot(t_deploy, np.rad2deg(theta_spring_1), 
            label=f"Hinge 1-2\nk={k_1[0]:.4f}, τmax={tau_max_1[0]:.4f}", color="blue")

    plt.plot(t_deploy, np.rad2deg(theta_spring_2), 
            label=f"Hinge 2-3\nk={k_2[0]:.4f}, τmax={tau_max_2[0]:.4f}", color="green")

    plt.plot(t_deploy, np.rad2deg(theta_spring_3), 
            label=f"Hinge 3-4\nk={k_3[0]:.4f}, τmax={tau_max_3[0]:.4f}", color="purple")

    plt.plot(t_deploy, np.rad2deg(theta_spring_4), 
            label=f"Hinge 4-5\nk={k_4[0]:.4f}, τmax={tau_max_4[0]:.4f}", color="orange")

    plt.xlabel("Time (s)")
    plt.ylabel("Deployment Angle (deg)")
    plt.title(f"Plot 6: Angle vs Time (Spring Driven, T = {deployment_time}s)", fontweight='bold')
    plt.grid()
    plt.legend(fontsize=9, loc="upper right")
    plt.tight_layout()
    plt.show()





# ---------------------------------------
# Plot 7: Overlay – Spring (Hooke's Law) vs Constant Torque
# ---------------------------------------

if plot_flags["plot_7_overlay_angle_spring_vs_const"]:
    plt.figure(figsize=(10,6))

    # Torsion Spring (Color – Hooke’s Law cosine model)
    plt.plot(t_deploy, np.rad2deg(theta_spring_root), label="Root Hinge (Spring)", color="red")
    plt.plot(t_deploy, np.rad2deg(theta_spring_1), label="Hinge 1-2 (Spring)", color="blue")
    plt.plot(t_deploy, np.rad2deg(theta_spring_2), label="Hinge 2-3 (Spring)", color="green")
    plt.plot(t_deploy, np.rad2deg(theta_spring_3), label="Hinge 3-4 (Spring)", color="purple")
    plt.plot(t_deploy, np.rad2deg(theta_spring_4), label="Hinge 4-5 (Spring)", color="orange")

    # Constant Torque (Black Dashed)
    plt.plot(t_deploy, np.rad2deg(theta_const_root), '--', color="black", alpha=0.7, label="Root Hinge (Const)")
    plt.plot(t_deploy, np.rad2deg(theta_const_1), '--', color="black", alpha=0.7)
    plt.plot(t_deploy, np.rad2deg(theta_const_2), '--', color="black", alpha=0.7)
    plt.plot(t_deploy, np.rad2deg(theta_const_3), '--', color="black", alpha=0.7)
    plt.plot(t_deploy, np.rad2deg(theta_const_4), '--', color="black", alpha=0.7)

    plt.xlabel("Time (s)")
    plt.ylabel("Deployment Angle (deg)")
    plt.title(f"Plot 7: Deployment Angle vs Time (Spring vs Constant Torque, T = {deployment_time}s)", fontweight='bold')
    plt.grid()
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()

# ---------------------------------------
# Make spring-driven angles the default for animation
# ---------------------------------------
theta_root = theta_spring_root
theta_1 = theta_spring_1
theta_2 = theta_spring_2
theta_3 = theta_spring_3
theta_4 = theta_spring_4


# ---------------------------------------
# Plot 8: Deployment Angle vs Time for Custom Springs
# ---------------------------------------

# Define your own torsional spring values (Nm/rad and Nm)
custom_springs = {
    "Root":   {"k": 0.0079,  "tau_max": 0.0373,  "I": I_root_eff, "theta_0": theta_0_root,  "color": "red"},
    "Hinge1": {"k": 0.0079,  "tau_max": 0.0373,  "I": I_1_eff,    "theta_0": theta_0_other, "color": "blue"},
    "Hinge2": {"k": 0.0079,  "tau_max": 0.0373,  "I": I_2_eff,    "theta_0": theta_0_other, "color": "green"},
    "Hinge3": {"k": 0.0079,  "tau_max": 0.0373,  "I": I_3_eff,    "theta_0": theta_0_other, "color": "purple"},
    "Hinge4": {"k": 0.0079,  "tau_max": 0.0373,  "I": I_4_eff,    "theta_0": theta_0_other, "color": "orange"},
}

if plot_flags["plot_8_user_defined_spring"]:
    plt.figure(figsize=(10,6))
    longest_deploy_time = 0

    for label, props in custom_springs.items():
        k = props["k"]
        tau_max = props["tau_max"]
        I = props["I"]
        theta_0 = props["theta_0"]
        color = props["color"]

        # Angular frequency and initial deployment curve
        omega = np.sqrt(k / I)
        theta_t = (tau_max / k) * np.cos(omega * t_deploy)

        # Find the first time angle reaches 0 or below
        deployed_index = np.where(theta_t <= 0)[0]

        if len(deployed_index) > 0:
            idx = deployed_index[0]
            theta_t[idx:] = 0  # Clamp to zero after that point
            deploy_time = t_deploy[idx]
            longest_deploy_time = max(longest_deploy_time, deploy_time)
            plt.axvline(deploy_time, color=color, linestyle='--', alpha=0.5)
            legend_label = (
                f"{label}\n"
                f"k = {k:.2f} Nm/rad, τₘₐₓ = {tau_max:.2f} Nm\n"
                f"t₀ = {deploy_time:.2f} s"
            )
        else:
            legend_label = (
                f"{label}\n"
                f"k = {k:.2f} Nm/rad, τₘₐₓ = {tau_max:.2f} Nm\n"
                f"t₀ > {t_deploy[-1]:.2f} s"
            )

        plt.plot(t_deploy, np.rad2deg(theta_t), label=legend_label, color=color)

    plt.xlabel("Time (s)")
    plt.ylabel("Deployment Angle (deg)")
    plt.title("Plot 8: Deployment with Manually Specified Springs", fontweight='bold')
    plt.grid()
    plt.xlim(0, longest_deploy_time + 2)  # Extend x-axis slightly past latest deploy
    plt.legend(fontsize=9)
    plt.tight_layout()
    plt.show()

# ---------------------------------------
# Export user-defined spring deployment angles (real spring response)
# ---------------------------------------

def compute_clamped_theta(k, tau_max, I, theta_0):
    omega = np.sqrt(k / I)
    theta = theta_0 * np.cos(omega * t_deploy)  # <- match animation logic
    deployed_idx = np.where(theta <= 0)[0]
    if len(deployed_idx) > 0:
        theta[deployed_idx[0]:] = 0
    return theta

theta_root_real = compute_clamped_theta(
    custom_springs["Root"]["k"], custom_springs["Root"]["tau_max"], I_root_eff, theta_0_root
)
theta_1_real = compute_clamped_theta(
    custom_springs["Hinge1"]["k"], custom_springs["Hinge1"]["tau_max"], I_1_eff, theta_0_other
)
theta_2_real = compute_clamped_theta(
    custom_springs["Hinge2"]["k"], custom_springs["Hinge2"]["tau_max"], I_2_eff, theta_0_other
)
theta_3_real = compute_clamped_theta(
    custom_springs["Hinge3"]["k"], custom_springs["Hinge3"]["tau_max"], I_3_eff, theta_0_other
)
theta_4_real = compute_clamped_theta(
    custom_springs["Hinge4"]["k"], custom_springs["Hinge4"]["tau_max"], I_4_eff, theta_0_other
)

