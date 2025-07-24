import numpy as np
import matplotlib.pyplot as plt

# === Skip first 15 episodes ===
skip = 15

# === Load sliced data ===
avg_reward = np.load("avg_reward_per_episode.npy")[skip:]
sinr = np.load("sinr_per_episode.npy")[skip:]
energy_efficiency = np.load("energy_efficiency_per_episode.npy")[skip:]
baseline_efficiency = np.load("energy_efficiency_baseline.npy")[skip:]
rl_energy = np.load("rl_energy_per_episode.npy")[skip:]
baseline_energy = np.load("baseline_energy_per_episode.npy")[skip:]
rl_power = np.load("rl_energy_per_episode.npy")[skip:]  # same file as energy
baseline_power = np.load("baseline_power_per_episode.npy")[skip:]
baseline_sinr = np.load("baseline_sinr_per_episode.npy")[skip:]

# === Adjust episode numbers ===
episodes = np.arange(skip + 1, skip + 1 + len(avg_reward))

# === Plot: Avg Reward ===
plt.figure(figsize=(10, 6))
plt.plot(episodes, avg_reward, label="Avg Total Reward")
plt.xlabel("Episode")
plt.ylabel("Reward")
plt.title("Avg Reward Per Episode (Post-Pretraining)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("posttrain_avg_reward.png")

# === Plot: SINR Comparison ===
plt.figure(figsize=(10, 6))
plt.plot(episodes, sinr, label="With RL (DDQN)")
plt.plot(episodes, baseline_sinr, label="Baseline", linestyle="--", color="orange")
plt.xlabel("Episode")
plt.ylabel("Global SINR (dB)")
plt.title("SINR Comparison: RL vs Baseline (Post-Pretraining)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("posttrain_sinr_comparison.png")

# === Plot: Energy Efficiency Comparison ===
plt.figure(figsize=(10, 6))
plt.plot(episodes, energy_efficiency, label="With RL (DDQN)")
plt.plot(episodes, baseline_efficiency, label="Baseline")
plt.xlabel("Episode")
plt.ylabel("Energy Efficiency (bits/Joule)")
plt.title("Energy Efficiency Comparison (Post-Pretraining)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("posttrain_energy_efficiency.png")

# === Plot: Energy Consumption Comparison ===
plt.figure(figsize=(10, 6))
plt.plot(episodes, rl_energy, label="With RL (DDQN)")
plt.plot(episodes, baseline_energy, label="Baseline")
plt.xlabel("Episode")
plt.ylabel("Total Energy (Joules)")
plt.title("Energy Consumption Comparison (Post-Pretraining)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("posttrain_energy_consumption.png")

# === Plot: Power Comparison ===
plt.figure(figsize=(10, 6))
plt.plot(episodes, rl_power, label="With RL (DDQN)")
plt.plot(episodes, baseline_power, label="Baseline")
plt.xlabel("Episode")
plt.ylabel("Total Power (Watts)")
plt.title("Power Consumption Comparison (Post-Pretraining)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("posttrain_power_comparison.png")

print("✅ Plots saved excluding first 15 episodes (treated as pre-training).")
