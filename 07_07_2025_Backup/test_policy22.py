import numpy as np
import tensorflow as tf
from ns3gym import ns3env
from fast_ddqn_multiagent import MultiAgentWrapper

# === Test Configuration ===
TEST_EPISODES = 100
MAX_STEPS = 1000
N_AGENTS = 3
STATE_DIM_PER_AGENT = 5
ACTION_DIM = 4
SEEDS = [1249, 2454, 3366]  # different test seeds
SIM_TIME = 10
PACKET_INTERVAL = 0.05
PACKET_SIZE_BYTES = 1024

def run_test_episode(wrapper, seed):
    env = ns3env.Ns3Env(port=5555, stepTime=0.01, startSim=True, simSeed=seed)
    obs = env.reset()
    agent_states = wrapper.split_obs(obs)

    total_energy = 0.0
    avg_sinr = 0.0
    step_count = 0
    active_ue_count = []
    sbs_states_per_step = []
    total_reward = 0.0

    while step_count < MAX_STEPS:
        actions = wrapper.act(agent_states, episode_num=TEST_EPISODES + 1)  # ε=0 now
        next_obs, reward, done, info = env.step(np.array(actions, dtype=np.uint32))

        info_str = info if isinstance(info, str) else info[0]
        info_parts = dict(item.split("=") for item in info_str.split(";"))
        energy = float(info_parts.get("total_energy", 0.0))
        sinr = float(info_parts.get("global_sinr", 0.0))
        active_ue = int(info_parts.get("active_ue", 0))
        total_reward += reward

        total_energy = energy
        avg_sinr += sinr
        active_ue_count.append(active_ue)
        agent_states = wrapper.split_obs(next_obs)
        sbs_states = [int(state[1]) for state in agent_states]
        sbs_states_per_step.append(sbs_states)
        step_count += 1

        if done:
            break

    env.close()
    return total_energy, avg_sinr / step_count if step_count > 0 else 0.0, active_ue_count, sbs_states_per_step, total_reward

def run_baseline_episode(seed):
    env = ns3env.Ns3Env(port=5556, stepTime=0.01, startSim=True, simSeed=seed)
    obs = env.reset()

    total_energy = 0.0
    avg_sinr = 0.0
    total_reward = 0.0
    step_count = 0

    while step_count < MAX_STEPS:
        actions = np.zeros(N_AGENTS, dtype=np.uint32)  # Always ACTIVE
        next_obs, reward, done, info = env.step(actions)

        info_str = info if isinstance(info, str) else info[0]
        info_parts = dict(item.split("=") for item in info_str.split(";"))
        energy = float(info_parts.get("total_energy", 0.0))
        sinr = float(info_parts.get("global_sinr", 0.0))

        total_energy = energy
        avg_sinr += sinr
        total_reward += reward
        step_count += 1

        if done:
            break

    env.close()
    return total_energy, avg_sinr / step_count if step_count > 0 else 0.0, total_reward


if __name__ == "__main__":
    wrapper = MultiAgentWrapper(N_AGENTS, STATE_DIM_PER_AGENT, ACTION_DIM)
    wrapper.load_all("trained_ddqn_agent")

    # Disable exploration
    for agent in wrapper.agents:
        agent.epsilon = 0.0

    test_energies = []
    test_sinrs = []
    test_efficiencies = []
    test_active_ue = []
    test_sbs_states = []
    test_rewards = []

    packets_per_ue = SIM_TIME / PACKET_INTERVAL
    total_packets = packets_per_ue * 30  # assuming 30 UEs
    total_bits = total_packets * PACKET_SIZE_BYTES * 8

    for seed in SEEDS:
        energy, sinr, active_ue, sbs_states, total_reward = run_test_episode(wrapper, seed)
        test_energies.append(energy)
        test_sinrs.append(sinr)
        test_efficiencies.append(total_bits / energy if energy > 0 else 0.0)
        test_active_ue.append(active_ue)
        test_sbs_states.append(sbs_states)
        test_rewards.append(total_reward)

    # Save results
    np.save("test_energy_per_episode.npy", np.array(test_energies))
    np.save("test_sinr_per_episode.npy", np.array(test_sinrs))
    np.save("test_efficiency_per_episode.npy", np.array(test_efficiencies))
    np.save("test_active_ue_per_step.npy", np.array(test_active_ue, dtype=object))
    np.save("test_total_reward.npy", np.array(test_rewards))

    # CSV log
    import csv
    with open("test_active_ue_and_sbs_state_per_step.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Seed", "Step", "Active_UEs", "SBS_0_State", "SBS_1_State", "SBS_2_State"])
        for i, seed in enumerate(SEEDS):
            for step in range(len(test_active_ue[i])):
                row = [seed, step, test_active_ue[i][step]] + test_sbs_states[i][step]
                writer.writerow(row)


    # === Run Baseline Test (always ACTIVE) ===
    baseline_energies = []
    baseline_sinrs = []
    baseline_rewards = []

    print("\n=== Running Baseline (Always ACTIVE) ===")
    for seed in SEEDS:
        energy, sinr, reward = run_baseline_episode(seed)
        baseline_energies.append(energy)
        baseline_sinrs.append(sinr)
        baseline_rewards.append(reward)
        print(f"Seed {seed} | Baseline Energy: {energy:.2f} J | SINR: {sinr:.2f} dB | Reward: {reward:.2f}")

    # Save baseline results
    np.save("baseline_energy.npy", np.array(baseline_energies))
    np.save("baseline_sinr.npy", np.array(baseline_sinrs))
    np.save("baseline_reward.npy", np.array(baseline_rewards))


    # Console Summary
    print("=== Test Results ===")
    for i, seed in enumerate(SEEDS):
        print(f"Seed {seed} | Energy: {test_energies[i]:.2f} J | SINR: {test_sinrs[i]:.2f} dB | EE: {test_efficiencies[i]:.2f} bits/J")
        print(f"Seed {seed} | Total Reward: {test_rewards[i]:.2f}")

    # Plot Reward
    import matplotlib.pyplot as plt
    try:
        print("📈 Plotting total reward per seed...")
        rewards = np.load("test_total_reward.npy")
        plt.figure(figsize=(8, 5))
        plt.bar(range(len(rewards)), rewards, color='purple')
        plt.title("Total Reward per Seed (Test Episodes)")
        plt.xlabel("Test Case Index")
        plt.ylabel("Cumulative Reward")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("test_reward_summary.png")
        print("✅ Reward plot saved as 'test_reward_summary.png'")
    except Exception as e:
        print("❌ Failed to plot total reward:", e)

    # === Comparison Plot: Test vs Baseline ===
    try:
        print("📊 Plotting Test vs Baseline Comparison...")
        import matplotlib.pyplot as plt

        test_energy = np.load("test_energy_per_episode.npy")
        test_sinr = np.load("test_sinr_per_episode.npy")
        test_efficiency = np.load("test_efficiency_per_episode.npy")
        test_reward = np.load("test_total_reward.npy")

        baseline_energy = np.load("baseline_energy.npy")
        baseline_sinr = np.load("baseline_sinr.npy")
        baseline_reward = np.load("baseline_reward.npy")

        # Extend arrays to include baseline mean
        energies = np.append(test_energy, np.mean(baseline_energy))
        sinrs = np.append(test_sinr, np.mean(baseline_sinr))
        efficiencies = np.append(test_efficiency, np.mean(test_efficiency))  # placeholder (baseline not available)
        rewards = np.append(test_reward, np.mean(baseline_reward))

        labels = [f"Seed {i+1}" for i in range(len(test_energy))] + ["Baseline"]

        fig, axs = plt.subplots(4, 1, figsize=(10, 16), constrained_layout=True)

        axs[0].bar(labels, energies, color='skyblue')
        axs[0].set_title("Total Energy Consumption (Test vs Baseline)")
        axs[0].set_ylabel("Energy (J)")

        axs[1].bar(labels, sinrs, color='orange')
        axs[1].set_title("Average SINR (Test vs Baseline)")
        axs[1].set_ylabel("SINR (dB)")

        axs[2].bar(labels, efficiencies, color='green')
        axs[2].set_title("Energy Efficiency (Test vs Baseline)")
        axs[2].set_ylabel("Bits/Joule")

        axs[3].bar(labels, rewards, color='purple')
        axs[3].set_title("Total Reward (Test vs Baseline)")
        axs[3].set_ylabel("Cumulative Reward")

        for ax in axs:
            ax.set_xlabel("Scenario")
            ax.grid(True)

        plt.savefig("test_vs_baseline_comparison.png")
        print("✅ Comparison plot saved as 'test_vs_baseline_comparison.png'")

    except Exception as e:
        print("❌ Failed to generate comparison plot:", e)

    # === Plot Active UE Count Per Step ===
    try:
        print("📈 Plotting Active UE Count per Hour...")
        test_active_ue = np.load("test_active_ue_per_step.npy", allow_pickle=True)

        plt.figure(figsize=(10, 6))
        for i, ue_counts in enumerate(test_active_ue):
            steps = len(ue_counts)
            hours = np.linspace(0, 24, steps)  # Map 0–1000 steps to 0–24 hours
            plt.plot(hours, ue_counts, label=f'Seed {i+1}', alpha=0.6)

        plt.title("Active UE Count over Simulated 24 Hours")
        plt.xlabel("Hour of Day")
        plt.ylabel("Number of Active UEs")
        plt.grid(True)
        plt.legend(loc="upper right", fontsize="small")
        plt.tight_layout()
        plt.savefig("active_ue_count_per_hour.png")
        print("✅ Saved: 'active_ue_count_per_hour.png'")
    except Exception as e:
        print("❌ Failed to plot Active UE Count:", e)

    # Run post-analysis script
    import subprocess
    try:
        print("📊 Running post-analysis script...")
        subprocess.run(["python3", "plot_sbs_analysis_multi_seed.py"], check=True)
        print("✅ Analysis script completed successfully.")
    except subprocess.CalledProcessError as e:
        print("❌ Failed to run plotting script:", e)
