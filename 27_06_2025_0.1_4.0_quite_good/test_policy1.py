
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
SEEDS = [342,718,901]  # different test seeds
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

    while step_count < MAX_STEPS:
        actions = wrapper.act(agent_states, episode_num=TEST_EPISODES+1)  # ε=0 now
        next_obs, reward, done, info = env.step(np.array(actions, dtype=np.uint32))

        info_str = info if isinstance(info, str) else info[0]
        info_parts = dict(item.split("=") for item in info_str.split(";"))
        energy = float(info_parts.get("total_energy", 0.0))
        sinr = float(info_parts.get("global_sinr", 0.0))

        total_energy = energy
        avg_sinr += sinr
        agent_states = wrapper.split_obs(next_obs)
        step_count += 1

        if done:
            break

    env.close()
    return total_energy, avg_sinr / step_count if step_count > 0 else 0.0

if __name__ == "__main__":
    wrapper = MultiAgentWrapper(N_AGENTS, STATE_DIM_PER_AGENT, ACTION_DIM)
    wrapper.load_all("trained_ddqn_agent")

    # Disable exploration
    for agent in wrapper.agents:
        agent.epsilon = 0.0

    test_energies = []
    test_sinrs = []
    test_efficiencies = []

    packets_per_ue = SIM_TIME / PACKET_INTERVAL
    total_packets = packets_per_ue * 30  # assuming 30 UEs
    total_bits = total_packets * PACKET_SIZE_BYTES * 8

    for seed in SEEDS:
        energy, sinr = run_test_episode(wrapper, seed)
        test_energies.append(energy)
        test_sinrs.append(sinr)
        test_efficiencies.append(total_bits / energy if energy > 0 else 0.0)

    np.save("test_energy_per_episode.npy", np.array(test_energies))
    np.save("test_sinr_per_episode.npy", np.array(test_sinrs))
    np.save("test_efficiency_per_episode.npy", np.array(test_efficiencies))

    print("=== Test Results ===")
    for i, seed in enumerate(SEEDS):
        print(f"Seed {seed} | Energy: {test_energies[i]:.2f} J | SINR: {test_sinrs[i]:.2f} dB | EE: {test_efficiencies[i]:.2f} bits/J")
