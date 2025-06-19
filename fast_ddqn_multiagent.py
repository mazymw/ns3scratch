
import numpy as np
import random
from collections import deque
import tensorflow as tf
from tensorflow import keras
from ns3gym import ns3env
import matplotlib.pyplot as plt
import os
import subprocess
import csv

# Hyperparameters
BATCH_SIZE = 128
TARGET_UPDATE_FREQ = 100
MEM_SIZE = 100000
EPISODES = 1
MAX_STEPS = 1000
N_AGENTS = 3
STATE_DIM_PER_AGENT = 5
ACTION_DIM = 4

# Environment parameters
NUM_UES = 30
SIM_TIME = 10
PACKET_INTERVAL = 0.05
PACKET_SIZE_BYTES = 1024

class FastDDQNAgent:
    def __init__(self, state_size, action_size, agent_id):
        self.state_size = state_size
        self.action_size = action_size
        self.agent_id = agent_id

        self.memory = deque(maxlen=MEM_SIZE)
        self.gamma = 0.98
        self.epsilon = 1.0
        self.epsilon_min = 0.05
        self.epsilon_decay = 0.9999700431
        self.lr = 0.0003

        self.model = self._build_model()
        self.target_model = self._build_model()
        self.update_target_network()
        self.train_steps = 0
        self.loss_history = []

    def _build_model(self):
        model = keras.Sequential([
            keras.layers.Dense(64, activation='relu', input_shape=(self.state_size,)),
            keras.layers.Dense(64, activation='relu'),
            keras.layers.Dense(self.action_size, activation='linear')
        ])
        model.compile(optimizer=keras.optimizers.Adam(learning_rate=self.lr), loss='mse')
        return model

    def update_target_network(self):
        self.target_model.set_weights(self.model.get_weights())

    def remember(self, state, action, reward, next_state, done):
        self.memory.append((state, int(action), reward, next_state, done))
 
    def act(self, state, episode_num):
        MASKING_THRESHOLD = 30  # Apply masking only for first 30 episodes

        is_transitioning = state[-1]  # last dim in your state
        current_state = int(state[1])  # current SBS state

        # Apply masking only during early episodes
        if episode_num <= MASKING_THRESHOLD and is_transitioning:
            valid_actions = [current_state]
        else:
            valid_actions = list(range(self.action_size))

        if np.random.rand() < self.epsilon:
            return random.choice(valid_actions)

        q_values = self.model.predict(state[np.newaxis], verbose=0)[0]
        
        # Apply masking to q_values
        masked_q = np.full_like(q_values, -np.inf)
        for a in valid_actions:
            masked_q[a] = q_values[a]

        return np.argmax(masked_q)


    def replay(self):
        if len(self.memory) < BATCH_SIZE:
            return None
        minibatch = random.sample(self.memory, BATCH_SIZE)
        states = np.array([s for s, _, _, _, _ in minibatch])
        actions = np.array([a for _, a, _, _, _ in minibatch])
        rewards = np.array([r for _, _, r, _, _ in minibatch])
        next_states = np.array([ns for _, _, _, ns, _ in minibatch])
        dones = np.array([d for _, _, _, _, d in minibatch])

        q_values = self.model.predict(states, verbose=0)
        next_q_values = self.model.predict(next_states, verbose=0)
        next_target_q_values = self.target_model.predict(next_states, verbose=0)
        max_actions = np.argmax(next_q_values, axis=1)
        target = q_values.copy()
        for i in range(BATCH_SIZE):
            if dones[i]:
                target[i, actions[i]] = rewards[i]
            else:
                target[i, actions[i]] = rewards[i] + self.gamma * next_target_q_values[i, max_actions[i]]

        loss = self.model.train_on_batch(states, target)
        self.loss_history.append(float(loss))
        self.train_steps += 1

        if self.train_steps % TARGET_UPDATE_FREQ == 0:
            self.update_target_network()

        if self.epsilon > self.epsilon_min:
            self.epsilon *= self.epsilon_decay
        return float(loss)

    # === FULL SAVE ===
    def save(self, path_prefix):
        self.model.save(f"{path_prefix}_model.h5")
        self.target_model.save(f"{path_prefix}_target_model.h5")
        state = {
            "epsilon": self.epsilon,
            "train_steps": self.train_steps
        }
        np.save(f"{path_prefix}_meta.npy", state)

    def load(self, path_prefix):
        self.model = keras.models.load_model(f"{path_prefix}_model.h5")
        self.target_model = keras.models.load_model(f"{path_prefix}_target_model.h5")
        state = np.load(f"{path_prefix}_meta.npy", allow_pickle=True).item()
        self.epsilon = state["epsilon"]
        self.train_steps = state["train_steps"]
        self.update_target_network()

    def save_loss_history(self, path_prefix):
        np.save(f"{path_prefix}_loss.npy", np.array(self.loss_history))


class MultiAgentWrapper:
    def __init__(self, n_agents, state_dim_per_agent, action_dim):
        self.n_agents = n_agents
        self.state_dim_per_agent = state_dim_per_agent
        self.action_dim = action_dim
        self.agents = [FastDDQNAgent(state_dim_per_agent, action_dim, i) for i in range(n_agents)]

    def split_obs(self, obs):
        return [obs[i*self.state_dim_per_agent:(i+1)*self.state_dim_per_agent] for i in range(self.n_agents)]
    
    def act(self, agent_states, episode_num):
        return [agent.act(agent_states[i], episode_num) for i, agent in enumerate(self.agents)]

    def remember(self, agent_states, actions, rewards, next_agent_states, done):
        for i, agent in enumerate(self.agents):
            agent.remember(agent_states[i], actions[i], rewards[i], next_agent_states[i], done)

    def replay(self):
        for agent in self.agents:
            agent.replay()

    def save_all(self, base_path="agent"):
        for i, agent in enumerate(self.agents):
            agent.save(f"{base_path}_{i}")

    def load_all(self, base_path="agent"):
        for i, agent in enumerate(self.agents):
            agent.load(f"{base_path}_{i}")

    def save_all_losses(self, base_path="agent"):
        for i, agent in enumerate(self.agents):
            agent.save_loss_history(f"{base_path}_{i}")

    @property
    def epsilons(self):
        return [a.epsilon for a in self.agents]

def simulate_baseline_energy(env, episodes, max_steps):
    print("\n=== Running Baseline (All SBS Always Active) ===")
    baseline_energy_per_episode = []

    for ep in range(1, episodes + 1):
        obs = env.reset()
        done = False
        step_count = 0
        total_energy = 0.0

        while not done and step_count < max_steps:
            action = [0] * N_AGENTS  # Force all SBS to SM3
            _, _, done, info = env.step(np.array(action, dtype=np.uint32))

            info_str = info if isinstance(info, str) else info[0]
            info_parts = dict(item.split("=") for item in info_str.split(";"))
            energy_value = float(info_parts.get("total_energy", 0.0))
            total_energy = energy_value
            step_count += 1

        print(f"Baseline Episode {ep}: Total Energy = {total_energy:.2f} J")
        baseline_energy_per_episode.append(total_energy)

    return baseline_energy_per_episode



if __name__ == "__main__":
    env = ns3env.Ns3Env(port=5555, stepTime=0.01, startSim=True, simSeed=1)
    wrapper = MultiAgentWrapper(N_AGENTS, STATE_DIM_PER_AGENT, ACTION_DIM)

    # === Resume control ===
    RESUME = False
    if RESUME and all(os.path.exists(f"trained_ddqn_agent_{i}_model.h5") for i in range(N_AGENTS)):
        wrapper.load_all("trained_ddqn_agent")
        print("✅ Successfully loaded previous model, epsilon and steps!")
    else:
        print("🚀 Training from scratch...")

    rewards_history, avg_rewards_per_episode = [], []
    step_rewards = []
    total_energy_per_episode = []
    avg_sbs_sinr_per_episode = []
    # avg_macro_sinr_per_episode = []
    sbs_state_log = {i: [] for i in range(N_AGENTS)} 
    print("==== Fast Multi-Agent DDQN Training Start ====")
    for ep in range(1, EPISODES+1):
        obs = env.reset()
        agent_states = wrapper.split_obs(obs)
        done = False
        episode_rewards = np.zeros(N_AGENTS)
        step_count = 0
        current_episode_energy = 0.0
        sbs_sinr_sum, macro_sinr_sum, sinr_step_count = 0.0, 0.0, 0

        while not done and step_count < MAX_STEPS:
            print("--------------------")
            print(f"Step {step_count+1} (Episode {ep})")
            print("--------------------")
            actions = wrapper.act(agent_states, ep)
            next_obs, reward, done, info = env.step(np.array(actions, dtype=np.uint32))
            for i in range(N_AGENTS):
                current_state = int(next_obs[i * STATE_DIM_PER_AGENT + 1])  # 2nd value
                sbs_state_log[i].append((ep, step_count, current_state))

            info_str = info if isinstance(info, str) else info[0]
            info_parts = dict(item.split("=") for item in info_str.split(";"))
            print(f"Step {step_count+1} | Reward: {reward} | Done: {done} | Info: {info_parts} | Next Obs: {next_obs}")
            
            energy_value = float(info_parts.get("total_energy", 0.0))
            global_sinr_value = float(info_parts.get("global_sinr", 0.0))
            current_episode_energy = energy_value
            sinr_step_count += 1
            sbs_sinr_sum += global_sinr_value

            total_step_reward = reward * N_AGENTS if isinstance(reward, (float, int)) else sum(reward)
            step_rewards.append(total_step_reward)
            print(f"Step result: next_obs={next_obs}, reward={reward}, done={done}, info={info}")

            rewards = [reward] * N_AGENTS if isinstance(reward, (float, int)) else list(reward)
            next_agent_states = wrapper.split_obs(next_obs)
            wrapper.remember(agent_states, actions, rewards, next_agent_states, done)
            wrapper.replay()
            agent_states = next_agent_states
            episode_rewards += rewards
            step_count += 1

        print(f"Episode {ep}: Reward: {episode_rewards} | Epsilon: {[round(e,2) for e in wrapper.epsilons]}")
        avg_reward = np.mean(episode_rewards)
        rewards_history.append(episode_rewards)
        avg_rewards_per_episode.append(avg_reward)
        total_energy_per_episode.append(current_episode_energy)
        avg_sbs_sinr_per_episode.append(sbs_sinr_sum / sinr_step_count)
        # avg_macro_sinr_per_episode.append(macro_sinr_sum / sinr_step_count)

        packets_per_ue = SIM_TIME / PACKET_INTERVAL
        total_packets = packets_per_ue * NUM_UES
        total_bits = total_packets * PACKET_SIZE_BYTES * 8
        ee_per_episode = [total_bits / e for e in total_energy_per_episode]

        # Plot reward
        plt.figure()
        plt.plot(avg_rewards_per_episode)
        plt.xlabel("Episode")
        plt.ylabel("Reward")
        plt.title("Avg Reward Per Episode")
        plt.grid(True)
        plt.savefig("avg_reward_per_episode_live.png")
        plt.close()

        # Plot SINR
        plt.figure()
        plt.plot(avg_sbs_sinr_per_episode)
        plt.xlabel("Episode")
        plt.ylabel("Global SINR (dB)")
        plt.title("Global Avg SINR Per Episode")
        plt.grid(True)
        plt.savefig("sinr_per_episode_live.png")
        plt.close()

        # Plot Energy Efficiency
        plt.figure()
        plt.plot(ee_per_episode)
        plt.xlabel("Episode")
        plt.ylabel("Energy Efficiency (bits/Joule)")
        plt.title("Energy Efficiency Per Episode")
        plt.grid(True)
        plt.savefig("energy_efficiency_per_episode_live.png")
        plt.close()

        print(f"Episode {ep}: Reward: {episode_rewards} | Epsilon: {[round(e,2) for e in wrapper.epsilons]}")

        if ep % 100 == 0:
            wrapper.save_all("trained_ddqn_agent")
            print("Saved model snapshot at episode", ep)

    env.close()

    wrapper.save_all("trained_ddqn_agent")
    print(" Final model saved.")
    wrapper.save_all_losses("trained_ddqn_agent")
    print("Loss histories saved.")

    with open("sbs_state_history.csv", "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["episode", "step"] + [f"SBS_{i}_state" for i in range(N_AGENTS)])

        # Transpose the log structure to group by (episode, step)
        max_len = max(len(sbs_state_log[i]) for i in sbs_state_log)
        for idx in range(max_len):
            try:
                ep, step = sbs_state_log[0][idx][:2]
                row = [ep, step] + [sbs_state_log[i][idx][2] for i in range(N_AGENTS)]
                writer.writerow(row)
            except IndexError:
                continue  # skip if any SBS has fewer entries (shouldn't happen)

    # === Simulate Baseline (always active SBS) ===
    print("\n[INFO] Launching baseline simulation...")
    os.environ["NS3_BASELINE"] = "1"
    env = ns3env.Ns3Env(port=5555, stepTime=0.01, startSim=True, simSeed=2)
    baseline_energy = simulate_baseline_energy(env, EPISODES, MAX_STEPS)
    env.close()
    os.environ["NS3_BASELINE"] = "0"

    # === Save energy results ===
    np.save("rl_energy_per_episode.npy", np.array(total_energy_per_episode))
    np.save("baseline_energy_per_episode.npy", np.array(baseline_energy))

    # === Calculate energy efficiency (bits/J) ===
    ee_rl = [total_bits / e for e in total_energy_per_episode]
    ee_baseline = [total_bits / e for e in baseline_energy]

    np.save("energy_efficiency_baseline.npy", np.array(ee_baseline))

    # === Plot: Energy Consumption Comparison ===
    plt.figure(figsize=(10,6))
    plt.plot(total_energy_per_episode, label="With RL (DDQN)")
    plt.plot(baseline_energy, label="Baseline (Always Active)")
    plt.xlabel("Episode")
    plt.ylabel("Total Energy (Joules)")
    plt.title("Energy Consumption Comparison")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("energy_comparison.png")
    plt.show()

    # === Plot: Energy Efficiency Comparison ===
    plt.figure(figsize=(10,6))
    plt.plot(ee_rl, label="With RL (DDQN)")
    plt.plot(ee_baseline, label="Baseline (Always Active)")
    plt.xlabel("Episode")
    plt.ylabel("Energy Efficiency (bits/Joule)")
    plt.title("Energy Efficiency Comparison")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("energy_efficiency_comparison.png")
    plt.show()

    packets_per_ue = SIM_TIME / PACKET_INTERVAL
    total_packets = packets_per_ue * NUM_UES
    total_bits = total_packets * PACKET_SIZE_BYTES * 8
    ee_per_episode = [total_bits / energy for energy in total_energy_per_episode]


    np.save("avg_reward_per_episode.npy", np.array(avg_rewards_per_episode))
    np.save("sinr_per_episode.npy", np.array(avg_sbs_sinr_per_episode))
    np.save("energy_efficiency_per_episode.npy", np.array(ee_per_episode))
    

    plt.figure(figsize=(10,6))
    plt.plot(range(1, len(ee_per_episode)+1), ee_per_episode, marker='o')
    plt.xlabel("Episode")
    plt.ylabel("Energy Efficiency (bits/Joule)")
    plt.title("Energy Efficiency Across Episodes")
    plt.grid(True)
    plt.savefig("energy_efficiency_per_episode.png")
    plt.show()

    plt.plot(range(1, len(avg_sbs_sinr_per_episode)+1), avg_sbs_sinr_per_episode, label="Global Avg SINR")
    plt.xlabel("Episode")
    plt.ylabel("Global SINR (dB)")
    plt.title("Global Average SINR Across Episodes")
    plt.legend()
    plt.grid(True)
    plt.savefig("sinr_per_episode.png")
    plt.show()
    # plt.figure(figsize=(10, 5))
    # plt.plot(step_rewards, label="Reward per Step")
    # plt.xlabel("Step")
    # plt.ylabel("Total Reward")
    # plt.title("Step-by-Step Reward During Training")
    # plt.grid(True)
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig("reward_per_step.png")
    # plt.show()
    plt.figure(figsize=(10, 6))
    plt.plot(avg_rewards_per_episode, label='Avg Total Reward')
    plt.xlabel('Episode')
    plt.ylabel('Reward')
    plt.title('Multi-Agent DDQN Learning Progress')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("avg_reward_per_episode.png")
    plt.show()

    # === Auto-run SBS state plot script from parent folder ===
    subprocess.run(["python3", "../plot_sbs_state_times.py"])


