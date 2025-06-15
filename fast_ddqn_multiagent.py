import numpy as np
import random
from collections import deque
import tensorflow as tf
from tensorflow import keras
from ns3gym import ns3env
import matplotlib.pyplot as plt

# Hyperparameters
BATCH_SIZE = 64
TARGET_UPDATE_FREQ = 100  # Hard copy target weights every N replay calls
MEM_SIZE = 30000
EPISODES = 100
MAX_STEPS = 1000
N_AGENTS = 3
STATE_DIM_PER_AGENT = 5
ACTION_DIM = 4

# === Environment parameters ===
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
        self.epsilon_decay =  0.99996946 # Decay rate for epsilon
        self.lr = 	0.0005

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

    def act(self, state):
        if np.random.rand() < self.epsilon:
            return random.randrange(self.action_size)
        q_values = self.model.predict(state[np.newaxis], verbose=0)
        return np.argmax(q_values[0])

    def replay(self):
        if len(self.memory) < BATCH_SIZE:
            return None
        minibatch = random.sample(self.memory, BATCH_SIZE)
        states = np.array([s for s, _, _, _, _ in minibatch])
        actions = np.array([a for _, a, _, _, _ in minibatch])
        rewards = np.array([r for _, _, r, _, _ in minibatch])
        next_states = np.array([ns for _, _, _, ns, _ in minibatch])
        dones = np.array([d for _, _, _, _, d in minibatch])

        # Vectorized Double DQN update
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

        # Hard update target network every N steps
        if self.train_steps % TARGET_UPDATE_FREQ == 0:
            self.update_target_network()

        # Epsilon decay
        if self.epsilon > self.epsilon_min:
            self.epsilon *= self.epsilon_decay
        return float(loss)

# ---- Multi-Agent Wrapper ----
class MultiAgentWrapper:
    def __init__(self, n_agents, state_dim_per_agent, action_dim):
        self.n_agents = n_agents
        self.state_dim_per_agent = state_dim_per_agent
        self.action_dim = action_dim
        self.agents = [FastDDQNAgent(state_dim_per_agent, action_dim, i) for i in range(n_agents)]

    def split_obs(self, obs):
        # Returns a list: each agent's state slice
        return [obs[i*self.state_dim_per_agent:(i+1)*self.state_dim_per_agent] for i in range(self.n_agents)]

    def act(self, agent_states):
        # Returns a list of actions, one per agent
        return [agent.act(agent_states[i]) for i, agent in enumerate(self.agents)]

    def remember(self, agent_states, actions, rewards, next_agent_states, done):
        for i, agent in enumerate(self.agents):
            agent.remember(agent_states[i], actions[i], rewards[i], next_agent_states[i], done)

    def replay(self):
        # Replay every step for all agents
        for agent in self.agents:
            agent.replay()

    @property
    def epsilons(self):
        return [a.epsilon for a in self.agents]

# ---- Main Training Loop ----

if __name__ == "__main__":
    env = ns3env.Ns3Env(port=5555, stepTime=0.01, startSim=True, simSeed=1)
    wrapper = MultiAgentWrapper(N_AGENTS, STATE_DIM_PER_AGENT, ACTION_DIM)
    rewards_history = []
    avg_rewards_per_episode = []
    step_rewards = []
    total_energy_per_episode = []
    avg_sbs_sinr_per_episode = []
    avg_macro_sinr_per_episode = []


    print("==== Fast Multi-Agent DDQN Training Start ====")
    for ep in range(1, EPISODES+1):
        obs = env.reset()
        agent_states = wrapper.split_obs(obs)
        done = False
        episode_rewards = np.zeros(N_AGENTS)
        step_count = 0
        current_episode_energy = 0.0
        sbs_sinr_sum = 0.0
        macro_sinr_sum = 0.0
        sinr_step_count = 0
        while not done and step_count < MAX_STEPS:
            print("--------------------")
            print(f"Step {step_count+1} (Episode {ep})")
            print("--------------------")
            actions = wrapper.act(agent_states)
            next_obs, reward, done, info = env.step(np.array(actions, dtype=np.uint32))
            # Parse ExtraInfo string
            info_str = info if isinstance(info, str) else info[0]
            info_parts = dict(item.split("=") for item in info_str.split(";"))

            # Energy parsing
            energy_value = float(info_parts.get("total_energy", 0.0))
            current_episode_energy = energy_value

            # SINR parsing
            sbs_sinr_value = float(info_parts.get("avg_sbs_sinr", 0.0))
            macro_sinr_value = float(info_parts.get("macro_sinr", 0.0))

            # Accumulate for this episode
            sbs_sinr_sum += sbs_sinr_value
            macro_sinr_sum += macro_sinr_value
            sinr_step_count += 1

            if isinstance(reward, (float, int)):
                total_step_reward = reward * N_AGENTS
            else:
                total_step_reward = sum(reward)

            step_rewards.append(total_step_reward)

            print(f"Step result: next_obs={next_obs}, reward={reward}, done={done}, info={info}")
            if isinstance(reward, (float, int)):
                rewards = [reward] * N_AGENTS
            else:
                rewards = list(reward) if hasattr(reward, '__len__') else [reward]*N_AGENTS
            next_agent_states = wrapper.split_obs(next_obs)
            wrapper.remember(agent_states, actions, rewards, next_agent_states, done)
            wrapper.replay()  # <-- Replay every step
            agent_states = next_agent_states
            episode_rewards += rewards
            step_count += 1
        print(f"Episode {ep}: Reward: {episode_rewards} | Epsilon: {[round(e,2) for e in wrapper.epsilons]}")
        rewards_history.append(episode_rewards)
        avg_reward = np.mean(episode_rewards)
        avg_rewards_per_episode.append(avg_reward)
        total_energy_per_episode.append(current_episode_energy)
        avg_sbs_sinr_per_episode.append(sbs_sinr_sum / sinr_step_count)
        avg_macro_sinr_per_episode.append(macro_sinr_sum / sinr_step_count)
        
        print(f"Episode {ep}: Reward: {episode_rewards} | Epsilon: {[round(e,2) for e in wrapper.epsilons]}")
        # Save model checkpoints every 100 episodes (optional)
        if ep % 100 == 0:
            for i, agent in enumerate(wrapper.agents):
                agent.model.save(f"fast_ddqn_agent_{i}_ep{ep}.h5")
    env.close()


    for i, agent in enumerate(wrapper.agents):
        agent.model.save(f"trained_ddqn_agent_{i}_exp10.h5")


     # === EE Calculation ===
    packets_per_ue = SIM_TIME / PACKET_INTERVAL
    total_packets = packets_per_ue * NUM_UES
    total_bits = total_packets * PACKET_SIZE_BYTES * 8  # bits

    ee_per_episode = [total_bits / energy for energy in total_energy_per_episode]

    # === EE Plot ===
    plt.figure(figsize=(10,6))
    plt.plot(range(1, len(ee_per_episode)+1), ee_per_episode, marker='o')
    plt.xlabel("Episode")
    plt.ylabel("Energy Efficiency (bits/Joule)")
    plt.title("Energy Efficiency Across Episodes")
    plt.grid(True)
    plt.savefig("energy_efficiency_per_episode.png")
    plt.show()

    plt.figure(figsize=(10,6))
    plt.plot(range(1, len(avg_sbs_sinr_per_episode)+1), avg_sbs_sinr_per_episode, label="SBS SINR")
    plt.plot(range(1, len(avg_macro_sinr_per_episode)+1), avg_macro_sinr_per_episode, label="Macro SINR")
    plt.xlabel("Episode")
    plt.ylabel("SINR (dB)")
    plt.title("Average SINR Across Episodes")
    plt.legend()
    plt.grid(True)
    plt.savefig("sinr_per_episode.png")
    plt.show()


    plt.figure(figsize=(10, 5))
    plt.plot(step_rewards, label="Reward per Step")
    plt.xlabel("Step")
    plt.ylabel("Total Reward")
    plt.title("Step-by-Step Reward During Training")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("reward_per_step.png")  # ✅ Save to file
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot(avg_rewards_per_episode, label='Avg Total Reward')
    plt.xlabel('Episode')
    plt.ylabel('Reward')
    plt.title('Multi-Agent DDQN Learning Progress')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("avg_reward_per_episode.png")  # ✅ Save to file
    plt.show()
