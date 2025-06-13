import numpy as np
import random
from collections import deque
import tensorflow as tf
from tensorflow import keras
from ns3gym import ns3env

DEBUG_VERBOSE = True  # Set to False if you want less frequent prints

# ---- DDQN Agent Definition ----

class DDQNAgent:
    def __init__(self, state_size, action_size, agent_id, gamma=0.99, lr=1e-3,
                 batch_size=64, mem_size=10000, tau=0.005, epsilon=1.0, epsilon_min=0.05, epsilon_decay=0.995):
        self.state_size = state_size
        self.action_size = action_size
        self.agent_id = agent_id

        self.gamma = gamma
        self.lr = lr
        self.batch_size = batch_size
        self.memory = deque(maxlen=mem_size)
        self.tau = tau

        self.epsilon = epsilon
        self.epsilon_min = epsilon_min
        self.epsilon_decay = epsilon_decay

        self.model = self._build_model()
        self.target_model = self._build_model()
        self.update_target_network()

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

    def soft_update_target_network(self):
        weights = np.array(self.model.get_weights(), dtype=object)
        target_weights = np.array(self.target_model.get_weights(), dtype=object)
        new_weights = self.tau * weights + (1 - self.tau) * target_weights
        self.target_model.set_weights(new_weights)

    def remember(self, state, action, reward, next_state, done):
        self.memory.append((state, int(action), reward, next_state, done))

    def act(self, state):
        if np.random.rand() < self.epsilon:
            chosen = random.randrange(self.action_size)
            if DEBUG_VERBOSE:
                print(f"[Agent {self.agent_id}] Random action: {chosen} (epsilon={self.epsilon:.2f})")
            return chosen
        q_values = self.model.predict(state[np.newaxis], verbose=0)
        chosen = np.argmax(q_values[0])
        if DEBUG_VERBOSE:
            print(f"[Agent {self.agent_id}] Greedy action: {chosen} (Q-values: {q_values[0]}, epsilon={self.epsilon:.2f})")
        return chosen

    def replay(self):
        if len(self.memory) < self.batch_size:
            if DEBUG_VERBOSE:
                print(f"[Agent {self.agent_id}] Not enough samples to train. {len(self.memory)}/{self.batch_size}")
            return

        indices = np.random.choice(len(self.memory), self.batch_size)
        minibatch = [self.memory[idx] for idx in indices]
        states, targets = [], []

        for state, action, reward, next_state, done in minibatch:
            q_update = reward
            if not done:
                next_action = np.argmax(self.model.predict(next_state[np.newaxis], verbose=0)[0])
                q_update += self.gamma * self.target_model.predict(next_state[np.newaxis], verbose=0)[0][next_action]
            q_values = self.model.predict(state[np.newaxis], verbose=0)[0]
            q_values[int(action)] = q_update
            states.append(state)
            targets.append(q_values)

        self.model.fit(np.array(states), np.array(targets), epochs=1, verbose=0)
        self.soft_update_target_network()
        if self.epsilon > self.epsilon_min:
            self.epsilon *= self.epsilon_decay
        if DEBUG_VERBOSE:
            print(f"[Agent {self.agent_id}] Performed replay/train step. Epsilon now: {self.epsilon:.2f}")

# ---- Multi-Agent Environment Wrapper ----

class MultiAgentNs3Env:
    def __init__(self, n_agents, env, state_dim_per_agent, action_dim):
        self.n_agents = n_agents
        self.env = env
        self.state_dim_per_agent = state_dim_per_agent
        self.action_dim = action_dim

    def split_obs(self, obs):
        return [obs[i*self.state_dim_per_agent:(i+1)*self.state_dim_per_agent] for i in range(self.n_agents)]

    def merge_actions(self, actions):
        return np.array(actions, dtype=np.uint32)

# ---- Main Training Loop ----

if __name__ == "__main__":
    PORT = 5555
    STEP_TIME = 1.0
    N_AGENTS = 3
    STATE_DIM_PER_AGENT = 3
    ACTION_DIM = 4
    EPISODES = 300

    env = ns3env.Ns3Env(port=PORT, stepTime=STEP_TIME, startSim=True, simSeed=1)
    print("Observation space:", env.observation_space)
    print("Action space:", env.action_space)
    wrapper = MultiAgentNs3Env(N_AGENTS, env, STATE_DIM_PER_AGENT, ACTION_DIM)

    agents = [DDQNAgent(STATE_DIM_PER_AGENT, ACTION_DIM, i) for i in range(N_AGENTS)]

    rewards_history = []

    print("==== Multi-Agent DDQN Training Start ====")
    for ep in range(1, EPISODES+1):
        print(f"\n--- Episode {ep} ---")
        obs = env.reset()
        agent_states = wrapper.split_obs(obs)
        done = False
        episode_rewards = np.zeros(N_AGENTS)
        t = 0

        while not done:
            print(f"\nStep {t+1}:")
            actions = np.array([agent.act(agent_states[i]) for i, agent in enumerate(agents)], dtype=np.float64)
            # actions = np.array(actions, dtype=np.uint32) 
            print(f"  Actions taken: {actions}")

            action_arr = wrapper.merge_actions(actions)
            next_obs, reward, done, info = env.step(action_arr)
        
            print(f"Step result: next_obs={next_obs}, reward={reward}, done={done}, info={info}")
            

            if isinstance(reward, (float, int)):
                rewards = [reward] * N_AGENTS
            else:
                rewards = list(reward) if hasattr(reward, '__len__') else [reward]*N_AGENTS

            next_agent_states = wrapper.split_obs(next_obs)

            for i in range(N_AGENTS):
                print(f"    [Agent {i}] State: {agent_states[i]}, Next State: {next_agent_states[i]}, Reward: {rewards[i]}")
                agents[i].remember(agent_states[i], actions[i], rewards[i], next_agent_states[i], done)

            agent_states = next_agent_states
            episode_rewards += rewards
            t += 1

            for agent in agents:
                agent.replay()

        print(f"\n>>> Episode {ep} finished in {t} steps | Rewards: {episode_rewards} | Epsilon: {[round(a.epsilon, 2) for a in agents]}")
        rewards_history.append(episode_rewards)

        if ep % 100 == 0:
            for i, agent in enumerate(agents):
                agent.model.save(f"ddqn_agent_{i}_ep{ep}.h5")
            print(f"[Checkpoint] Saved models at episode {ep}")

    env.close()

    print("\n==== Training Finished ====")

    try:
        import matplotlib.pyplot as plt
        rewards_history = np.array(rewards_history)
        plt.plot(np.mean(rewards_history, axis=1), label='Average Total Reward')
        plt.xlabel('Episode')
        plt.ylabel('Reward')
        plt.title('Multi-Agent DDQN Learning Curve')
        plt.legend()
        plt.show()
    except ImportError:
        pass