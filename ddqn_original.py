import random
import numpy as np
import tensorflow as tf
from collections import deque
from ns3gym import ns3env
import os

# # Define Q-Network using Keras
# def build_q_network(obs_dim, action_dim):
#     model = tf.keras.Sequential([
#         tf.keras.layers.Input(shape=(obs_dim,)),
#         tf.keras.layers.Dense(128, activation='relu'),
#         tf.keras.layers.Dense(128, activation='relu'),
#         tf.keras.layers.Dense(action_dim)
#     ])
#     model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3), loss='mse')
#     return model

# # DDQN Agent
# class DDQNAgent:
#     def __init__(self, obs_dim, action_dim, id):
#         self.obs_dim = obs_dim
#         self.action_dim = action_dim
#         self.id = id

#         self.gamma = 0.99
#         self.epsilon = 1.0
#         self.epsilon_decay = 0.995
#         self.epsilon_min = 0.01
#         self.batch_size = 64
#         self.memory = deque(maxlen=10000)

#         self.q_net = build_q_network(obs_dim, action_dim)
#         self.target_net = build_q_network(obs_dim, action_dim)
#         self.update_target_network()

#     def act(self, state):
#         if np.random.rand() < self.epsilon:
#             return random.randint(0, self.action_dim - 1)
#         q_values = self.q_net.predict(np.expand_dims(state, axis=0), verbose=0)
#         return np.argmax(q_values[0])

#     def memorize(self, state, action, reward, next_state, done):
#         self.memory.append((state, action, reward, next_state, done))

#     def replay(self):
#         if len(self.memory) < self.batch_size:
#             return

#         batch = random.sample(self.memory, self.batch_size)
#         states, actions, rewards, next_states, dones = zip(*batch)

#         states = np.array(states)
#         next_states = np.array(next_states)
#         rewards = np.array(rewards)
#         dones = np.array(dones).astype(float)

#         q_values = self.q_net.predict(states, verbose=0)
#         q_next = self.q_net.predict(next_states, verbose=0)
#         q_target_next = self.target_net.predict(next_states, verbose=0)

#         for i in range(self.batch_size):
#             next_action = np.argmax(q_next[i])
#             q_values[i, actions[i]] = rewards[i] + (1 - dones[i]) * self.gamma * q_target_next[i, next_action]

#         self.q_net.fit(states, q_values, epochs=1, verbose=0)

#         if self.epsilon > self.epsilon_min:
#             self.epsilon *= self.epsilon_decay

#     def update_target_network(self):
#         self.target_net.set_weights(self.q_net.get_weights())

#     def save_model(self, path):
#         self.q_net.save(path)


# Initialize ns3-gym environment
env = ns3env.Ns3Env(port=5555)
print("Environment initialized.")
obs_dim = env.observation_space.shape[0]
action_dim = env.action_space.n
num_agents = env.num_agents

agents = [DDQNAgent(obs_dim, action_dim, id=i) for i in range(num_agents)]
num_episodes = 1000

for episode in range(num_episodes):
    state = env.reset()
    done = [False] * num_agents
    total_rewards = [0 for _ in range(num_agents)]

    while not all(done):
        actions = [agents[i].act(state[i]) for i in range(num_agents)]
        next_state, reward, done, _ = env.step(actions)

        for i in range(num_agents):
            agents[i].memorize(state[i], actions[i], reward[i], next_state[i], done[i])
            agents[i].replay()
            total_rewards[i] += reward[i]

        state = next_state

    for agent in agents:
        agent.update_target_network()

    print(f"Episode {episode+1}/{num_episodes} - Total Reward: {[round(r, 2) for r in total_rewards]}")

    if (episode + 1) % 100 == 0:
        os.makedirs("checkpoints", exist_ok=True)
        for agent in agents:
            agent.save_model(f"checkpoints/ddqn_agent_{agent.id}_episode_{episode+1}")
