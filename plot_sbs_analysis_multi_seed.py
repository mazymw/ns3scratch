import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load CSV
df = pd.read_csv("test_active_ue_and_sbs_state_per_step.csv")

# Extract sorted unique seed values
SEEDS = sorted(df["Seed"].unique().tolist())

# Define bright RGB colors for first 3 seeds
bright_seed_colors = [
    (1.0, 0.0, 0.0),  # Bright Red
    (0.0, 1.0, 0.0),  # Bright Green
    (0.0, 0.0, 1.0),  # Bright Blue
]

# Assign colors to seeds dynamically
seed_colors = {
    seed: bright_seed_colors[i] if i < len(bright_seed_colors) else "black"
    for i, seed in enumerate(SEEDS)
}

# === Bright colors for SBS and state labels ===
bright_colors = {
    "SBS_0": (1.0, 0.4, 0.4),  # Bright Red
    "SBS_1": (0.4, 1.0, 0.4),  # Bright Green
    "SBS_2": (0.4, 0.6, 1.0),  # Bright Blue
    "ACTIVE": (1.0, 0.0, 0.0),
    "SM1":    (0.0, 1.0, 0.0),
    "SM2":    (0.0, 0.0, 1.0),
    "SM3":    (1.0, 1.0, 0.0),
}

# === 1. State Distribution (averaged across seeds) ===
sbs_states = ["SBS_0_State", "SBS_1_State", "SBS_2_State"]
all_counts = []

for sbs in sbs_states:
    temp = df.groupby(["Seed", sbs]).size().unstack(fill_value=0)
    temp.columns.name = "State"
    temp = temp.rename(columns={0: "ACTIVE", 1: "SM1", 2: "SM2", 3: "SM3"})
    temp["SBS"] = sbs
    all_counts.append(temp.reset_index())

count_df = pd.concat(all_counts)
mean_counts = count_df.groupby("SBS")[["ACTIVE", "SM1", "SM2", "SM3"]].mean()

sbs_color_list = [bright_colors["SBS_0"], bright_colors["SBS_1"], bright_colors["SBS_2"]]

plt.figure(figsize=(10, 6))
ax = mean_counts.T.plot(kind="bar", color=sbs_color_list)
plt.title("Average SBS State Distribution Over Seeds")
plt.ylabel("Avg Number of Steps")
plt.xlabel("State")
plt.grid(axis="y")
plt.xticks(rotation=0)
plt.legend(title="Base Station")
plt.tight_layout()
plt.savefig("sbs_state_distribution_avg.png")

# === 2. Boxplot: Active UEs by SBS state ===
long_df = pd.melt(
    df,
    id_vars=["Seed", "Step", "Active_UEs"],
    value_vars=sbs_states,
    var_name="SBS",
    value_name="State"
)
state_map = {0: "ACTIVE", 1: "SM1", 2: "SM2", 3: "SM3"}
long_df["State_Label"] = long_df["State"].map(state_map)

box_palette = {state: bright_colors[state] for state in ["ACTIVE", "SM1", "SM2", "SM3"]}

plt.figure(figsize=(10, 6))
sns.boxplot(
    x="State_Label", y="Active_UEs", data=long_df,
    order=["ACTIVE", "SM1", "SM2", "SM3"],
    palette=box_palette
)
plt.title("Distribution of Active UEs by SBS State")
plt.xlabel("SBS State")
plt.ylabel("Number of Active UEs")
plt.grid(True)
plt.tight_layout()
plt.savefig("ue_count_by_state_multi_seed.png")

# === 3. State Over Time: One File per SBS per Seed ===
sbs_ids = [0, 1, 2]

for seed in SEEDS:
    seed_group = df[df["Seed"] == seed]
    for i in sbs_ids:
        plt.figure(figsize=(10, 5))
        hours = np.linspace(0, 24, len(seed_group))
        plt.plot(hours, seed_group[f"SBS_{i}_State"],
                 drawstyle="steps-post",
                 color=seed_colors.get(seed, "black"),
                 label=f"Seed {seed}")
        plt.title(f"SBS {i} State Over Time (Seed {seed})")
        plt.ylabel("State")
        plt.xlabel("Hour of Day")
        plt.yticks([0, 1, 2, 3], ["ACTIVE", "SM1", "SM2", "SM3"])
        plt.gca().invert_yaxis()
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"sbs{i}_state_seed{seed}.png")
        plt.close()

# === 4. Combined Plot: All Seeds per SBS ===
for i in sbs_ids:
    plt.figure(figsize=(12, 6))
    for seed in SEEDS:
        seed_group = df[df["Seed"] == seed]
        hours = np.linspace(0, 24, len(seed_group))
        plt.plot(hours, seed_group[f"SBS_{i}_State"],
                 drawstyle="steps-post",
                 color=seed_colors.get(seed, "black"),
                 label=f"Seed {seed}",
                 alpha=0.7)
    plt.title(f"SBS {i} State Over Time (All Seeds)")
    plt.ylabel("State")
    plt.xlabel("Hour of Day")
    plt.yticks([0, 1, 2, 3], ["ACTIVE", "SM1", "SM2", "SM3"])
    plt.gca().invert_yaxis()
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"sbs{i}_state_all_seeds.png")
    plt.close()
