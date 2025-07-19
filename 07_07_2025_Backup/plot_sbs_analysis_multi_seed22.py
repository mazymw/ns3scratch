
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load CSV
df = pd.read_csv("test_active_ue_and_sbs_state_per_step.csv")

# === 1. State Distribution (per seed averaged) ===
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

plt.figure(figsize=(10, 6))
mean_counts.T.plot(kind="bar")
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

plt.figure(figsize=(10, 6))
sns.boxplot(x="State_Label", y="Active_UEs", data=long_df, order=["ACTIVE", "SM1", "SM2", "SM3"])
plt.title("Distribution of Active UEs by SBS State")
plt.xlabel("SBS State")
plt.ylabel("Number of Active UEs")
plt.grid(True)
plt.tight_layout()
plt.savefig("ue_count_by_state_multi_seed.png")

# === 3. State Over Time: Seed Overlay for All SBS ===
sbs_ids = [0, 1, 2]

for i in sbs_ids:
    plt.figure(figsize=(12, 6))
    for seed, group in df.groupby("Seed"):
        plt.plot(group["Step"], group[f"SBS_{i}_State"], drawstyle="steps-post", alpha=0.5, label=f"Seed {seed}")
    plt.title(f"SBS {i} State Over Time by Seed")
    plt.ylabel("State")
    plt.xlabel("Step")
    plt.yticks([0, 1, 2, 3], ["ACTIVE", "SM1", "SM2", "SM3"])
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"sbs{i}_state_over_time_by_seed.png")
