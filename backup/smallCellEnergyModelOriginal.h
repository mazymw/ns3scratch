#ifndef SMALL_CELL_ENERGY_MODEL_H
#define SMALL_CELL_ENERGY_MODEL_H

#include "ns3/core-module.h"
#include "ns3/energy-module.h"
#include "ns3/object.h"
#include "ns3/log.h"
#include "ns3/nstime.h"

namespace ns3 {

class SmallCellEnergyModel : public Object
{
public:
  // TypeId declaration for SmallCellEnergyModel
  static TypeId GetTypeId(void)
  {
    static TypeId tid = TypeId("ns3::SmallCellEnergyModel")
                          .SetParent<Object>()
                          .AddConstructor<SmallCellEnergyModel>()
                          .AddAttribute("State", "Current state of the small cell",
                                        StringValue("active"),
                                        MakeStringAccessor(&SmallCellEnergyModel::m_currentState),
                                        MakeStringChecker());
    return tid;
  }

  // Constructor: Set default values
  SmallCellEnergyModel()
    : m_currentState("active"),
      m_activePower(0.5),
      m_sleepState1Power(0.1),
      m_sleepState2Power(0.02),
      m_remainingEnergy(100.0) // initial energy in joules
  {}

  // Set the state of the small cell
  void SetState(const std::string& state)
  {
    m_currentState = state;
  }

  // Set initial energy
  void SetInitialEnergy(double energy)
  {
    m_remainingEnergy = energy;
  }

  // Get the remaining energy
  double GetRemainingEnergy() const
  {
    return m_remainingEnergy;
  }

  // Get the power consumption depending on the state
  double GetPowerConsumption() const
  {
    if (m_currentState == "active")
      return m_activePower;
    else if (m_currentState == "sleep1")
      return m_sleepState1Power;
    else if (m_currentState == "sleep2")
      return m_sleepState2Power;
    return 0.0;
  }

  // Update the energy consumption over a period of time
  void UpdateEnergyConsumption(Time duration)
  {
    double power = GetPowerConsumption(); // watts
    double energyUsed = power * duration.GetSeconds(); // joules
    m_remainingEnergy -= energyUsed;

    NS_LOG_UNCOND("[" << Simulator::Now().GetSeconds()
                      << "s] State: " << m_currentState
                      << ", Power: " << power
                      << "W, Energy Left: " << m_remainingEnergy << "J");
  }

private:
  std::string m_currentState; // Current state of the small cell (active, sleep1, sleep2)
  double m_activePower; // Power consumption in active state
  double m_sleepState1Power; // Power consumption in sleep state 1
  double m_sleepState2Power; // Power consumption in sleep state 2
  double m_remainingEnergy; // Remaining energy in joules
};

} // namespace ns3

#endif // SMALL_CELL_ENERGY_MODEL_H
