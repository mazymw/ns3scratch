#ifndef SMALL_CELL_ENERGY_MODEL_H
#define SMALL_CELL_ENERGY_MODEL_H

#include "ns3/core-module.h"
#include "ns3/energy-module.h"
#include "ns3/object.h"
#include "ns3/log.h"
#include "ns3/nstime.h"
#include "ns3/energy-source.h"
#include "ns3/device-energy-model.h"
#include "ns3/simulator.h"
#include <map>
#include <string>

// void LogSinrDuringTransition(uint32_t sbsNodeId, std::string label, ns3::Time endTime);

namespace ns3 {

class SmallCellEnergyModel : public DeviceEnergyModel
{
public:
  enum SmallCellState
  {
    ACTIVE = 0,
    SM1 = 1,
    SM2 = 2,
    SM3 = 3
  };

  static TypeId GetTypeId(void)
  {
    static TypeId tid = TypeId("ns3::SmallCellEnergyModel")
                          .SetParent<DeviceEnergyModel>()
                          .AddConstructor<SmallCellEnergyModel>()
                          .AddAttribute("InitialEnergy", "Initial energy in joules",
                                        DoubleValue(1000.0),
                                        MakeDoubleAccessor(&SmallCellEnergyModel::m_remainingEnergy),
                                        MakeDoubleChecker<double>())
                          .AddAttribute("ActiveTxPower", "Transmission power in ACTIVE state (W)",
                                        DoubleValue(30.0),
                                        MakeDoubleAccessor(&SmallCellEnergyModel::m_activeTxPower),
                                        MakeDoubleChecker<double>())
                          .AddTraceSource("TotalPower", "Total power consumption",
                                        MakeTraceSourceAccessor(&SmallCellEnergyModel::m_totalPower),
                                        "ns3::TracedValueCallback::Double")
                          .AddTraceSource("TxPower", "Transmission power",
                                        MakeTraceSourceAccessor(&SmallCellEnergyModel::m_txPower),
                                        "ns3::TracedValueCallback::Double")
                          .AddTraceSource("EnergyConsumed", "Total energy consumed",
                                        MakeTraceSourceAccessor(&SmallCellEnergyModel::m_energyConsumedTrace),
                                        "ns3::TracedValueCallback::Double");
    return tid;
  }

  SmallCellEnergyModel()
    : m_currentState(ACTIVE),
      m_remainingEnergy(1000.0),
      m_energyConsumed(0.0),
      m_activeTxPower(33.0),
      m_transitioning(false),
      m_transitionEndTime(Seconds(0))
  {
    // Initialize power consumption maps (total power, tx power)
    m_totalPowerMap[ACTIVE] = 20.7;    
    m_totalPowerMap[SM1] = 8.22;     
    m_totalPowerMap[SM2] = 3.55;      
    m_totalPowerMap[SM3] = 3.36;        
    
    // Transmission power for each state
    m_txPowerMap[ACTIVE] = m_activeTxPower;  // Uses attribute value
    m_txPowerMap[SM1] = 0;
    m_txPowerMap[SM2] = 0;
    m_txPowerMap[SM3] = 0;

    // State transition delays (keeping your original values)
    m_activationDelay[SM1] = Seconds(0.00355);
    m_activationDelay[SM2] = Seconds(0.05);
    m_activationDelay[SM3] = Seconds(0.5);
    
    // Initialize traced values
    m_totalPower = m_totalPowerMap[ACTIVE];
    m_txPower = m_txPowerMap[ACTIVE];
    m_energyConsumedTrace = 0.0;
    m_activeEnteredTime = Seconds(0);
  }

  void SetNodeId(uint32_t id) { m_nodeId = id; }
  uint32_t GetNodeId() const { return m_nodeId; }

void SetState(SmallCellState state)
{
    // If already transitioning, ignore further state change requests
    if (m_transitioning) {
        std::cout << Simulator::Now().GetSeconds() << "s: [SBS " << m_nodeId
                  << "] Ignoring SetState(" << StateToString(state)
                  << ") because SBS is still transitioning!" << std::endl;
        return;
    }

    // Prevent switching out of ACTIVE before 1.5s after transition to ACTIVE
    const double minActiveDuration = 0.1; // seconds
    if (m_currentState == ACTIVE && state != ACTIVE) {
        Time now = Simulator::Now();
        Time timeInActive = now - m_activeEnteredTime;
        if (timeInActive.GetSeconds() < minActiveDuration) {
            double remaining = minActiveDuration - timeInActive.GetSeconds();
            std::cout << now.GetSeconds() << "s: [SBS " << m_nodeId
                      << "] Refusing to leave ACTIVE: Only "
                      << timeInActive.GetSeconds() << "s elapsed, need "
                      << minActiveDuration << "s. Will retry in "
                      << remaining << "s." << std::endl;
            return;
        }
    }

    if (m_currentState != state) {
        std::cout << Simulator::Now().GetSeconds() << "s: [SBS " << m_nodeId
                  << "] Changing state from "
                  << StateToString(m_currentState) << " to " << StateToString(state) << std::endl;

        if (state == ACTIVE && m_currentState != ACTIVE) {
            std::cout << Simulator::Now().GetSeconds() << "s: Scheduling activation delay of "
                      << m_activationDelay[m_currentState].GetSeconds()
                      << "s before entering ACTIVE." << std::endl;

            // === Added SINR Logging ===
            // Time transitionEnd = Simulator::Now() + m_activationDelay[m_currentState];
            // Simulator::ScheduleNow(&LogSinrDuringTransition, m_nodeId, "During Transition", transitionEnd + MilliSeconds(2));

            m_transitioning = true;
            m_transitionEndTime = Simulator::Now() + m_activationDelay[m_currentState];
            Simulator::Schedule(m_activationDelay[m_currentState],
                                &SmallCellEnergyModel::FinishTransition, this);
        } else {
            m_currentState = state;
            UpdatePowerValues();
            m_transitioning = false;
            m_transitionEndTime = Seconds(0);
            if (state == ACTIVE) {
                m_activeEnteredTime = Simulator::Now();
            }
        }
    }
}


  void FinishTransition()
  {
      std::cout << Simulator::Now().GetSeconds() << "s: Finished transition, now ACTIVE." << std::endl;
      m_currentState = ACTIVE;
      UpdatePowerValues();
      m_transitioning = false;
      m_transitionEndTime = Seconds(0);
      // Now SBS is truly ACTIVE, record the time
      m_activeEnteredTime = Simulator::Now();
}

  SmallCellState GetState() const
  {
    return m_currentState;
  }

  // Returns true if in the middle of a transition to ACTIVE
  bool IsTransitioning() const
  {
    return m_transitioning;
  }

  // Returns the simulation time when the current transition will finish
  Time GetTransitionEndTime() const
  {
    return m_transitionEndTime;
  }

  double GetTotalPowerConsumption() const
  {
    return m_totalPower;
  }

  double GetTransmissionPower() const
  {
    return m_txPower;
  }

  void UpdateEnergyConsumption(Time duration)
  {
    double energyUsed = m_totalPower * duration.GetSeconds();
    m_remainingEnergy -= energyUsed;
    m_energyConsumed += energyUsed;
    m_energyConsumedTrace = m_energyConsumed;

    if (m_remainingEnergy < 0)
    {
      m_remainingEnergy = 0;
      HandleEnergyDepletion();
    }

    NS_LOG_INFO("[" << Simulator::Now().GetSeconds() << "s] State: " << StateToString(m_currentState)
                  << ", Total Power: " << m_totalPower << "W"
                  << ", Tx Power: " << m_txPower << "W"
                  << ", Energy Left: " << m_remainingEnergy << "J");
  }

  // === DeviceEnergyModel overrides ===
  virtual void SetEnergySource(Ptr<EnergySource> source) override
  {
    m_energySource = source;
  }

  virtual double GetTotalEnergyConsumption() const override
  {
    return m_energyConsumed;
  }

  virtual void ChangeState(int newState) override
  {
    SetState(static_cast<SmallCellState>(newState));
  }

  virtual void HandleEnergyDepletion() override
  {
    NS_LOG_WARN("Energy depleted at " << Simulator::Now().GetSeconds() << "s");
    SetState(SM3);  // Go to deepest sleep state when depleted
  }

  virtual void HandleEnergyRecharged() override
  {
    NS_LOG_INFO("Energy recharged at " << Simulator::Now().GetSeconds() << "s");
    SetState(ACTIVE);
  }

  virtual void HandleEnergyChanged() override
  {
    NS_LOG_DEBUG("Energy level changed at " << Simulator::Now().GetSeconds() << "s");
  }

private:
  void UpdatePowerValues()
  {
    m_totalPower = m_totalPowerMap[m_currentState];
    m_txPower = m_txPowerMap[m_currentState];
  }

  static std::string StateToString(SmallCellState state)
  {
    switch (state)
    {
      case ACTIVE: return "ACTIVE";
      case SM1: return "SM1";
      case SM2: return "SM2";
      case SM3: return "SM3";
      default: return "UNKNOWN";
    }
  }

  SmallCellState m_currentState;
  std::map<SmallCellState, double> m_totalPowerMap;
  std::map<SmallCellState, double> m_txPowerMap;
  std::map<SmallCellState, Time> m_activationDelay;
  
  // Member variables
  double m_remainingEnergy;
  double m_energyConsumed;          // Regular double for internal tracking
  double m_activeTxPower;
  uint32_t m_nodeId;
  TracedValue<double> m_totalPower; // Traced value for total power
  TracedValue<double> m_txPower;    // Traced value for transmission power
  TracedValue<double> m_energyConsumedTrace; // Traced value for energy consumed
  Ptr<EnergySource> m_energySource;

  // --- New transition state variables ---
  bool m_transitioning;
  Time m_transitionEndTime;
  Time m_activeEnteredTime;
};

} // namespace ns3

#endif // SMALL_CELL_ENERGY_MODEL_H
