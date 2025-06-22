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
                          .AddTraceSource("EnergyDrained", "Traces energy consumption",
                                          MakeTraceSourceAccessor(&SmallCellEnergyModel::m_energyDrained),
                                          "ns3::TracedValueCallback::Double");
    return tid;
  }

  SmallCellEnergyModel()
    : m_currentState(ACTIVE),
      m_remainingEnergy(1000.0)
  {
    m_powerMap[ACTIVE] = 30;
    m_powerMap[SM1] = 0;
    m_powerMap[SM2] = 0;
    m_powerMap[SM3] = 0;

    m_activationDelay[SM1] = Seconds(35.5 / 1000.0);
    m_activationDelay[SM2] = MilliSeconds(0.5);
    m_activationDelay[SM3] = MilliSeconds(5);
  }

  void SetState(SmallCellState state)
  {
    m_currentState = state;
  }

  SmallCellState GetState() const
  {
    return m_currentState;
  }

  double GetPowerConsumption() const
  {
    return m_powerMap.at(m_currentState);
  }

  void UpdateEnergyConsumption(Time duration)
  {
    double power = GetPowerConsumption();
    double energyUsed = power * duration.GetSeconds();
    m_remainingEnergy -= energyUsed;
    m_energyDrained = energyUsed;

    if (m_remainingEnergy < 0)
      m_remainingEnergy = 0;

    NS_LOG_UNCOND("[" << Simulator::Now().GetSeconds() << "s] State: " << m_currentState
                      << ", Power: " << power << "W, Energy Left: " << m_remainingEnergy << "J");
  }



  
  // === Implementing all pure virtual methods from DeviceEnergyModel ===

  virtual void SetEnergySource(Ptr<EnergySource> source) override
  {
    m_energySource = source;
  }

  Ptr<EnergySource> GetEnergySource() const 
  {
    return m_energySource;
  }

  virtual double GetTotalEnergyConsumption() const override
  {
    return 1000.0 - m_remainingEnergy;
  }

  virtual void ChangeState(int newState) override
  {
    m_currentState = static_cast<SmallCellState>(newState);
  }

  virtual void HandleEnergyDepletion() override
  {
    NS_LOG_UNCOND("Energy depleted at " << Simulator::Now().GetSeconds() << "s");
  }

  virtual void HandleEnergyRecharged() override
  {
    NS_LOG_UNCOND("Energy recharged at " << Simulator::Now().GetSeconds() << "s");
  }

  virtual void HandleEnergyChanged() override
  {
    NS_LOG_UNCOND("Energy level changed at " << Simulator::Now().GetSeconds() << "s");
  }

private:
  SmallCellState m_currentState;
  std::map<SmallCellState, double> m_powerMap;
  std::map<SmallCellState, Time> m_activationDelay;
  double m_remainingEnergy;
  TracedValue<double> m_energyDrained;
  Ptr<EnergySource> m_energySource;
};

} // namespace ns3

#endif // SMALL_CELL_ENERGY_MODEL_H
