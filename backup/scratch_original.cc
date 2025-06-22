#include <fstream>
#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/lte-module.h"
#include "ns3/mobility-module.h"
#include "ns3/applications-module.h"
#include "ns3/config-store.h"
#include "ns3/point-to-point-module.h"
#include "ns3/netanim-module.h"
#include "smallCellEnergyModel.h"
#include "ns3/basic-energy-source.h"
#include "ns3/energy-source-container.h"
#include "ns3/device-energy-model-container.h"
#include "ns3/flow-monitor-helper.h"
// #include "smallCellEnergyModelOriginal.h"
#include "ns3/opengym-module.h"
#include <map>
#include <vector>
#include <unordered_map>
#include <cmath>

// #include "ns3/config-store.h"
// #include "ns3/string.h"
// #include "ns3/config.h"

using namespace ns3;
// using namespace ns3::energy;

// Global tracking maps
std::map<uint32_t, uint32_t> sbsToUeCount;         // SBS NodeId -> UE count
std::map<uint64_t, uint32_t> ueToSbsMap;           // IMSI -> SBS NodeId
std::map<uint16_t, uint32_t> cellIdToNodeId;       // CellId -> SBS NodeId
std::unordered_map<uint32_t, Ptr<SmallCellEnergyModel>> sbsEnergyModels;
std::map<std::pair<uint16_t, uint16_t>, uint64_t> rntiToImsiMap;

std::map<uint64_t, ApplicationContainer> ueServerApps; // Map IMSI to server apps
std::map<uint64_t, ApplicationContainer> ueClientApps; // Map IMSI to client apps

// Global map to store UE IMSI and the corresponding SBS NodeId
std::map<uint32_t, std::vector<uint64_t>> sbsToUeMap;
std::map<uint64_t, bool> ueActivityMap;

NS_LOG_COMPONENT_DEFINE("SmallCellSleep");
NS_OBJECT_ENSURE_REGISTERED(SmallCellEnergyModel);

// === Callback for connection established ===
void NotifyConnectionEstablishedUe(std::string context, uint64_t imsi,
                                   uint16_t cellId, uint16_t rnti)
{
    if (cellIdToNodeId.find(cellId) == cellIdToNodeId.end()) {
        std::cout << Simulator::Now().GetSeconds()
                  << "s: Unknown CellId " << cellId
                  << " during connection for IMSI " << imsi << std::endl;
        return;
    }

    uint32_t sbsNodeId = cellIdToNodeId[cellId];

    // Avoid duplicates
    if (ueToSbsMap.count(imsi)) {
        uint32_t prevSbs = ueToSbsMap[imsi];
        auto& list = sbsToUeMap[prevSbs];
        list.erase(std::remove(list.begin(), list.end(), imsi), list.end());
    }

    // Update maps
    sbsToUeMap[sbsNodeId].push_back(imsi);
    ueToSbsMap[imsi] = sbsNodeId;

    rntiToImsiMap[std::make_pair(cellId, rnti)] = imsi;
    std::cout << Simulator::Now().GetSeconds()
              << "s: UE with IMSI " << imsi
              << " connected to NodeId: " << sbsNodeId
              << " (CellId " << cellId << ", RNTI " << rnti << ")" << std::endl;
}


// void NotifyConnectionEstablishedUe(std::string context, uint64_t imsi, uint16_t cellId, uint16_t rnti)
// {
//     // NS_LOG_UNCOND("UE with IMSI " << imsi << " connected to CellId: " << cellId);

//     // Find the corresponding SBS NodeId based on the CellId
//     if (cellIdToNodeId.find(cellId) == cellIdToNodeId.end()) {
        
//         NS_LOG_UNCOND("Unknown CellId " << cellId);
//         return;
//     }

//     uint32_t sbsNodeId = cellIdToNodeId[cellId];
//     NS_LOG_UNCOND("UE with IMSI " << imsi << " connected to NodeId: " << sbsNodeId);
//     // Add the IMSI of the connected UE to the map for the corresponding SBS NodeId
//     sbsToUeMap[sbsNodeId].push_back(imsi);
// }



// void NotifyConnectionEstablishedUe(std::string context, uint64_t imsi, uint16_t cellId, uint16_t rnti)
// {
//     NS_LOG_UNCOND("UE with IMSI " << imsi << " connected to CellId: " << cellId);

//     if (cellIdToNodeId.find(cellId) == cellIdToNodeId.end()) {
//         NS_LOG_UNCOND("Unknown CellId " << cellId);
//         return;
//     }

//     uint32_t newSbsId = cellIdToNodeId[cellId];

//     if (ueToSbsMap.count(imsi) > 0) {
//         uint32_t oldSbsId = ueToSbsMap[imsi];
//         if (oldSbsId == newSbsId) {
//             return; // No change, avoid double-counting
//         }
//         sbsToUeCount[oldSbsId] = std::max(0U, sbsToUeCount[oldSbsId] - 1);
//     }

//     ueToSbsMap[imsi] = newSbsId;
//     sbsToUeCount[newSbsId]++;
// }

// Function to print SBS -> UE IMSI mapping
void PrintSbsUeMappings()
{
    NS_LOG_UNCOND("=== SBS -> UE IMSI Mappings at " << Simulator::Now().GetSeconds() << "s ===");

    for (const auto& entry : sbsToUeMap)
    {
        uint32_t sbsNodeId = entry.first;
        const std::vector<uint64_t>& ueList = entry.second;

        NS_LOG_UNCOND("SBS NodeId " << sbsNodeId << " has the following UE(s) connected:");
        for (uint64_t imsi : ueList)
        {
            NS_LOG_UNCOND("  - UE IMSI: " << imsi);
        }
    }

    Simulator::Schedule(Seconds(1.0), &PrintSbsUeMappings);
}

// void PrintSbsUeMappings()
// {
//     NS_LOG_UNCOND("=== SBS -> UE IMSI Mappings at " << Simulator::Now().GetSeconds() << "s ===");
    
//     // Iterate through the SBS to UE map and print the IMSIs connected to each SBS
//     for (const auto& entry : sbsToUeMap)
//     {
//         uint32_t sbsNodeId = entry.first;
//         const std::vector<uint64_t>& ueImsiList = entry.second;

//         NS_LOG_UNCOND("SBS NodeId " << sbsNodeId << " has the following UE(s) connected:");
//         for (const auto& imsi : ueImsiList)
//         {
//             NS_LOG_UNCOND("  - UE IMSI: " << imsi);
//         }
//     }

//     // Schedule the function to be called every second
//     Simulator::Schedule(Seconds(1.0), &PrintSbsUeMappings);
// }

// === Periodic print of SBS UE counts ===
// void PrintSbsConnections()
// {
//     NS_LOG_UNCOND("=== SBS Connection Summary at " << Simulator::Now().GetSeconds() << "s ===");
//     for (const auto& pair : sbsToUeCount)
//     {
//         NS_LOG_UNCOND("SBS NodeId " << pair.first << " has " << pair.second << " UE(s) connected.");
//     }
//     Simulator::Schedule(Seconds(1.0), &PrintSbsConnections);
// }
void PrintSbsConnections()
{
    std::cout << "=== SBS Connection Summary at " << Simulator::Now().GetSeconds() << "s ===" << std::endl;
    for (const auto& pair : sbsToUeMap)
    {
        uint32_t sbsNodeId = pair.first;
        const std::vector<uint64_t>& ueImsis = pair.second;

        std::cout << "SBS NodeId " << sbsNodeId << " has the following UE(s) connected:" << std::endl;
        for (uint64_t imsi : ueImsis)
        {
            std::cout << "  - UE IMSI: " << imsi << std::endl;
        }
    }
    Simulator::Schedule(Seconds(0.5), &PrintSbsConnections);
}


// === Periodic energy monitor ===
// void PrintRemainingEnergy(std::vector<Ptr<SmallCellEnergyModelOriginal>> energyModels)
// {
//     NS_LOG_UNCOND("Time: " << Simulator::Now().GetSeconds() << "s - Remaining energy:");
//     for (size_t i = 0; i < energyModels.size(); ++i)
//     {
//         NS_LOG_UNCOND("  SBS[" << i << "]: " << energyModels[i]->GetRemainingEnergy() << " J");
//     }
//     Simulator::Schedule(Seconds(1.0), &PrintRemainingEnergy, energyModels);
// }

// void PrintUePositions(NodeContainer ueNodes)
// {
//     NS_LOG_UNCOND("=== UE Positions at " << Simulator::Now().GetSeconds() << "s ===");
//     for (uint32_t i = 0; i < ueNodes.GetN(); ++i)
//     {
//         Ptr<Node> ue = ueNodes.Get(i);
//         Ptr<MobilityModel> mobility = ue->GetObject<MobilityModel>();
//         if (mobility)
//         {
//             Vector pos = mobility->GetPosition();
//             NS_LOG_UNCOND("  UE[" << i + 1 << "] Position: (" << pos.x << ", " << pos.y << ", " << pos.z << ")");
//         }
//     }
//     Simulator::Schedule(Seconds(1.0), &PrintUePositions, ueNodes);
// }

void PrintUePositions(NodeContainer ueNodes)
{
    std::cout << "=== UE Positions at " << Simulator::Now().GetSeconds() << "s ===" << std::endl;
    for (uint32_t i = 0; i < ueNodes.GetN(); ++i)
    {
        Ptr<Node> ue = ueNodes.Get(i);
        Ptr<MobilityModel> mobility = ue->GetObject<MobilityModel>();
        if (mobility)
        {
            Vector pos = mobility->GetPosition();
            std::cout << "  UE[" << i + 1 << "] Position: (" << pos.x << ", " << pos.y << ", " << pos.z << ")" << std::endl;
        }
    }

    Simulator::Schedule(Seconds(1.0), &PrintUePositions, ueNodes);
}


void PrintSbsPositions(NodeContainer smallCellEnbs)
{
    std::cout << "=== SBS Positions at " << Simulator::Now().GetSeconds() << "s ===" << std::endl;

    for (uint32_t i = 0; i < smallCellEnbs.GetN(); ++i)
    {
        Ptr<Node> sbs = smallCellEnbs.Get(i);
        uint32_t nodeId = sbs->GetId();  // Get the NodeId of the SBS
        Ptr<MobilityModel> sbsMobility = sbs->GetObject<MobilityModel>();

        if (sbsMobility)
        {
            Vector sbsPos = sbsMobility->GetPosition();
            std::cout << "  SBS NodeId " << nodeId 
                      << " Position: (" << sbsPos.x << ", " 
                      << sbsPos.y << ", " << sbsPos.z << ")" << std::endl;
        }
    }

    // Optional: reschedule if you want to monitor over time
    // Simulator::Schedule(Seconds(1.0), &PrintSbsPositions, smallCellEnbs);
}

void LogUePositions(Ptr<Node> ueNode, uint32_t imsi) {
    Ptr<MobilityModel> mob = ueNode->GetObject<MobilityModel>();
    Vector pos = mob->GetPosition();

    std::ofstream logFile("ue_positions.log", std::ios_base::app);  // Open in append mode
    if (logFile.is_open()) {
        logFile << "Time " << Simulator::Now().GetSeconds() << "s - "
                << "UE[" << imsi << "] Position: ("
                << pos.x << ", " << pos.y << ", " << pos.z << ")" << std::endl;
        logFile.close();
    }
    
    Simulator::Schedule(Seconds(1.0), &LogUePositions, ueNode, imsi);
}

void NotifyHandoverEndOkUe(std::string context, uint64_t imsi,
    uint16_t cellId, uint16_t rnti)
{
    if (cellIdToNodeId.find(cellId) == cellIdToNodeId.end()) {
    NS_LOG_UNCOND("Unknown CellId " << cellId << " during handover for IMSI " << imsi);
    return;
    }

    uint32_t newSbsNodeId = cellIdToNodeId[cellId];

    // Remove from old SBS
    if (ueToSbsMap.count(imsi)) {
    uint32_t oldSbsNodeId = ueToSbsMap[imsi];
    auto& oldList = sbsToUeMap[oldSbsNodeId];
    oldList.erase(std::remove(oldList.begin(), oldList.end(), imsi), oldList.end());
    }

    // Add to new SBS
    sbsToUeMap[newSbsNodeId].push_back(imsi);
    ueToSbsMap[imsi] = newSbsNodeId;

    std::cout << Simulator::Now().GetSeconds()
    << "s: UE IMSI " << imsi
    << " successfully handed over to NodeId " << newSbsNodeId
    << " (CellId " << cellId << ", RNTI " << rnti << ")" << std::endl;


    Simulator::Schedule(MilliSeconds(10), [=]() {
        rntiToImsiMap[std::make_pair(cellId, rnti)] = imsi;
    });
}

void UpdateEnergy(Ptr<SmallCellEnergyModel> model)
{
  Time interval = Seconds(1); // or change as needed
  model->UpdateEnergyConsumption(interval);
  Simulator::Schedule(interval, &UpdateEnergy, model);
}

// void
// NotifyHandoverEndOkUe(std::string context, uint64_t imsi, uint16_t cellId, uint16_t rnti)
// {
//     if (cellIdToNodeId.find(cellId) == cellIdToNodeId.end()) {
//         NS_LOG_UNCOND("Unknown CellId " << cellId << " during handover for IMSI " << imsi);
//         return;
//     }

//     uint32_t sbsNodeId = cellIdToNodeId[cellId];

//     std::cout << Simulator::Now().GetSeconds()
//               << "s: UE IMSI " << imsi
//               << " successfully handed over to NodeId " << sbsNodeId
//               << " (CellId " << cellId << ", RNTI " << rnti << ")" << std::endl;
// }

void PrintSbsUeCounts()
{
    std::cout << "=== SBS UE Count Summary at " << Simulator::Now().GetSeconds() << "s ===" << std::endl;

    for (const auto& entry : sbsToUeMap)
    {
        uint32_t sbsNodeId = entry.first;
        size_t ueCount = entry.second.size();

        std::cout << "SBS NodeId " << sbsNodeId << " has " << ueCount << " UE(s) connected." << std::endl;
    }

    // Schedule this function to run again in 1 second
    Simulator::Schedule(Seconds(1.0), &PrintSbsUeCounts);
}


std::string StateToString(SmallCellEnergyModel::SmallCellState state)
{
    switch (state)
    {
        case SmallCellEnergyModel::ACTIVE: return "ACTIVE";
        case SmallCellEnergyModel::SM1:    return "SM1";
        case SmallCellEnergyModel::SM2:    return "SM2";
        case SmallCellEnergyModel::SM3:    return "SM3";
        default:                           return "UNKNOWN";
    }
}

void UpdateSbsSleepStates()
{
    // uint32_t sleepThreshold = 1; // Minimum active UEs to stay awake

    for (const auto& entry : sbsToUeMap)
    {
        uint32_t sbsNodeId = entry.first;
        const auto& ueList = entry.second;
        
        if (sbsEnergyModels.find(sbsNodeId) != sbsEnergyModels.end())
        {
            Ptr<SmallCellEnergyModel> scEnergyModel = sbsEnergyModels[sbsNodeId];
            
            // Count ACTIVE UEs connected to this SBS
            uint32_t activeUeCount = 0;
            for (uint64_t imsi : ueList) {
                if (ueActivityMap[imsi]) {
                    activeUeCount++;
                }
            }

            // Set state based on active UEs
            // if (activeUeCount == 0) {
            //     scEnergyModel->SetState(SmallCellEnergyModel::SM3); // Deep sleep
            // } 
            // else if (activeUeCount <= sleepThreshold) {
            //     scEnergyModel->SetState(SmallCellEnergyModel::SM1); // Light sleep
            // }
            // else {
                scEnergyModel->SetState(SmallCellEnergyModel::ACTIVE);
            // }

            // Update SBS transmission power
            Ptr<Node> sbsNode = NodeList::GetNode(sbsNodeId);
            Ptr<LteEnbNetDevice> enbDev = sbsNode->GetDevice(0)->GetObject<LteEnbNetDevice>();
            Ptr<LteEnbPhy> enbPhy = enbDev->GetPhy();
            
            // Get power from energy model (convert W to dBm if needed)
            double txPowerDbm = scEnergyModel->GetTransmissionPower();
            // double txPowerDbm = 10 * log10(txPowerW) + 30; // Convert W to dBm
            enbPhy->SetTxPower(txPowerDbm);

            std::cout << Simulator::Now().GetSeconds()
                      << "s: SBS " << sbsNodeId << " has " << activeUeCount << " active UEs"
                      << ", State: " << StateToString(scEnergyModel->GetState())
                      << ", Tx Power: " << txPowerDbm << " dBm" << std::endl;
        }
    }

    Simulator::Schedule(Seconds(0.5), &UpdateSbsSleepStates);
}


void SinrTraceCallback(std::string context, double rsrp, double sinr)
{
    // Optional: parse IMSI from context string if needed
    std::cout << "RSRP: " << rsrp << " dBm, SINR: " << sinr << " dB" << std::endl;
}


void
ReportTraceCallback(uint16_t cellId, uint16_t rnti, double rsrp, double sinr, uint8_t componentCarrierId)
{
    uint32_t sbsNodeId = cellIdToNodeId[cellId];
    double timeNow = Simulator::Now().GetSeconds();

    auto key = std::make_pair(cellId, rnti);
    auto it = rntiToImsiMap.find(key);
    uint64_t imsi = (it != rntiToImsiMap.end()) ? it->second : 0;  // 0 = unknown

    std::cout << "Time=" << timeNow << "s"
              << ", NodeId=" << sbsNodeId
              << ", RNTI=" << rnti
              << ", IMSI=" << imsi
              << ", RSRP=" << rsrp
              << ", SINR=" << sinr
              << ", CCId=" << static_cast<uint32_t>(componentCarrierId) << std::endl;
}

void PrintRadioParameters(NetDeviceContainer enbDevs, NetDeviceContainer ueDevs)
{
    std::cout << "\n==== Effective Radio Parameters ====\n";

    // Example: Print DL Tx Power from the first eNB
    if (enbDevs.GetN() > 0)
    {
        Ptr<LteEnbNetDevice> enbDev = DynamicCast<LteEnbNetDevice>(enbDevs.Get(0));
        Ptr<LteEnbPhy> phy = enbDev->GetPhy();
        std::cout << "eNB DL Transmission Power: " << phy->GetTxPower() << " dBm" << std::endl;
    }

    // Example: Print UE noise figure from the first UE
    if (ueDevs.GetN() > 0)
    {
        Ptr<LteUeNetDevice> ueDev = DynamicCast<LteUeNetDevice>(ueDevs.Get(0));
        Ptr<LteUePhy> phy = ueDev->GetPhy();
        std::cout << "UE Noise Figure: " << phy->GetNoiseFigure() << " dB" << std::endl;
    }

    std::cout << "=====================================\n";
}

void PrintUeToSbsDistances(NodeContainer ueNodes,
                           std::map<uint16_t, uint32_t> cellIdToNodeId)
{
    for (uint32_t i = 0; i < ueNodes.GetN(); ++i)
    {
        Ptr<Node> ue = ueNodes.Get(i);
        Ptr<MobilityModel> ueMobility = ue->GetObject<MobilityModel>();
        Vector uePos = ueMobility->GetPosition();

        Ptr<NetDevice> ueDevice = ue->GetDevice(0);
        Ptr<LteUeNetDevice> lteUeDev = DynamicCast<LteUeNetDevice>(ueDevice);

        uint64_t imsi = lteUeDev->GetImsi(); 
        uint16_t cellId = lteUeDev->GetRrc()->GetCellId();

        if (cellIdToNodeId.find(cellId) != cellIdToNodeId.end())
        {
            uint32_t sbsNodeId = cellIdToNodeId[cellId];
            Ptr<Node> sbsNode = NodeList::GetNode(sbsNodeId);
            Ptr<MobilityModel> sbsMobility = sbsNode->GetObject<MobilityModel>();
            Vector sbsPos = sbsMobility->GetPosition();

            double distance = std::sqrt(std::pow(uePos.x - sbsPos.x, 2) +
                                        std::pow(uePos.y - sbsPos.y, 2) +
                                        std::pow(uePos.z - sbsPos.z, 2));

            std::cout << "UE IMSI " << imsi
                      << " connected to SBS NodeId " << sbsNodeId
                      << " (CellId " << cellId << ")"
                      << " at distance: " << distance << " m" << std::endl;
        }
        else
        {
            std::cout << "UE IMSI " << imsi
                      << " is connected to unknown CellId " << cellId << std::endl;
        }
    }
}

std::vector<Vector> GetSbsPositions(NodeContainer smallCellEnbs)
{
    std::vector<Vector> sbsPositions;

    for (uint32_t i = 0; i < smallCellEnbs.GetN(); ++i)
    {
        Ptr<Node> sbs = smallCellEnbs.Get(i);
        Ptr<MobilityModel> sbsMobility = sbs->GetObject<MobilityModel>();

        if (sbsMobility)
        {
            Vector pos = sbsMobility->GetPosition();
            sbsPositions.push_back(pos);
        }
    }

    return sbsPositions;
}

void ToggleUeApplications(uint64_t imsi, bool activate)
{
    // Get the applications for this UE
    auto serverApp = ueServerApps.find(imsi);
    auto clientApp = ueClientApps.find(imsi);
    
    if (serverApp != ueServerApps.end() && clientApp != ueClientApps.end()) {
        if (activate) {
            serverApp->second.Start(Seconds(0.1));
            clientApp->second.Start(Seconds(0.1));
        } else {
            serverApp->second.Stop(Seconds(0));
            clientApp->second.Stop(Seconds(0));
        }
    }
}

// Modify the ScheduleUeActivity function to be more comprehensive
void UpdateUeActivity(uint64_t imsi, double simulationTime) 
{
    Ptr<UniformRandomVariable> rand = CreateObject<UniformRandomVariable>();
    double nextCheck = rand->GetValue(1, 3); // Check every 1-3 seconds

    // Capture 'rand' in the lambda capture list
    Simulator::Schedule(Seconds(nextCheck), [imsi, simulationTime, rand]() {
        double now = Simulator::Now().GetSeconds();
        double simTimeHours = fmod((now / simulationTime) * 24.0, 24.0);

        double activeProbability;
        if (simTimeHours >= 0 && simTimeHours < 6) activeProbability = 0.1;
        else if (simTimeHours >= 6 && simTimeHours < 9) activeProbability = 0.4;
        else if (simTimeHours >= 9 && simTimeHours < 17) activeProbability = 0.8;
        else if (simTimeHours >= 17 && simTimeHours < 22) activeProbability = 0.5;
        else activeProbability = 0.2;

        bool newStatus = (rand->GetValue() < activeProbability);
        ueActivityMap[imsi] = newStatus;

        // Toggle applications based on new status
        ToggleUeApplications(imsi, newStatus);

        std::cout << "Time=" << now << "s (Hour=" << simTimeHours
                  << "): UE IMSI=" << imsi
                  << (newStatus ? " is now ACTIVE" : " is now INACTIVE") << std::endl;

        UpdateUeActivity(imsi, simulationTime);
    });
}

// double GetUserRequestProbability(double timeSeconds) {
//     // Simulate 1 full day in 60s: Map time 0–60s to 0–24h
//     double hour = (timeSeconds / 60.0) * 24.0;

//     if (hour >= 0 && hour < 6) return 0.01;   // Night: almost no requests
//     if (hour >= 6 && hour < 12) return 0.3;   // Morning: increasing
//     if (hour >= 12 && hour < 18) return 0.6;  // Afternoon: peak
//     if (hour >= 18 && hour < 22) return 0.4;  // Evening: still high
//     return 0.1;                               // Late evening
// }

uint32_t CountActiveUesNearSbs(Ptr<Node> sbsNode, NodeContainer ueNodes, double radius)
{
    Ptr<MobilityModel> sbsMobility = sbsNode->GetObject<MobilityModel>();
    uint32_t count = 0;

    for (uint32_t i = 0; i < ueNodes.GetN(); ++i)
    {
        Ptr<Node> ueNode = ueNodes.Get(i);
        Ptr<MobilityModel> ueMobility = ueNode->GetObject<MobilityModel>();

        Ptr<NetDevice> ueDev = ueNode->GetDevice(0);
        Ptr<LteUeNetDevice> ueLteDev = ueDev->GetObject<LteUeNetDevice>();
        uint64_t imsi = ueLteDev->GetImsi();

        if (ueActivityMap[imsi]) // active
        {
            double dist = ueMobility->GetDistanceFrom(sbsMobility);
            if (dist <= radius)
                count++;
        }
    }

    return count;
}

// Add this helper function to get active UE count per SBS
uint32_t GetActiveUeCount(uint32_t sbsNodeId)
{
    uint32_t count = 0;
    if (sbsToUeMap.find(sbsNodeId) != sbsToUeMap.end()) {
        for (uint64_t imsi : sbsToUeMap[sbsNodeId]) {
            if (ueActivityMap[imsi]) {
                count++;
            }
        }
    }
    return count;
}

// Add this to print active UE counts
void PrintActiveUeCounts()
{
    std::cout << "=== Active UE Counts at " << Simulator::Now().GetSeconds() << "s ===" << std::endl;
    for (const auto& entry : sbsToUeMap) {
        uint32_t sbsNodeId = entry.first;
        std::cout << "SBS " << sbsNodeId << ": " 
                  << GetActiveUeCount(sbsNodeId) << " active UEs" << std::endl;
    }
    Simulator::Schedule(Seconds(1.0), &PrintActiveUeCounts);
}


// === START OF RL ===

class LteGymEnv : public OpenGymEnv
{
public:
    LteGymEnv(const std::vector<Ptr<Node>>& sbsNodes,
              const std::map<uint32_t, Ptr<SmallCellEnergyModel>>& energyModels);

    // OpenGym interface implementation
    Ptr<OpenGymSpace> GetObservationSpace() override;
    Ptr<OpenGymSpace> GetActionSpace() override;
    Ptr<OpenGymDataContainer> GetObservation() override;
    float GetReward() override;
    bool GetGameOver() override;
    std::string GetExtraInfo() override;
    bool ExecuteActions(Ptr<OpenGymDataContainer> action) override;

private:
    std::vector<Ptr<Node>> m_sbsNodes;
    std::map<uint32_t, Ptr<SmallCellEnergyModel>> m_energyModels;

    std::map<uint32_t, Time> m_lastUpdateTimes;
    std::map<uint32_t, double> m_lastEnergyConsumptions;
    std::map<uint32_t, uint32_t> m_activeUeCounts;

    void CollectEnvData();
};

LteGymEnv::LteGymEnv(const std::vector<Ptr<Node>>& sbsNodes,
                     const std::map<uint32_t, Ptr<SmallCellEnergyModel>>& energyModels)
    : m_sbsNodes(sbsNodes), m_energyModels(energyModels)
{
    for (const auto& [id, model] : energyModels) {
        m_lastUpdateTimes[id] = Seconds(0);
        m_lastEnergyConsumptions[id] = 0.0;
        m_activeUeCounts[id] = 0;
    }
}

Ptr<OpenGymSpace> LteGymEnv::GetObservationSpace()
{
    uint32_t dim = m_energyModels.size() * 3;  // 3 values per SBS
    std::vector<uint32_t> shape = {dim};
    float low = 0.0;
    float high = 100.0;  // Adjust as needed

    return CreateObject<OpenGymBoxSpace>(low, high, shape, "float32");
}

void LteGymEnv::CollectEnvData()
{
    Time now = Simulator::Now();

    for (const auto& [sbsId, model] : m_energyModels) {
        Time interval = now - m_lastUpdateTimes[sbsId];
        m_lastUpdateTimes[sbsId] = now;

        // Count active UEs
        uint32_t activeUes = 0;
        if (sbsToUeMap.find(sbsId) != sbsToUeMap.end()) {
            for (uint64_t imsi : sbsToUeMap[sbsId]) {
                if (ueActivityMap[imsi]) {
                    activeUes++;
                }
            }
        }
        m_activeUeCounts[sbsId] = activeUes;

        // Energy consumption
        m_lastEnergyConsumptions[sbsId] = model->GetTotalEnergyConsumption();
    }
}

Ptr<OpenGymDataContainer> LteGymEnv::GetObservation()
{
    CollectEnvData();
    std::vector<double> obs;

    for (const auto& [sbsId, model] : m_energyModels) {
        obs.push_back((double)m_activeUeCounts[sbsId]);
        obs.push_back((double)model->GetState());
        obs.push_back(m_lastEnergyConsumptions[sbsId]);
    }

    Ptr<OpenGymBoxContainer<double>> container = CreateObject<OpenGymBoxContainer<double>>(std::vector<uint32_t>{(uint32_t)obs.size()});
    container->SetData(obs);
    return container;
}

float LteGymEnv::GetReward()
{
    CollectEnvData();
    float totalReward = 0.0;

    for (const auto& [sbsId, model] : m_energyModels) {
        double energyPenalty = -m_lastEnergyConsumptions[sbsId] * 0.1;
        double qosReward = m_activeUeCounts[sbsId] * 5.0;
        totalReward += energyPenalty + qosReward;
    }

    return totalReward;
}

bool LteGymEnv::ExecuteActions(Ptr<OpenGymDataContainer> action)
{
    Ptr<OpenGymBoxContainer<uint32_t>> box = DynamicCast<OpenGymBoxContainer<uint32_t>>(action);
    if (!box) return false;

    std::vector<uint32_t> actions = box->GetData();
    if (actions.size() != m_energyModels.size()) return false;

    size_t i = 0;
    for (const auto& [sbsId, model] : m_energyModels) {
        model->SetState((SmallCellEnergyModel::SmallCellState)actions[i]);
        ++i;
    }

    return true;
}

// LteGymEnv::~LteGymEnv()
// {
//     // Optional: cleanup code
// }



Ptr<OpenGymSpace> LteGymEnv::GetActionSpace()
{
    uint32_t numStates = 4;  // e.g., 4 power states
    std::vector<uint32_t> shape = {1};
    return CreateObject<OpenGymDiscreteSpace>(numStates);
}

bool LteGymEnv::GetGameOver()
{
    return false;  // Or custom logic to stop the episode
}

std::string LteGymEnv::GetExtraInfo()
{
    return "Extra info if needed";
}





// === END OF RL ===

// === Main Simulation ===
int main(int argc, char *argv[])
{
    uint16_t numSmallCells = 3;
    uint16_t numUesPerCell = 3;
    uint16_t simulationTime = 15;
    uint32_t openGymPort = 5555;
    uint32_t simSeed = 0;
    CommandLine cmd;
    cmd.AddValue("openGymPort", "OpenGym port", openGymPort);
    cmd.AddValue("simSeed", "Simulation seed", simSeed);
    cmd.Parse(argc, argv);

    NS_LOG_UNCOND("Simulation is running.");
    

    Ptr<LteHelper> lteHelper = CreateObject<LteHelper>();
    Ptr<PointToPointEpcHelper> epcHelper = CreateObject<PointToPointEpcHelper>();
    lteHelper->SetEpcHelper(epcHelper);

    lteHelper->SetAttribute("PathlossModel", StringValue("ns3::LogDistancePropagationLossModel"));
    lteHelper->SetHandoverAlgorithmType("ns3::A3RsrpHandoverAlgorithm");
    lteHelper->SetHandoverAlgorithmAttribute("Hysteresis", DoubleValue(1.0));
    lteHelper->SetHandoverAlgorithmAttribute("TimeToTrigger", TimeValue(MilliSeconds(128)));
    lteHelper->SetAttribute("FadingModel", StringValue("ns3::TraceFadingLossModel"));

    Config::SetDefault("ns3::LteUePhy::RsrpSinrSamplePeriod", UintegerValue(1000));

    Config::SetDefault("ns3::LogDistancePropagationLossModel::Exponent", DoubleValue(3.5));
    Config::SetDefault("ns3::TraceFadingLossModel::TraceFilename",
                   StringValue("src/lte/model/fading-traces/fading_trace_EVA_60kmph.fad"));
    Config::SetDefault("ns3::TraceFadingLossModel::WindowSize", TimeValue(Seconds(0.5)));
    Config::SetDefault("ns3::TraceFadingLossModel::RbNum", UintegerValue(100));
    // Config::SetDefault("ns3::LteSpectrumPhy::NoiseFigure", DoubleValue(9.0));
    // Config::SetDefault("ns3::TraceFadingLossModel::SubBands", UintegerValue(1));
    Ptr<Node> pgw = epcHelper->GetPgwNode();

    NodeContainer macroEnb;
    macroEnb.Create(1);

    NodeContainer smallCellEnbs;
    smallCellEnbs.Create(numSmallCells);

    NodeContainer ueNodes;
    ueNodes.Create(numUesPerCell * numSmallCells);

    // === Mobility Setup ===
    MobilityHelper mobility;
    // MobilityHelper ueMobility;
    // Set static position for the first UE
    // Ptr<ListPositionAllocator> staticAlloc = CreateObject<ListPositionAllocator>();
    // staticAlloc->Add(Vector(250.0, 244.0, 0.0));
    // MobilityHelper staticMobility;
    // staticMobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    // staticMobility.SetPositionAllocator(staticAlloc);
    // staticMobility.Install(ueNodes.Get(0));  

    // NodeContainer mobileUes;
    // for (uint32_t i = 1; i < ueNodes.GetN(); ++i) {
    //     mobileUes.Add(ueNodes.Get(i));
    // }
    // Create Position Allocator separately
    Ptr<PositionAllocator> uePositionAlloc = CreateObject<RandomRectanglePositionAllocator>();
    uePositionAlloc->SetAttribute("X", StringValue("ns3::UniformRandomVariable[Min=0.0|Max=780.0]"));
    uePositionAlloc->SetAttribute("Y", StringValue("ns3::UniformRandomVariable[Min=0.0|Max=240.0]"));

    // Configure the UE MobilityHelper with RandomWaypointMobilityModel
    MobilityHelper ueMobility;
    ueMobility.SetMobilityModel("ns3::RandomWaypointMobilityModel",
        "Speed", StringValue("ns3::UniformRandomVariable[Min=1.0|Max=50.0]"),
        "Pause", StringValue("ns3::ConstantRandomVariable[Constant=0.5]"),
        "PositionAllocator", PointerValue(uePositionAlloc));
    ueMobility.SetPositionAllocator(uePositionAlloc);
    ueMobility.Install(ueNodes);

    for (uint32_t i = 0; i < ueNodes.GetN(); ++i) {
        Ptr<Node> ueNode = ueNodes.Get(i);
        uint32_t imsi = i + 1;  // adjust if needed
        LogUePositions(ueNode, imsi);
    }


    

    // Macro eNB
    mobility.SetPositionAllocator("ns3::GridPositionAllocator",
        "MinX", DoubleValue(390.0), "MinY", DoubleValue(0.0),
        "DeltaX", DoubleValue(20.0), "GridWidth", UintegerValue(1),
        "LayoutType", StringValue("RowFirst"));
    mobility.Install(macroEnb);

    // Small cells
    mobility.SetPositionAllocator("ns3::GridPositionAllocator",
        "MinX", DoubleValue(140.0), "MinY", DoubleValue(100.0),
        "DeltaX", DoubleValue(250.0), "GridWidth", UintegerValue(numSmallCells),
        "LayoutType", StringValue("RowFirst"));
    mobility.Install(smallCellEnbs);

    // === LTE Device Installation ===
    NetDeviceContainer macroEnbLteDevs = lteHelper->InstallEnbDevice(macroEnb);
    NetDeviceContainer smallCellLteDevs = lteHelper->InstallEnbDevice(smallCellEnbs);
    NetDeviceContainer ueLteDevs = lteHelper->InstallUeDevice(ueNodes);

    NodeContainer allEnbs;
    allEnbs.Add(macroEnb);
    allEnbs.Add(smallCellEnbs);
    lteHelper->AddX2Interface(allEnbs); // Add X2 interface for all eNBs

    PrintRadioParameters(smallCellLteDevs, ueLteDevs);
    for (uint32_t i = 0; i < smallCellEnbs.GetN(); ++i)
    {
        Ptr<Node> scNode = smallCellEnbs.Get(i);

        // Create and install energy source
        Ptr<BasicEnergySource> energySource = CreateObject<BasicEnergySource>();
        energySource->SetInitialEnergy(1000.0);
        scNode->AggregateObject(energySource);

        // Create and install your custom energy model
        Ptr<SmallCellEnergyModel> scEnergyModel = CreateObject<SmallCellEnergyModel>();
        scEnergyModel->SetEnergySource(energySource);
        energySource->AppendDeviceEnergyModel(scEnergyModel);


        uint32_t nodeId = scNode->GetId();
        sbsEnergyModels[nodeId] = scEnergyModel;
        
        Simulator::Schedule(Seconds(1.0), &UpdateEnergy, scEnergyModel);
    }
    
    // Build CellId -> NodeId map
    // for (uint16_t i = 0; i < smallCellLteDevs.GetN(); ++i)
    // {
    //     Ptr<LteEnbNetDevice> enbDev = smallCellLteDevs.Get(i)->GetObject<LteEnbNetDevice>();
    //     if (enbDev)
    //     {
    //         uint16_t cellId = enbDev->GetCellId();
    //         uint32_t nodeId = smallCellEnbs.Get(i)->GetId();
    //         cellIdToNodeId[cellId] = nodeId;
    //     }
    // }
    // Register ALL eNBs (macro + small cell) in the cellIdToNodeId map
    for (uint16_t i = 0; i < allEnbs.GetN(); ++i)
    {
        Ptr<Node> enbNode = allEnbs.Get(i);
        Ptr<NetDevice> dev = enbNode->GetDevice(0); // Assume LTE device is at index 0
        Ptr<LteEnbNetDevice> enbDev = dev->GetObject<LteEnbNetDevice>();
        if (enbDev)
        {
            uint16_t cellId = enbDev->GetCellId();
            uint32_t nodeId = enbNode->GetId();
            cellIdToNodeId[cellId] = nodeId;
        }
    }


    // for (const auto& pair : cellIdToNodeId)
    // {
    //     std::cout << "Mapped CellId " << pair.first << " to NodeId " << pair.second << std::endl;
    // }

    // === Internet stack ===
    InternetStackHelper internet;
    internet.Install(ueNodes);
    internet.Install(macroEnb);
    internet.Install(smallCellEnbs);
    Ipv4InterfaceContainer ueIpIface = epcHelper->AssignUeIpv4Address(NetDeviceContainer(ueLteDevs));

    // // === Initial attachment (round-robin) ===
    // for (uint16_t i = 0; i < ueNodes.GetN(); ++i)
    // {
    //     Ptr<NetDevice> ueDev = ueLteDevs.Get(i);
    //     uint16_t smallCellIndex = i % numSmallCells;
    //     Ptr<NetDevice> enbDev = smallCellLteDevs.Get(smallCellIndex);
    
    //     lteHelper->Attach(ueDev, enbDev);
    // }

    for (uint16_t i = 0; i < ueNodes.GetN(); ++i)
    {
        Ptr<Node> ueNode = ueNodes.Get(i);
        Ptr<MobilityModel> ueMobility = ueNode->GetObject<MobilityModel>();

        double minDistance = std::numeric_limits<double>::max();
        Ptr<NetDevice> closestEnbDev = nullptr;

        for (uint16_t j = 0; j < smallCellEnbs.GetN(); ++j)
        {
            Ptr<Node> enbNode = smallCellEnbs.Get(j);
            Ptr<MobilityModel> enbMobility = enbNode->GetObject<MobilityModel>();

            double distance = ueMobility->GetDistanceFrom(enbMobility);
            if (distance < minDistance)
            {
                minDistance = distance;
                closestEnbDev = smallCellLteDevs.Get(j);
            }
        }

        lteHelper->Attach(ueLteDevs.Get(i), closestEnbDev);
    }


    for (uint32_t i = 0; i < ueLteDevs.GetN(); ++i)
    {
        Ptr<LteUeNetDevice> ueDevice = ueLteDevs.Get(i)->GetObject<LteUeNetDevice>();
        Ptr<LteUePhy> uePhy = ueDevice->GetPhy();
        
        
        uePhy->TraceConnectWithoutContext("ReportCurrentCellRsrpSinr", MakeCallback(&ReportTraceCallback));
    }
   

    // === Attach callback for handover + connection tracking ===
    for (uint32_t i = 0; i < ueLteDevs.GetN(); ++i)
    {
        Ptr<NetDevice> dev = ueLteDevs.Get(i);
        Ptr<LteUeNetDevice> ueLte = dev->GetObject<LteUeNetDevice>();
        Ptr<LteUeRrc> rrc = ueLte->GetRrc();
        rrc->TraceConnect("ConnectionEstablished", "", MakeCallback(&NotifyConnectionEstablishedUe));
        rrc->TraceConnect("HandoverEndOk", "", MakeCallback(&NotifyHandoverEndOkUe));

    }

    

    for (uint32_t i = 0; i < ueNodes.GetN(); ++i)
    {
        Ptr<NetDevice> ueDev = ueLteDevs.Get(i);
        Ptr<LteUeNetDevice> ueLteDev = ueDev->GetObject<LteUeNetDevice>();
        uint64_t imsi = ueLteDev->GetImsi();

        ueActivityMap[imsi] = false; // default to inactive
        
        UpdateUeActivity(imsi, simulationTime);
    
        // Log initial position
        LogUePositions(ueNodes.Get(i), imsi);
    }

    uint16_t port = 5000;
    ApplicationContainer serverApps, clientApps;


    for (uint32_t i = 0; i < ueNodes.GetN(); ++i)
    {
        Ptr<LteUeNetDevice> ueLteDev = ueLteDevs.Get(i)->GetObject<LteUeNetDevice>();
        uint64_t imsi = ueLteDev->GetImsi();

        // UDP Server
        UdpServerHelper udpServer(port);
        ueServerApps[imsi] = udpServer.Install(ueNodes.Get(i));
        ueServerApps[imsi].Start(Seconds(1.0));
        ueServerApps[imsi].Stop(Seconds(simulationTime));

        // UDP Client
        UdpClientHelper udpClient(ueIpIface.GetAddress(i), port);
        udpClient.SetAttribute("MaxPackets", UintegerValue(32000));
        udpClient.SetAttribute("Interval", TimeValue(MilliSeconds(50)));
        udpClient.SetAttribute("PacketSize", UintegerValue(1024));
        ueClientApps[imsi] = udpClient.Install(pgw);
        ueClientApps[imsi].Start(Seconds(1.1)); // Start slightly after server
        ueClientApps[imsi].Stop(Seconds(simulationTime));
    }

    serverApps.Start(Seconds(1.0));
    clientApps.Start(Seconds(1.0));
    serverApps.Stop(Seconds(simulationTime));
    clientApps.Stop(Seconds(simulationTime));



    FlowMonitorHelper flowmonHelper;
    Ptr<FlowMonitor> monitor = flowmonHelper.InstallAll();  

    // === Schedule Logging ===
    // Simulator::Schedule(Seconds(1.0), &PrintRemainingEnergy, energyModels); 
    Simulator::Schedule(Seconds(1.0), &PrintSbsPositions, smallCellEnbs);
    Simulator::Schedule(Seconds(1.0), &PrintSbsPositions, macroEnb);
    Simulator::Schedule(Seconds(0.5), &PrintSbsConnections);
    Simulator::Schedule(Seconds(1.0), &PrintSbsUeMappings);
    Simulator::Schedule(Seconds(1.0), &PrintSbsUeCounts);
    Simulator::Schedule(Seconds(0.5), &UpdateSbsSleepStates);
    // Simulator::Schedule(Seconds(1.0), &PrintUeToSbsDistances, ueNodes, cellIdToNodeId);
    // Simulator::Schedule(Seconds(1.0), &PrintUePositions, ueNodes);


    // OPENGYM INTERFACE //
    std::vector<Ptr<Node>> sbsNodes;
    std::map<uint32_t, Ptr<SmallCellEnergyModel>> energyModels;

    for (uint32_t i = 0; i < smallCellEnbs.GetN(); ++i) {
        Ptr<Node> sbsNode = smallCellEnbs.Get(i);
        uint32_t nodeId = sbsNode->GetId();
        Ptr<SmallCellEnergyModel> model = sbsEnergyModels[nodeId];

        sbsNodes.push_back(sbsNode);
        energyModels[nodeId] = model;
    }

    // Create single environment and interface
    Ptr<LteGymEnv> gymEnv = CreateObject<LteGymEnv>(sbsNodes, energyModels);
    Ptr<OpenGymInterface> gymInterface = CreateObject<OpenGymInterface>(openGymPort);
    gymInterface->Notify(gymEnv);
    


    Simulator::Stop(Seconds(simulationTime));
    Simulator::Run();
    

    monitor->CheckForLostPackets();
    monitor->SerializeToXmlFile("flowmon-results.xml", true, true); 
    Simulator::Destroy();

    return 0;
}