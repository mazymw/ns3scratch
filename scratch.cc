
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
#include <functional>
#include <stdexcept>
#include <numeric>
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
std::map<Ptr<UdpServer>, Ptr<Node>> udpServerToUeNode;
std::set<uint32_t> macroNodeIds;
NetDeviceContainer* globalEnbDevs;
std::map<uint32_t, double> sbsSinrAverage;
double macroSinrAverage = 0.0;
NetDeviceContainer* globalUeDevs;

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

static std::ofstream&
GetRsrpSinrLog(const std::string& filename = "rsrp-sinr-log.csv",
               bool               append    = false)
{
    static std::ofstream logFile;
    static bool          init = false;

    if (!init)
    {
        std::ios_base::openmode mode = std::ios::out |
                                       (append ? std::ios::app
                                               : std::ios::trunc);

        logFile.open(filename, mode);
        if (!logFile.is_open())
        {
            throw std::runtime_error("Cannot open log file: " + filename);
        }

        if (!append)      // header once
            logFile << "time_s,nodeId,rnti,imsi,rsrp_w,sinr_lin,ccId\n";

        init = true;
    }
    return logFile;
}

static std::ofstream&
GetDebugLog(const std::string& filename = "sinr-debug.log",
            bool               append    = false)
{
    static std::ofstream logFile;
    static bool          init = false;

    if (!init)
    {
        std::ios_base::openmode mode = std::ios::out |
                                       (append ? std::ios::app
                                               : std::ios::trunc);

        logFile.open(filename, mode);
        if (!logFile.is_open())
            throw std::runtime_error("Cannot open debug log: " + filename);

        init = true;
    }
    return logFile;          // same instance on every call
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

    // std::cout << Simulator::Now().GetSeconds()
    // << "s: UE IMSI " << imsi
    // << " successfully handed over to NodeId " << newSbsNodeId
    // << " (CellId " << cellId << ", RNTI " << rnti << ")" << std::endl;


    Simulator::Schedule(MilliSeconds(10), [=]() {
        rntiToImsiMap[std::make_pair(cellId, rnti)] = imsi;
    });
}

// void UpdateEnergy(Ptr<SmallCellEnergyModel> model)
// {
//   Time interval = Seconds(1); // or change as needed
//   model->UpdateEnergyConsumption(interval);
//   Simulator::Schedule(interval, &UpdateEnergy, model);
// }

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


// void
// ReportTraceCallback(uint16_t cellId, uint16_t rnti, double rsrp, double sinr, uint8_t componentCarrierId)
// {
//     uint32_t sbsNodeId = cellIdToNodeId[cellId];
//     double timeNow = Simulator::Now().GetSeconds();

//     auto key = std::make_pair(cellId, rnti);
//     auto it = rntiToImsiMap.find(key);
//     uint64_t imsi = (it != rntiToImsiMap.end()) ? it->second : 0;  // 0 = unknown

//     std::cout << "Time=" << timeNow << "s"
//               << ", NodeId=" << sbsNodeId
//               << ", RNTI=" << rnti
//               << ", IMSI=" << imsi
//               << ", RSRP=" << rsrp
//               << ", SINR=" << sinr
//               << ", CCId=" << static_cast<uint32_t>(componentCarrierId) << std::endl;
// }

void ReportTraceCallback(uint16_t cellId,
                    uint16_t rnti,
                    double   rsrp,   // linear W
                    double   sinr,   // linear ratio
                    uint8_t  componentCarrierId)
{
    uint32_t nodeId = cellIdToNodeId[cellId];
    double   t      = Simulator::Now().GetSeconds();

    auto key     = std::make_pair(cellId, rnti);
    auto it      = rntiToImsiMap.find(key);
    uint64_t imsi = (it != rntiToImsiMap.end()) ? it->second : 0;

    GetRsrpSinrLog()               // <-- log entry
        << t        << ','
        << nodeId   << ','
        << rnti     << ','
        << imsi     << ','
        << rsrp     << ','
        << sinr     << ','
        << uint32_t(componentCarrierId)
        << '\n';
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
    double nextCheck = rand->GetValue(0.1, 0.3); // Randomly check every 0.1 to 0.3 seconds

    // Capture 'rand' in the lambda capture list
    Simulator::Schedule(Seconds(nextCheck), [imsi, simulationTime, rand]() {
        double now = Simulator::Now().GetSeconds();
        double simTimeHours = fmod((now / simulationTime) * 24.0, 24.0);
        // std::cout << "SimTime: " << now << "s  (Hour of day: " << simTimeHours << ")" << std::endl;
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

        // std::cout << "Time=" << now << "s (Hour=" << simTimeHours
        //           << "): UE IMSI=" << imsi
        //           << (newStatus ? " is now ACTIVE" : " is now INACTIVE") << std::endl;

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

// Utility: dBm to Watt
double dBmToWatt(double dBm) {
    return std::pow(10.0, (dBm - 30.0) / 10.0);
}

// Utility: Watt to dB
double WattTodB(double watt) {
    return 10.0 * std::log10(watt);
}

// Path loss model (simple log-distance; tweak as needed)
double CalculateRxPowerDbm(double txPowerDbm, double distance, double pathLossExp = 3.5) {
    if (distance < 1.0) distance = 1.0;
    return txPowerDbm - 10 * pathLossExp * std::log10(distance);
}



bool IsMacroNode(uint32_t nodeId) {
    extern std::set<uint32_t> macroNodeIds;
    return macroNodeIds.find(nodeId) != macroNodeIds.end();
}


// double CalculateSinr(
//     Ptr<Node> ueNode,
//     Ptr<Node> candidateEnbNode,
//     double    candidateTxPowerDbm,
//     double    pathLossExponent,
//     NetDeviceContainer& enbDevs,
//     bool enableLog = true)  // <-- new param
// {
//     auto&  dbg   = GetDebugLog();
//     double tNow  = Simulator::Now().GetSeconds();

//     Ptr<MobilityModel> ueMob = ueNode->GetObject<MobilityModel>();
//     Ptr<MobilityModel> candidateEnbMob = candidateEnbNode->GetObject<MobilityModel>();
//     if (enableLog) {
//         Vector uePos = ueMob->GetPosition();
//         Vector enbPos = candidateEnbMob->GetPosition();

//         // dbg << "    [DEBUG] UE Position:  (" << uePos.x << ", " << uePos.y << ", " << uePos.z << ")\n";
//         // dbg << "    [DEBUG] eNB Position: (" << enbPos.x << ", " << enbPos.y << ", " << enbPos.z << ")\n";
//     }
//     double distance = ueMob->GetDistanceFrom(candidateEnbMob);

//     double rxSignalDbm  = CalculateRxPowerDbm(candidateTxPowerDbm, distance,  pathLossExponent);
//     double rxSignalWatt = dBmToWatt(rxSignalDbm);

//     double interferenceWatt = 0.0;
//     for (uint32_t i = 0; i < enbDevs.GetN(); ++i) {
//         Ptr<LteEnbNetDevice> enbDev = DynamicCast<LteEnbNetDevice>(enbDevs.Get(i));
//         if (!enbDev) continue;

//         Ptr<Node> enbNode = enbDev->GetNode();
//         uint32_t enbNodeId = enbNode->GetId();

//         if (enbNodeId == candidateEnbNode->GetId()) continue; // skip self

//         if (IsMacroNode(enbNodeId)) {
//             // dbg << "    [DEBUG] Skipping interference from macro eNB NodeId=" << enbNodeId << "\n";
//             continue;
//         }// <<< SKIP MACRO NODES
//         Ptr<MobilityModel> otherEnbMob = enbDev->GetNode()->GetObject<MobilityModel>();
//         double otherDist = ueMob->GetDistanceFrom(otherEnbMob);

//         Ptr<LteEnbPhy> enbPhy = enbDev->GetPhy();
//         double otherTxPowerDbm = enbPhy ? enbPhy->GetTxPower() : 44.77;
//         double rxInterfDbm = CalculateRxPowerDbm(otherTxPowerDbm, otherDist,  pathLossExponent);
//         double rxInterfWatt = dBmToWatt(rxInterfDbm);

//         interferenceWatt += rxInterfWatt;

//         if (enableLog) {
//             Vector interfererPos = otherEnbMob->GetPosition();

//             dbg << "    [" << tNow << "s]  Interferer eNB NodeId=" << enbNodeId
//                 << " | Tx=" << otherTxPowerDbm << " dBm"
//                 << " | Pos=(" << interfererPos.x << ", " << interfererPos.y << ", " << interfererPos.z << ")"
//                 << " | Dist=" << otherDist << " m"
//                 << " | Rx=" << rxInterfDbm << " dBm (" << rxInterfWatt << " W)\n";
//         }
//     }

//     double noiseWatt = 1.38e-23 * 290 * 20e6;
//     double denom = interferenceWatt + noiseWatt;
//     if (denom <= 0.0) denom = 1e-12;
//     double sinr = rxSignalWatt / denom;
//     double sinrDb = 10.0 * std::log10(sinr);

//     uint64_t imsi = 0;
//     Ptr<NetDevice> ueNetDev = ueNode->GetDevice(0);
//     Ptr<LteUeNetDevice> ueLteDev = ueNetDev->GetObject<LteUeNetDevice>();
//     if (ueLteDev) {
//         imsi = ueLteDev->GetImsi();
//     }

//     if (enableLog) {
//         dbg << "\n[" << tNow << "s]  Calculating SINR  UE(IMSI=" << imsi
//             << ")  Candidate eNB(NodeId=" << candidateEnbNode->GetId() << ")\n";
//         dbg << "    UE pos: " << ueMob->GetPosition()
//             << " | eNB pos: " << candidateEnbMob->GetPosition()
//             << " | Distance: " << distance << " m\n";
//         dbg << "    Candidate TxPower: " << candidateTxPowerDbm << " dBm\n";
//         dbg << "    [" << tNow << "s]  Signal Power from Serving eNB: "
//             << rxSignalDbm << " dBm (" << rxSignalWatt << " W)\n";
//         dbg << "    [" << tNow << "s]  Total Interference: " << interferenceWatt << " W\n";
//         dbg << "    [" << tNow << "s]  Noise Power: " << noiseWatt << " W\n";
//         dbg << "    [" << tNow << "s]  SINR: " << sinr << " (linear)  = "
//             << sinrDb << " dB\n";
//     }

//     return sinrDb;
// }

double CalculateSinr(
    Ptr<Node> ueNode,
    Ptr<Node> candidateEnbNode,
    double candidateTxPowerDbm,
    double pathLossExponent,
    NetDeviceContainer& enbDevs,
    bool enableLog = true  // <-- still exists to maintain function signature
) {
    Ptr<MobilityModel> ueMob = ueNode->GetObject<MobilityModel>();
    Ptr<MobilityModel> candidateEnbMob = candidateEnbNode->GetObject<MobilityModel>();
    double distance = ueMob->GetDistanceFrom(candidateEnbMob);

    double rxSignalDbm  = CalculateRxPowerDbm(candidateTxPowerDbm, distance, pathLossExponent);
    double rxSignalWatt = dBmToWatt(rxSignalDbm);

    double interferenceWatt = 0.0;
    int activeSmallCellCount = 0;
    for (uint32_t i = 0; i < enbDevs.GetN(); ++i) {
        Ptr<LteEnbNetDevice> enbDev = DynamicCast<LteEnbNetDevice>(enbDevs.Get(i));
        if (!enbDev) continue;

        Ptr<Node> enbNode = enbDev->GetNode();
        uint32_t enbNodeId = enbNode->GetId();

        if (enbNodeId == candidateEnbNode->GetId()) continue;
        if (IsMacroNode(enbNodeId)) continue;

        Ptr<MobilityModel> otherEnbMob = enbNode->GetObject<MobilityModel>();
        double otherDist = ueMob->GetDistanceFrom(otherEnbMob);

        Ptr<LteEnbPhy> enbPhy = enbDev->GetPhy();
        double otherTxPowerDbm = enbPhy ? enbPhy->GetTxPower() : 44.77;

        if (otherTxPowerDbm > 0)  // treat anything higher than "off" threshold as active
            ++activeSmallCellCount;


        double rxInterfDbm = CalculateRxPowerDbm(otherTxPowerDbm, otherDist, pathLossExponent);
        double rxInterfWatt = dBmToWatt(rxInterfDbm);

        interferenceWatt += rxInterfWatt;
    }

        // Impose minimum interference if fewer than 2 active small cells
    if (activeSmallCellCount <= 1 && interferenceWatt < 1e-9) {
        interferenceWatt = 1e-9;
        // if (enableLog) {
                   
        //     std::cout << "[INFO] Enforcing interference floor: 1e-9 W due to low active SBS count (" 
        //               << activeSmallCellCount << ")" << std::endl;
        // }
    }

    double noiseWatt = 1.38e-23 * 290 * 20e6;
    double denom = interferenceWatt + noiseWatt;
    if (denom <= 0.0) denom = 1e-12;

    double sinr = rxSignalWatt / denom;
    double sinrDb = 10.0 * std::log10(sinr);
    return sinrDb;
}





void CustomSinrHandover(
    Ptr<LteHelper> lteHelper, 
    NodeContainer& ueNodes, 
    NetDeviceContainer& enbDevs, 
    NetDeviceContainer& ueDevs, 
    double pathLossExponent, 
    double sinrHysteresisDb = 2.0 // Minimum SINR improvement needed (dB)
) {
    NS_ASSERT(lteHelper != nullptr);

    if (enbDevs.GetN() == 0 || ueDevs.GetN() == 0) {
        std::cout << "[WARNING] Empty device containers, skipping SINR-based handover decision" << std::endl;
        return;
    }

    for (uint32_t i = 0; i < ueDevs.GetN(); ++i) {
        // std::cout << "\n[DEBUG] Processing UE " << i << std::endl;

        Ptr<LteUeNetDevice> ueDevice = DynamicCast<LteUeNetDevice>(ueDevs.Get(i));
        if (!ueDevice) continue;

        Ptr<LteUeRrc> ueRrc = ueDevice->GetRrc();
        if (!ueRrc || ueRrc->GetState() != LteUeRrc::CONNECTED_NORMALLY) continue;

        Ptr<Node> ueNode = ueDevice->GetNode();
        uint32_t ueNodeId = ueNode->GetId();

        // Find serving eNB (by NodeId, not CellId)
        Ptr<LteEnbNetDevice> servingEnb = nullptr;
        double servingSinrDb = -std::numeric_limits<double>::infinity();
        uint32_t servingEnbNodeId = 0;

        Ptr<LteEnbNetDevice> bestEnbDev = nullptr;
        double bestSinrDb = -std::numeric_limits<double>::infinity();
        uint32_t bestEnbNodeId = 0;

        Ptr<LteEnbNetDevice> targetEnbDevice = nullptr;
        uint32_t targetEnbNodeId = 0;

        // Find serving eNB
        for (uint32_t j = 0; j < enbDevs.GetN(); ++j) {
            Ptr<LteEnbNetDevice> enbDevice = DynamicCast<LteEnbNetDevice>(enbDevs.Get(j));
            if (!enbDevice) continue;
            Ptr<Node> enbNode = enbDevice->GetNode();
            if (ueRrc->GetCellId() == enbDevice->GetCellId()) { // can add a map for robustness
                servingEnb = enbDevice;
                servingEnbNodeId = enbNode->GetId();
                Ptr<LteEnbPhy> enbPhy = enbDevice->GetPhy();
                double txPowerDbm = enbPhy ? enbPhy->GetTxPower() : 44.77;
                servingSinrDb = CalculateSinr(ueNode, enbNode, txPowerDbm, pathLossExponent, enbDevs, true); 
                bestSinrDb = servingSinrDb;
                bestEnbDev = servingEnb;
                bestEnbNodeId = servingEnbNodeId;
                // std::cout << "[DEBUG] Found serving eNB (NodeId: " << servingEnbNodeId << ", CellId: " << enbDevice->GetCellId()
                //           << ") SINR: " << servingSinrDb << " dB" << std::endl;
                break;
            }
        }
        if (!servingEnb) continue;

        // Check all neighbor eNBs by NodeId
        for (uint32_t j = 0; j < enbDevs.GetN(); ++j) {
            Ptr<LteEnbNetDevice> enbDevice = DynamicCast<LteEnbNetDevice>(enbDevs.Get(j));
            if (!enbDevice) continue;
            Ptr<Node> enbNode = enbDevice->GetNode();
            uint32_t enbNodeId = enbNode->GetId();
            if (enbNodeId == servingEnbNodeId) continue; // skip serving eNB

            Ptr<LteEnbPhy> enbPhy = enbDevice->GetPhy();
            double txPowerDbm = enbPhy ? enbPhy->GetTxPower() : 44.77;
            double neighborSinrDb = CalculateSinr(ueNode, enbNode, txPowerDbm, pathLossExponent, enbDevs, false); 

            // std::cout << "[DEBUG] Neighbor eNB (NodeId: " << enbNodeId << ", CellId: " << enbDevice->GetCellId()
            //           << ") SINR: " << neighborSinrDb << " dB" << std::endl;

            if (neighborSinrDb > bestSinrDb) {
                bestSinrDb = neighborSinrDb;
                bestEnbDev = enbDevice;
                bestEnbNodeId = enbNodeId;
            }
        }

        uint64_t imsi = ueDevice->GetObject<LteUeNetDevice>()->GetImsi();
        // Handover decision: compare NodeIds
        // std::cout << "\n[Handover Check] UE NodeId=" << ueNodeId << ", IMSI=" << ueRrc->GetUeDevice()->GetImsi() << std::endl;
        // std::cout << "  Serving eNB NodeId: " << servingEnbNodeId << ", SINR: " << servingSinrDb << " dB" << std::endl;
        // std::cout << "  Best Candidate eNB NodeId: " << bestEnbNodeId << ", SINR: " << bestSinrDb << " dB" << std::endl;
        // std::cout << "  SINR Improvement: " << (bestSinrDb - servingSinrDb) << " dB (Threshold: " << sinrHysteresisDb << " dB)" << std::endl;

        // if (bestEnbNodeId != servingEnbNodeId) {
        //     std::cout << "  ✅ Candidate is different from current serving eNB.\n";
        // } else {
        //     std::cout << "  ❌ Candidate is the same as serving eNB. No handover.\n";
        // }

        // if ((bestSinrDb - servingSinrDb) > sinrHysteresisDb) {
        //     std::cout << "  ✅ SINR improvement exceeds hysteresis threshold.\n";
        // } else {
        //     std::cout << "  ❌ SINR improvement too small for handover.\n";
        // }

        // if (ueRrc->GetState() == LteUeRrc::CONNECTED_NORMALLY) {
        //     std::cout << "  ✅ UE is CONNECTED_NORMALLY.\n";
        // } else {
        //     std::cout << "  ❌ UE is NOT in CONNECTED_NORMALLY state.\n";
        // }


        // if (imsi != 0) {
        //     std::cout << "  ✅ UE IMSI is valid (" << imsi << ").\n";
        // } else {
        //     std::cout << "  ❌ UE IMSI is 0. Handover not possible.\n";
        // }


        if (bestEnbNodeId != servingEnbNodeId &&
            (bestSinrDb - servingSinrDb) > sinrHysteresisDb &&
            ueRrc->GetState() == LteUeRrc::CONNECTED_NORMALLY &&
            ueRrc->GetRnti() != 0)
        {
            // std::cout << "  🚀 Triggering SINR-based handover from NodeId " << servingEnbNodeId
            //         << " to NodeId " << bestEnbNodeId << "\n" << std::endl;
            lteHelper->HandoverRequest(Seconds(0), ueDevice, servingEnb, bestEnbDev);
        }
    }

    // Schedule the next handover decision
    const double HANDOVER_CHECK_INTERVAL = 0.001;
    Simulator::Schedule(
        Seconds(HANDOVER_CHECK_INTERVAL),
        &CustomSinrHandover,
        lteHelper, ueNodes, enbDevs, ueDevs, pathLossExponent, sinrHysteresisDb
    );
    // std::cout << "\n[INFO] SINR-based handover decision completed. Next check scheduled in "
    //           << HANDOVER_CHECK_INTERVAL << " seconds.\n" << std::endl;
}


void
InitializeSbsTxPower ()
{
  for (const auto &entry : sbsEnergyModels)
    {
      uint32_t                    nodeId = entry.first;
      Ptr<SmallCellEnergyModel>   model  = entry.second;

      // 1. Get the power the model thinks it should radiate (W)
      double txPowerW   = model->GetTransmissionPower ();        // 30 W default
      // 2. Convert W → dBm for LteEnbPhy::SetTxPower()
      double txPowerDbm = (txPowerW > 0)
                            ? txPowerW
                            : 0;                           // “off” flag

      // 3. Push it into the PHY
      Ptr<Node>            sbsNode = NodeList::GetNode (nodeId);
      Ptr<LteEnbNetDevice> enbDev  =
          sbsNode->GetDevice (0)->GetObject<LteEnbNetDevice> ();
      Ptr<LteEnbPhy>       enbPhy  = enbDev->GetPhy ();
      enbPhy->SetTxPower (txPowerDbm);

      std::cout << "0s: [SBS " << nodeId << "] Initial TxPower set to "
                << txPowerDbm << " dBm"  << std::endl;
    }
}

static std::ofstream&
GetSbsSinrLog(const std::string& filename = "sbs_avg_sinr.csv", bool append = false)
{
    static std::ofstream logFile;
    static bool initialized = false;

    if (!initialized) {
        std::ios_base::openmode mode = std::ios::out | (append ? std::ios::app : std::ios::trunc);
        logFile.open(filename, mode);
        if (!logFile.is_open()) {
            throw std::runtime_error("Cannot open SINR log file: " + filename);
        }

        // Write CSV header
        logFile << "time_s,sbs_id,avg_sinr_db,valid_ue_count\n";
        initialized = true;
    }

    return logFile;
}

void SampleSbsSinr(NetDeviceContainer& enbDevs, NetDeviceContainer& ueDevs, double pathLossExponent) {
    double tNow = Simulator::Now().GetSeconds();
    // auto& logFile = GetSbsSinrLog();

    for (const auto& [sbsId, ueList] : sbsToUeMap) {
        double sum = 0.0;
        uint32_t count = 0;

        Ptr<Node> sbsNode = NodeList::GetNode(sbsId);
        if (!sbsNode) {
            continue;
        }

        Ptr<LteEnbNetDevice> enbDev = sbsNode->GetDevice(0)->GetObject<LteEnbNetDevice>();
        double txPowerDbm = enbDev->GetPhy()->GetTxPower();

        for (uint64_t imsi : ueList) {
            //Only consider active UEs
            if (!ueActivityMap[imsi]) continue;

            for (uint32_t i = 0; i < ueDevs.GetN(); ++i) {
                Ptr<LteUeNetDevice> ueDev = ueDevs.Get(i)->GetObject<LteUeNetDevice>();
                if (ueDev->GetImsi() != imsi) continue;

                Ptr<Node> ueNode = ueDev->GetNode();
                double sinr = CalculateSinr(ueNode, sbsNode, txPowerDbm, pathLossExponent, *globalEnbDevs, false);

                if (!std::isnan(sinr) && !std::isinf(sinr)) {
                    sum += sinr;
                    count++;
                }
                break;
            }
        }

        double avg = (count > 0) ? sum / count : 0.0;
        sbsSinrAverage[sbsId] = avg;

        // Log to file
        // logFile << tNow << "," << sbsId << "," << avg << "," << count << "\n";
    }
}

void SampleMacroSinr(NetDeviceContainer& ueDevs, uint32_t macroNodeId, double pathLossExponent) {
    double tNow = Simulator::Now().GetSeconds();

    Ptr<Node> macroNode = NodeList::GetNode(macroNodeId);
    if (!macroNode) return;

    Ptr<LteEnbNetDevice> enbDev = macroNode->GetDevice(0)->GetObject<LteEnbNetDevice>();
    double txPowerDbm = enbDev->GetPhy()->GetTxPower();

    double sum = 0.0;
    uint32_t count = 0;

    for (uint32_t i = 0; i < ueDevs.GetN(); ++i) {
        Ptr<LteUeNetDevice> ueDev = ueDevs.Get(i)->GetObject<LteUeNetDevice>();
        Ptr<Node> ueNode = ueDev->GetNode();

        // Only consider UEs actually connected to macro
        uint64_t imsi = ueDev->GetImsi();
        if (ueToSbsMap.count(imsi) && ueToSbsMap[imsi] == macroNodeId) {
            if (!ueActivityMap[imsi]) continue; 
            double sinr = CalculateSinr(ueNode, macroNode, txPowerDbm, pathLossExponent, *globalEnbDevs, false);
            if (!std::isnan(sinr) && !std::isinf(sinr)) {
                sum += sinr;
                count++;
            }
        }
    }

    macroSinrAverage = (count > 0) ? sum / count : 0.0;
    std::cout << "[Macro BS] Avg SINR at " << tNow << "s = " << macroSinrAverage << " dB (from " << count << " UEs)\n";
}

// void StartSinrSampling(Time duration)
// {
//     double interval = 0.01;
//     Time startTime = Simulator::Now();
//     Time endTime = startTime + duration;

//     // Schedule SINR sampling every 1ms within the duration window
//     for (Time t = Seconds(0.0); t < duration; t += Seconds(interval)) {
//         Simulator::Schedule(t, &SampleSbsSinr, std::ref(*globalEnbDevs), std::ref(*globalUeDevs), 3.5);
//     }

//     // Schedule clearing SINR samples right after the 0.2s window
//     Simulator::Schedule(duration, [](){
//         sbsSinrSamples.clear();
//     });
// }

// static std::ofstream&
// GetSinrActivationLog(const std::string& filename = "sinr_activation_log.csv",
//                      bool               append    = false)
// {
//     static std::ofstream logFile;
//     static bool          initialized = false;

//     if (!initialized)
//     {
//         std::ios_base::openmode mode = std::ios::out |
//                                        (append ? std::ios::app : std::ios::trunc);
//         logFile.open(filename, mode);
//         if (!logFile.is_open())
//         {
//             throw std::runtime_error("Cannot open SINR activation log file: " + filename);
//         }

//         // Write header
//         logFile << "time_s,label,sbsNodeId,ueImsi,distance_m,sinr_db\n";
//         initialized = true;
//     }

//     return logFile;
// }

// void LogSinrDuringTransition(uint32_t sbsNodeId, std::string label, Time endTime)
// {
//     double now = Simulator::Now().GetSeconds();
//     auto& logFile = GetSinrActivationLog();

//     Ptr<Node> sbsNode = NodeList::GetNode(sbsNodeId);
//     if (!sbsNode) {
//         // std::cout << now << "s: [ERROR] SBS NodeId " << sbsNodeId << " not found!\n";
//         return;
//     }

//     Ptr<MobilityModel> sbsMob = sbsNode->GetObject<MobilityModel>();
//     if (!sbsMob) {
//         // std::cout << now << "s: [ERROR] SBS MobilityModel not found for NodeId " << sbsNodeId << "\n";
//         return;
//     }

//     Ptr<SmallCellEnergyModel> energyModel = sbsEnergyModels[sbsNodeId];
//     if (!energyModel) {
//         std::cout << now << "s: [ERROR] No energy model found for SBS NodeId " << sbsNodeId << "\n";
//         return;
//     }

//     double txPowerDbm = energyModel->GetTransmissionPower();

//     uint32_t totalUes = 0, consideredUes = 0;
//     for (uint32_t i = 0; i < NodeList::GetNNodes(); ++i)
//     {
//         Ptr<Node> ue = NodeList::GetNode(i);
//         Ptr<NetDevice> dev = ue->GetDevice(0);
//         Ptr<LteUeNetDevice> ueDev = DynamicCast<LteUeNetDevice>(dev);
//         if (!ueDev) continue;

//         totalUes++;

//         Ptr<MobilityModel> ueMob = ue->GetObject<MobilityModel>();
//         if (!ueMob) continue;

//         double dist = ueMob->GetDistanceFrom(sbsMob);

//         double sinrDb = CalculateSinr(
//             ue, sbsNode, txPowerDbm, 3.5, *globalEnbDevs, false);

//         if (std::isnan(sinrDb) || std::isinf(sinrDb)) {
//             std::cout << now << "s: [WARNING] Invalid SINR for UE IMSI " << ueDev->GetImsi()
//                       << " at distance " << dist << " m\n";
//             continue;
//         }

//         consideredUes++;
//         logFile << now << "," << label << "," << sbsNodeId << ","
//                 << ueDev->GetImsi() << "," << dist << "," << sinrDb
//                 << "," << txPowerDbm << "\n";
//         logFile.flush();
//     }

//     if (Simulator::Now() + MilliSeconds(1) < endTime) {
//         Simulator::Schedule(MilliSeconds(1), &LogSinrDuringTransition,
//                             sbsNodeId, label, endTime);
//     } else {
//         std::cout << now << "s: End of transition window reached — stopping SINR logging.\n";
//     }
// }



// === START OF RL ===

class LteGymEnv : public OpenGymEnv
{
public:
    LteGymEnv(const std::vector<Ptr<Node>>& sbsNodes,
              const std::map<uint32_t, Ptr<SmallCellEnergyModel>>& energyModels, double simTime);

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
    double m_simulationTime;

    std::map<uint32_t, Time> m_lastUpdateTimes;
    std::map<uint32_t, double> m_lastEnergyConsumptions;
    std::map<uint32_t, uint32_t> m_activeUeCounts;
    std::map<uint32_t, SmallCellEnergyModel::SmallCellState> m_lastStates;
    double m_lastAvgMacroSinr;
    std::map<uint32_t, double> m_lastSbsSinrAverage;
    double m_globalAvgSINR = 0.0;

    void CollectEnvData();
};

LteGymEnv::LteGymEnv(const std::vector<Ptr<Node>>& sbsNodes,
                     const std::map<uint32_t, Ptr<SmallCellEnergyModel>>& energyModels,
                     double simTime)
    : m_sbsNodes(sbsNodes), m_energyModels(energyModels), m_simulationTime(simTime)
{
    for (const auto& [id, model] : energyModels) {
        m_lastUpdateTimes[id] = Seconds(0);
        m_lastEnergyConsumptions[id] = 0.0;
        m_activeUeCounts[id] = 0;
        m_lastStates[id] = model->GetState(); 
    }
}

Ptr<OpenGymSpace> LteGymEnv::GetObservationSpace()
{
    std::cout << "GetObservationSpace() called" << std::endl;
    uint32_t dim = m_energyModels.size() * 3;  // 3 values per SBS
    std::vector<uint32_t> shape = {dim};
    float low = 0.0;
    float high = 100000.0;  // Adjust as needed
    return CreateObject<OpenGymBoxSpace>(low, high, shape, "float32");
}

Ptr<OpenGymSpace> LteGymEnv::GetActionSpace()
{
    uint32_t numSbs = 3;  
    std::vector<uint32_t> shape = {numSbs};  

    // Print action space for debugging
    std::cout << "Action Space: [";
    for (size_t i = 0; i < shape.size(); ++i) {
        std::cout << shape[i];
        if (i < shape.size() - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    float low = 0;   // Minimum action value (ACTIVE)
    float high = 3;  // Maximum action value (SM3)
    std::string dtype = TypeNameGet<float>(); // or "float32"

    Ptr<OpenGymBoxSpace> space = CreateObject<OpenGymBoxSpace>(low, high, shape, dtype);
    return space;
}


Ptr<OpenGymDataContainer> LteGymEnv::GetObservation()
{
    std::cout << "GetObservation() called" << std::endl;

    CollectEnvData();  // Collect fresh data every step

    std::vector<double> obs;
    for (const auto& [sbsId, model] : m_energyModels) {
        obs.push_back((double)m_activeUeCounts[sbsId]);
        obs.push_back((double)model->GetState());
        obs.push_back(m_lastEnergyConsumptions[sbsId]);  // Instantaneous power
        obs.push_back(m_globalAvgSINR);
        obs.push_back(model->IsTransitioning() ? 1.0 : 0.0);
    }

    // Add global SINR as final part of observation:
   

    Ptr<OpenGymBoxContainer<double>> container = CreateObject<OpenGymBoxContainer<double>>(std::vector<uint32_t>{(uint32_t)obs.size()});
    container->SetData(obs);
    return container;
}

void LteGymEnv::CollectEnvData()
{
    Time now = Simulator::Now();
    double simTimeHours = fmod((now.GetSeconds() / m_simulationTime) * 24.0, 24.0);

    SampleSbsSinr(*globalEnbDevs, *globalUeDevs, 3.5);
    SampleMacroSinr(*globalUeDevs, /*macroNodeId=*/3, 3.5);

    std::cout << "[LOG] SimTime: " << now.GetSeconds()
              << "s  (Hour of day: " << simTimeHours << ")" << std::endl;

    uint32_t totalActiveUEs = 0;
    double totalSINRWeightedSum = 0.0;

    for (const auto& [sbsId, model] : m_energyModels) {

        uint32_t activeUes = 0;
        std::cout << "  [SBS " << sbsId << "] UEs mapped: ";
        if (sbsToUeMap.find(sbsId) != sbsToUeMap.end()) {
            for (uint64_t imsi : sbsToUeMap[sbsId]) {
                std::cout << imsi << (ueActivityMap[imsi] ? "(active)" : "(inactive)") << " ";
                if (ueActivityMap[imsi]) {
                    activeUes++;
                }
            }
        } else {
            std::cout << "None";
        }
        std::cout << "| Active UEs: " << activeUes << std::endl;
        m_activeUeCounts[sbsId] = activeUes;

        // Energy: use instantaneous power
        double power = model->GetTotalPowerConsumption();
        m_lastEnergyConsumptions[sbsId] = power;

        // SINR for this SBS
        double sinrAvg = 0.0;
        if (sbsSinrAverage.count(sbsId)) {
            sinrAvg = sbsSinrAverage[sbsId];
        }
        std::cout << "  [SBS " << sbsId << "] Avg SINR this step: " << sinrAvg << " dB" << std::endl;

        // === Accumulate for global SINR calculation ===
        totalActiveUEs += activeUes;
        totalSINRWeightedSum += sinrAvg * activeUes;
    }

    // Compute global average SINR across all active UEs:
    if (totalActiveUEs > 0)
        m_globalAvgSINR = totalSINRWeightedSum / totalActiveUEs;
    else
        m_globalAvgSINR = 0.0;
}


float LteGymEnv::GetReward()
{
    std::cout << "GetReward() called" << std::endl;

    float totalReward = 0.0;

    for (const auto& [sbsId, model] : m_energyModels)
    {
        // === ENERGY ===
        double power = model->GetTotalPowerConsumption();  // W
        double energyScore = -1.0 * power;

        // === SINR ===
        double sinrDb = 0.0;
        if (sbsSinrAverage.count(sbsId)) {
            sinrDb = sbsSinrAverage[sbsId];
        }
        m_lastSbsSinrAverage[sbsId] = sinrDb; 
        double qosScore = 3.0 * sinrDb;

        // === SWITCHING COST ===
        SmallCellEnergyModel::SmallCellState currentState = model->GetState();
        SmallCellEnergyModel::SmallCellState previousState = m_lastStates[sbsId];
        double switchingPenalty = 0.0;

        if (currentState != previousState) {
            switchingPenalty = 0.5;  // You can adjust this weight
        }

        double sbsReward = energyScore + qosScore - switchingPenalty;

        totalReward += sbsReward;

        // === UPDATE previous state for next step ===
        m_lastStates[sbsId] = currentState;

        // Debug print (optional)
        std::cout << "  [SBS " << sbsId << "] Power=" << power
                  << "W, SINR=" << sinrDb
                  << "dB, SwitchingPenalty=" << switchingPenalty
                  << ", Reward=" << sbsReward << std::endl;
    }

    // === Macro SINR contribution ===
    std::cout << "  [Macro BS] SINR = " << macroSinrAverage << " dB" << std::endl;
    m_lastAvgMacroSinr = macroSinrAverage;
    totalReward += 1.0 * macroSinrAverage;

    // Clear SINR buffer after use
    sbsSinrAverage.clear();
    macroSinrAverage = 0.0;  

    return totalReward;
}




bool LteGymEnv::ExecuteActions(Ptr<OpenGymDataContainer> action)
{
    std::cout << "ExecuteActions() called" << std::endl;
    std::cout << action << std::endl;
    Ptr<OpenGymBoxContainer<float>> box = DynamicCast<OpenGymBoxContainer<float>>(action);

    if (!box) {
        std::cout << "  [ERROR] Box container is NULL!" << std::endl;
        return false;
    }

    std::vector<float> actions = box->GetData();
    std::cout << "  Action vector size: " << actions.size() << std::endl;
    std::cout << "  Actions: ";
    for (auto a : actions) std::cout << a << " ";
    std::cout << std::endl;

    if (actions.size() != m_energyModels.size()) {
        std::cout << "  [ERROR] Actions size mismatch! Actions: " << actions.size() << ", Expected: " << m_energyModels.size() << std::endl;
        return false;
    }

    size_t i = 0;
    for (const auto& [sbsId, model] : m_energyModels) {
        model->SetState((SmallCellEnergyModel::SmallCellState)actions[i]);
        double txPowerW = model->GetTransmissionPower();
        double txPowerDbm = (txPowerW > 0) ? txPowerW : 0;

        Ptr<Node> sbsNode = NodeList::GetNode(sbsId);
        Ptr<LteEnbNetDevice> enbDev = sbsNode->GetDevice(0)->GetObject<LteEnbNetDevice>();
        Ptr<LteEnbPhy> enbPhy = enbDev->GetPhy();
        enbPhy->SetTxPower(txPowerDbm);

        std::cout << Simulator::Now().GetSeconds()
              << "s: [SBS " << sbsId << "] State: " << StateToString(model->GetState())
              << ", Tx Power: " << txPowerDbm << " dBm" << std::endl;

        ++i;
    }
    return true;
}
 


bool LteGymEnv::GetGameOver()
{
    std::cout << "GetGameOver() called" << std::endl;
    // End episode when simulation time is up (step-based)
    double now = Simulator::Now().GetSeconds();
    return now >= m_simulationTime;
}


std::string LteGymEnv::GetExtraInfo()
{
    // === Total Energy ===
    double totalEnergy = 0.0;
    for (const auto& [sbsId, model] : m_energyModels)
    {
        totalEnergy += model->GetTotalEnergyConsumption();
    }

    // === Global UE-Weighted SINR calculation ===
    double totalSbsSinr = 0.0;
    uint32_t totalUEs = 0;

    for (const auto& [sbsId, sinr] : m_lastSbsSinrAverage)
    {
        uint32_t activeUes = m_activeUeCounts[sbsId];  // <-- reuse the latest active UE count per SBS
        totalSbsSinr += sinr * activeUes;
        totalUEs += activeUes;
    }

    double avgSbsSinr = (totalUEs > 0) ? (totalSbsSinr / totalUEs) : 0.0;

    double avgMacroSinr = m_lastAvgMacroSinr;

    // === Pack everything into ExtraInfo string ===
    std::ostringstream oss;
    oss << "total_energy=" << totalEnergy
        << ";avg_sbs_sinr=" << avgSbsSinr
        << ";macro_sinr=" << avgMacroSinr;

    return oss.str();
}



// === END RL ENV DEFINITION ===

// === STEP CALLBACK FOR RL ===

void StepCallback(Ptr<OpenGymInterface> gymInterface, double stepTime, double simulationTime)
{
    double now = Simulator::Now().GetSeconds();
    // End when simulation time reached
    if (now >= simulationTime) {
        Simulator::Stop();
        return;
    }
    // Optionally, print debug info:
    // std::cout << "NS-3 RL Step at " << now << "s" << std::endl;

    // Schedule the next step (always stepTime ahead)
    Simulator::Schedule(Seconds(stepTime), &StepCallback, gymInterface, stepTime, simulationTime);
}


void ScheduleNextStateRead(double envStepTime, Ptr<OpenGymInterface> openGym)
{
  Simulator::Schedule (Seconds(envStepTime), &ScheduleNextStateRead, envStepTime, openGym);
  openGym->NotifyCurrentState();
}





// === END OF RL ===

// === Main Simulation ===
int main(int argc, char *argv[])
{
    uint16_t numSmallCells = 3;
    uint16_t numUesPerCell = 10;
    double simulationTime = 10; // seconds (should be double for step arithmetic)
    double stepTime = 0.01; // seconds
    double pathLossExponent = 3.5; // Path loss exponent for SINR calculation, 
    double sinrHysteresisDb = 2.0; // Minimum SINR improvement needed for handover (dB)
    // double handoverCheckInterval = 0.2; // seconds
    uint32_t openGymPort = 5555;
    uint32_t simSeed = 1;
    uint32_t testArg = 0;
    CommandLine cmd;
    cmd.AddValue("openGymPort", "OpenGym port", openGymPort);
    cmd.AddValue("simSeed", "Simulation seed", simSeed);
    cmd.AddValue ("simTime", "Simulation time in seconds. Default: 10s", simulationTime);
    cmd.AddValue("stepTime", "RL Step time (seconds)", stepTime);
    cmd.Parse(argc, argv);

    std::cout << "Simulation is running.";
    std::cout << "Ns3Env parameters:" << std::endl;
    std::cout << "--simulationTime: " << simulationTime << std::endl;
    std::cout << "--openGymPort: " << openGymPort << std::endl;
    // std::cout << "--envStepTime: " << envStepTime << std::endl;
    std::cout << "--seed: " << simSeed << std::endl;
    std::cout << "--testArg: " << testArg << std::endl;
    std::cout << "--testArg: " << testArg << std::endl;

    RngSeedManager::SetSeed (1);
    RngSeedManager::SetRun (simSeed);
    

    Ptr<LteHelper> lteHelper = CreateObject<LteHelper>();
    Ptr<PointToPointEpcHelper> epcHelper = CreateObject<PointToPointEpcHelper>();
    lteHelper->SetEpcHelper(epcHelper);

    // lteHelper->SetAttribute("PathlossModel", StringValue("ns3::LogDistancePropagationLossModel"));
    // lteHelper->SetHandoverAlgorithmType("ns3::A3RsrpHandoverAlgorithm");
    // lteHelper->SetHandoverAlgorithmAttribute("Hysteresis", DoubleValue(0));
    // lteHelper->SetHandoverAlgorithmAttribute("TimeToTrigger", TimeValue(MilliSeconds(0)));

    // lteHelper->SetHandoverAlgorithmType("ns3::A2A4RsrqHandoverAlgorithm");
    // lteHelper->SetHandoverAlgorithmAttribute("ServingCellThreshold", UintegerValue(28));  // Choose 0–34; lower is more aggressive
    // lteHelper->SetHandoverAlgorithmAttribute("NeighbourCellOffset", UintegerValue(1));    // Lower is more aggressive

    // lteHelper->SetAttribute("FadingModel", StringValue("ns3::TraceFadingLossModel"));

    Config::SetDefault("ns3::LteUePhy::RsrpSinrSamplePeriod", UintegerValue(1));

    // Config::SetDefault("ns3::LogDistancePropagationLossModel::Exponent", DoubleValue(3.2));
    // Config::SetDefault("ns3::TraceFadingLossModel::TraceFilename",
    //                StringValue("src/lte/model/fading-traces/fading_trace_EVA_60kmph.fad"));
    // Config::SetDefault("ns3::TraceFadingLossModel::WindowSize", TimeValue(Seconds(0.5)));
    // Config::SetDefault("ns3::TraceFadingLossModel::RbNum", UintegerValue(100));
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
    uePositionAlloc->SetAttribute("X", StringValue("ns3::UniformRandomVariable[Min=0.0|Max=1000.0]"));
    uePositionAlloc->SetAttribute("Y", StringValue("ns3::UniformRandomVariable[Min=800.0|Max=1000.0]"));

    // Configure the UE MobilityHelper with RandomWaypointMobilityModel
    // MobilityHelper ueMobility;
    // ueMobility.SetMobilityModel("ns3::RandomWaypointMobilityModel",
    //     "Speed", StringValue("ns3::UniformRandomVariable[Min=1.0|Max=50.0]"),
    //     "Pause", StringValue("ns3::ConstantRandomVariable[Constant=0.5]"),
    //     "PositionAllocator", PointerValue(uePositionAlloc));
    // ueMobility.SetPositionAllocator(uePositionAlloc);
    // ueMobility.Install(ueNodes);

    MobilityHelper ueMobility;
    ueMobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    ueMobility.SetPositionAllocator(uePositionAlloc); // for non-manual UEs
    ueMobility.Install(ueNodes);

    // Override position of UE[0] and UE[1]
    for (uint32_t i = 0; i < ueNodes.GetN(); ++i) {
        Ptr<Node> ueNode = ueNodes.Get(i);
        Ptr<MobilityModel> mob = ueNode->GetObject<MobilityModel>();

        if (i == 0) {
            // IMSI = 1 (first UE)
            mob->SetPosition(Vector(250.0, 850.0, 0.0));
        } else if (i == 1) {
            // IMSI = 2 (second UE)
            mob->SetPosition(Vector(375.0, 800.0, 0.0));
        }

        // Log IMSI + position
        uint32_t imsi = i + 1;
        // LogUePositions(ueNode, imsi);
    }


    

    // Macro eNB
    mobility.SetPositionAllocator("ns3::GridPositionAllocator",
        "MinX", DoubleValue(500.0), "MinY", DoubleValue(0.0),
        "DeltaX", DoubleValue(20.0), "GridWidth", UintegerValue(1),
        "LayoutType", StringValue("RowFirst"));
    mobility.Install(macroEnb);

    // Small cells
    mobility.SetPositionAllocator("ns3::GridPositionAllocator",
        "MinX", DoubleValue(250.0), "MinY", DoubleValue(750.0),
        "DeltaX", DoubleValue(250.0), "GridWidth", UintegerValue(numSmallCells),
        "LayoutType", StringValue("RowFirst"));
    mobility.Install(smallCellEnbs);

    // === LTE Device Installation 
    NetDeviceContainer macroEnbLteDevs = lteHelper->InstallEnbDevice(macroEnb);
    NetDeviceContainer smallCellLteDevs = lteHelper->InstallEnbDevice(smallCellEnbs);
    NetDeviceContainer ueLteDevs = lteHelper->InstallUeDevice(ueNodes);

    NetDeviceContainer enbDevs;
    enbDevs.Add(macroEnbLteDevs);
    enbDevs.Add(smallCellLteDevs);
    globalEnbDevs = &enbDevs;

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

        scEnergyModel->SetNodeId(nodeId);
        
        // Simulator::Schedule(Seconds(1.0), &UpdateEnergy, scEnergyModel);
    }


    for (uint32_t i = 0; i < macroEnbLteDevs.GetN (); ++i)
    {
        Ptr<LteEnbNetDevice> enbDev = macroEnbLteDevs.Get (i)->GetObject<LteEnbNetDevice> ();
        Ptr<LteEnbPhy>       enbPhy = enbDev->GetPhy ();

        enbPhy->SetTxPower (30.0);   // value is in dBm
    }
    InitializeSbsTxPower ();
    
    for (uint32_t i = 0; i < macroEnb.GetN(); ++i) {
        macroNodeIds.insert(macroEnb.Get(i)->GetId());
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
        
        
        // uePhy->TraceConnectWithoutContext("ReportCurrentCellRsrpSinr", MakeCallback(&ReportTraceCallback));
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
        // LogUePositions(ueNodes.Get(i), imsi);
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



    Simulator::Schedule(
        Seconds(0.001),
        &CustomSinrHandover,
        lteHelper, ueNodes, enbDevs, ueLteDevs, pathLossExponent, sinrHysteresisDb
    );

    FlowMonitorHelper flowmonHelper;
    Ptr<FlowMonitor> monitor = flowmonHelper.InstallAll();  



    globalUeDevs = &ueLteDevs;
    // === Schedule Logging ===
    // Simulator::Schedule(Seconds(1.0), &PrintRemainingEnergy, energyModels);
    
    
    // Simulator::Schedule(Seconds(1.0), &PrintSbsPositions, smallCellEnbs);
    // Simulator::Schedule(Seconds(1.0), &PrintSbsPositions, macroEnb);
    // Simulator::Schedule(Seconds(0.5), &PrintSbsConnections);
    // Simulator::Schedule(Seconds(1.0), &PrintSbsUeMappings);
    // Simulator::Schedule(Seconds(1.0), &PrintSbsUeCounts);
    // Simulator::Schedule(Seconds(0.5), &UpdateSbsSleepStates);


    // Simulator::Schedule(Seconds(1.0), &PrintUeToSbsDistances, ueNodes, cellIdToNodeId);
    // Simulator::Schedule(Seconds(1.0), &PrintUePositions, ueNodes);


    // === OPENGYM INTERFACE ===
    std::vector<Ptr<Node>> sbsNodes;
    std::map<uint32_t, Ptr<SmallCellEnergyModel>> energyModels;

    for (uint32_t i = 0; i < smallCellEnbs.GetN(); ++i) {
        Ptr<Node> sbsNode = smallCellEnbs.Get(i);
        uint32_t nodeId = sbsNode->GetId();
        Ptr<SmallCellEnergyModel> model = sbsEnergyModels[nodeId];
        sbsNodes.push_back(sbsNode);
        energyModels[nodeId] = model;
    }

    // Create environment and interface
    Ptr<OpenGymInterface> openGym = CreateObject<OpenGymInterface> (openGymPort);
    Ptr<LteGymEnv> lteEnv = CreateObject<LteGymEnv>(sbsNodes, energyModels, simulationTime);

    openGym->SetGetActionSpaceCb(MakeCallback (&LteGymEnv::GetActionSpace, lteEnv) );
    openGym->SetGetObservationSpaceCb( MakeCallback (&LteGymEnv::GetObservationSpace, lteEnv) );
    openGym->SetGetGameOverCb( MakeCallback (&LteGymEnv::GetGameOver, lteEnv) );
    openGym->SetGetObservationCb( MakeCallback (&LteGymEnv::GetObservation, lteEnv) );
    openGym->SetGetRewardCb( MakeCallback (&LteGymEnv::GetReward, lteEnv) );
    openGym->SetGetExtraInfoCb( MakeCallback (&LteGymEnv::GetExtraInfo, lteEnv) );
    openGym->SetExecuteActionsCb( MakeCallback (&LteGymEnv::ExecuteActions, lteEnv) );
    Simulator::Schedule (Seconds(0.0), &ScheduleNextStateRead, stepTime, openGym);
    // Simulator::Schedule(Seconds(0.0), &SampleSbsSinr, std::ref(enbDevs), std::ref(ueLteDevs), pathLossExponent);
    


    // GetRsrpSinrLog();
    GetDebugLog();


    std::cout <<"Simulation start"<< std::endl;
    Simulator::Stop (Seconds (simulationTime));
    Simulator::Run ();
    std::cout << "Simulation stop";

    openGym->NotifySimulationEnd();
    for (const auto& [sbsId, model] : sbsEnergyModels)
    {
        model->FlushFinalState();
    }
    Simulator::Destroy ();

    // Ptr<LteGymEnv> gymEnv = CreateObject<LteGymEnv>(sbsNodes, energyModels, simulationTime);
    // Ptr<OpenGymInterface> gymInterface = CreateObject<OpenGymInterface>(openGymPort);
    // gymInterface->Notify(gymEnv);

    // // === RL STEP-BASED SCHEDULING ===
    // // Schedule the RL step callback (starts at t=0.0, repeats every stepTime)
    // Simulator::Schedule(Seconds(0.0), &StepCallback, gymInterface, stepTime, simulationTime);
    // Simulator::Run();
    // std::cout << "=== Simulation finished at " << Simulator::Now().GetSeconds() << "s ===" << std::endl;
    // gymInterface->NotifySimulationEnd();
    // Simulator::Destroy();

    for (const auto& [sbsId, model] : sbsEnergyModels)
    {
        model->ExportStateTimeToCsv(1);
    }
    // monitor->CheckForLostPackets();
    // monitor->SerializeToXmlFile("flowmon-results.xml", true, true); 


    return 0;
}