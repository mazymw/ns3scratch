
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
// #include "smallCellEnergyModelOriginal.h"
#include <map>
#include <vector>
#include <unordered_map>
#include <cmath>

// #include "ns3/config-store.h"
// #include "ns3/string.h"
// #include "ns3/config.h"

using namespace ns3;
using namespace ns3::energy;

// Global tracking maps
std::map<uint32_t, uint32_t> sbsToUeCount;         // SBS NodeId -> UE count
std::map<uint64_t, uint32_t> ueToSbsMap;           // IMSI -> SBS NodeId
std::map<uint16_t, uint32_t> cellIdToNodeId;       // CellId -> SBS NodeId
std::unordered_map<uint32_t, Ptr<SmallCellEnergyModel>> sbsEnergyModels;
std::map<std::pair<uint16_t, uint16_t>, uint64_t> rntiToImsiMap;

// Global map to store UE IMSI and the corresponding SBS NodeId
std::map<uint32_t, std::vector<uint64_t>> sbsToUeMap;

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
    uint32_t sleepThreshold = 2;

    for (const auto& entry : sbsToUeMap)
    {
        uint32_t sbsNodeId = entry.first;
        size_t numConnectedUes = entry.second.size();

        if (sbsEnergyModels.find(sbsNodeId) != sbsEnergyModels.end())
        {
            Ptr<SmallCellEnergyModel> scEnergyModel = sbsEnergyModels[sbsNodeId];

            if (numConnectedUes < sleepThreshold)
            {
                scEnergyModel->SetState(SmallCellEnergyModel::SM2);  // Or SM3
            }
            else
            {
                scEnergyModel->SetState(SmallCellEnergyModel::ACTIVE);
            }

            auto state = scEnergyModel->GetState();
            double power = scEnergyModel->GetPowerConsumption();

            Ptr<Node> sbsNode = NodeList::GetNode(sbsNodeId);

            // Get the LteEnbNetDevice (assumes only one device)
            Ptr<LteEnbNetDevice> enbDev = sbsNode->GetDevice(0)->GetObject<LteEnbNetDevice>();
            Ptr<LteEnbPhy> enbPhy = enbDev->GetPhy();
            double txPowerDbm = power;
            enbPhy->SetTxPower(txPowerDbm);
 
            std::cout << Simulator::Now().GetSeconds()
                      << "s: SBS " << sbsNodeId
                      << " is in state " << StateToString(state)
                      << ", consuming power = " << txPowerDbm << " dbm" << std::endl;
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



// === Main Simulation ===
int main(int argc, char *argv[])
{
    uint16_t numSmallCells = 3;
    uint16_t numUesPerCell = 3;
    uint16_t simulationTime = 15;
    CommandLine cmd;
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

    Config::SetDefault("ns3::LteUePhy::RsrpSinrSamplePeriod", UintegerValue(100));

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

        // === Get Positions of SBSs ===
    // std::vector<Vector> sbsPositions;
    // for (uint32_t i = 0; i < smallCellEnbs.GetN(); ++i) {
    //     Ptr<Node> enb = smallCellEnbs.Get(i);
    //     Ptr<MobilityModel> mob = enb->GetObject<MobilityModel>();
    //     if (mob) {
    //         sbsPositions.push_back(mob->GetPosition());
    //     }
    // }

    // Ptr<UniformRandomVariable> randX = CreateObject<UniformRandomVariable>();
    // randX->SetAttribute("Min", DoubleValue(0.0));
    // randX->SetAttribute("Max", DoubleValue(600.0));

    // Ptr<UniformRandomVariable> randY = CreateObject<UniformRandomVariable>();
    // randY->SetAttribute("Min", DoubleValue(0.0));
    // randY->SetAttribute("Max", DoubleValue(150.0));

    // // === Set UE Positions Ensuring Minimum Distance from SBSs ===
    // double minDistance = 10.0; // Minimum distance in meters

    // for (uint32_t i = 0; i < ueNodes.GetN(); ++i) {
    //     Ptr<Node> ue = ueNodes.Get(i);
    //     Vector uePos;
    //     bool valid = false;

    //     while (!valid) {
    //         double x = randX->GetValue();
    //         double y = randY->GetValue();
    //         uePos = Vector(x, y, 0.0);

    //         valid = true;
    //         for (const Vector& sbsPos : sbsPositions) {
    //             if (CalculateDistance(uePos, sbsPos) < minDistance) {
    //                 valid = false;
    //                 break;
    //             }
    //         }
    //     }

    //     Ptr<MobilityModel> mob = CreateObject<ConstantPositionMobilityModel>();
    //     mob->SetPosition(uePos);
    //     ue->AggregateObject(mob);

    //     // Optional: Log UE position with IMSI = i+1
    //     LogUePositions(ue, i + 1);
    // }
    

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
    // Config::SetDefault("ns3::LteUePhy::RsrpSinrSamplePeriod", TimeValue(MilliSeconds(1)));
    // Config::Connect("/NodeList/4/DeviceList/4/LteUePhy/ReportCurrentCellRsrpSinr", MakeCallback(&SinrTraceCallback));
    // for (uint16_t i = 0; i < ueNodes.GetN(); ++i)
    // {
    //     uint32_t nodeId = ueNodes.Get(i)->GetId();

    //     std::ostringstream path;
    //     path << "/NodeList/" << nodeId
    //         << "/DeviceList/0/LteUeNetDevice/Phy/ReportCurrentCellRsrpSinr";

    //     std::cout << "Scheduling connect to: " << path.str() << std::endl;

    //     Simulator::Schedule(Seconds(1.0), &Config::Connect,
    //         path.str(), MakeCallback(&SinrTraceCallback));
    // }
    
    // for (uint16_t i = 0; i < ueNodes.GetN(); ++i)
    // {
    //     uint32_t nodeId = ueNodes.Get(i)->GetId();

    //     std::ostringstream path;
    //     path << "/NodeList/" << nodeId
    //         << "/DeviceList/0/LteUeNetDevice/ComponentCarrierMapUe/0/ComponentCarrierUe/LteUePhy/ReportCurrentCellRsrpSinr";

    //     std::cout << "Trying to connect to: " << path.str() << std::endl;

    //     Config::Connect(path.str(), MakeCallback(&SinrTraceCallback));
    // }

    // === Attach callback for handover + connection tracking ===
    for (uint32_t i = 0; i < ueLteDevs.GetN(); ++i)
    {
        Ptr<NetDevice> dev = ueLteDevs.Get(i);
        Ptr<LteUeNetDevice> ueLte = dev->GetObject<LteUeNetDevice>();
        Ptr<LteUeRrc> rrc = ueLte->GetRrc();
        rrc->TraceConnect("ConnectionEstablished", "", MakeCallback(&NotifyConnectionEstablishedUe));
        rrc->TraceConnect("HandoverEndOk", "", MakeCallback(&NotifyHandoverEndOkUe));

    }

    // // === Energy model attachment ===
    // std::vector<Ptr<SmallCellEnergyModelOriginal>> energyModels(llCells);
    // for (uint16_t i = 0; i < numSmallCells; ++i)
    // {
    //     energyModels[i] = CreateObject<SmallCellEnergyModelOriginal>();
    //     energyModels[i]->SetState("active");
    //     energyModels[i]->SetInitialEnergy(100.0);
    //     smallCellEnbs.Get(i)->AggregateObject(energyModels[i]);

    //     Simulator::Schedule(Seconds(i + 1.0), &SmallCellEnergyModelOriginal::UpdateEnergyConsumption,
    //                         energyModels[i], Seconds(1.0));
    // }

    // Config::Connect("/NodeList/*/DeviceList/*/LteEnbRrc/HandoverStart",
    //     MakeCallback(&HandoverStartCallback));

    // Config::Connect("/NodeList/*/DeviceList/*/LteEnbRrc/HandoverEndOk",
    //         MakeCallback(&HandoverEndOkCallback));



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

    Simulator::Stop(Seconds(simulationTime));
    Simulator::Run();
    Simulator::Destroy();

    return 0;
}