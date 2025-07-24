// Function to print SBS -> UE IMSI mapping
// void PrintSbsUeMappings()
// {
//     NS_LOG_UNCOND("=== SBS -> UE IMSI Mappings at " << Simulator::Now().GetSeconds() << "s ===");

//     for (const auto& entry : sbsToUeMap)
//     {
//         uint32_t sbsNodeId = entry.first;
//         const std::vector<uint64_t>& ueList = entry.second;

//         NS_LOG_UNCOND("SBS NodeId " << sbsNodeId << " has the following UE(s) connected:");
//         for (uint64_t imsi : ueList)
//         {
//             NS_LOG_UNCOND("  - UE IMSI: " << imsi);
//         }
//     }

//     Simulator::Schedule(Seconds(1.0), &PrintSbsUeMappings);
// }

// void PrintSbsConnections()
// {
//     std::cout << "=== SBS Connection Summary at " << Simulator::Now().GetSeconds() << "s ===" << std::endl;
//     for (const auto& pair : sbsToUeMap)
//     {
//         uint32_t sbsNodeId = pair.first;
//         const std::vector<uint64_t>& ueImsis = pair.second;

//         std::cout << "SBS NodeId " << sbsNodeId << " has the following UE(s) connected:" << std::endl;
//         for (uint64_t imsi : ueImsis)
//         {
//             std::cout << "  - UE IMSI: " << imsi << std::endl;
//         }
//     }
//     Simulator::Schedule(Seconds(0.5), &PrintSbsConnections);
// }

// void PrintUePositions(NodeContainer ueNodes)
// {
//     std::cout << "=== UE Positions at " << Simulator::Now().GetSeconds() << "s ===" << std::endl;
//     for (uint32_t i = 0; i < ueNodes.GetN(); ++i)
//     {
//         Ptr<Node> ue = ueNodes.Get(i);
//         Ptr<MobilityModel> mobility = ue->GetObject<MobilityModel>();
//         if (mobility)
//         {
//             Vector pos = mobility->GetPosition();
//             std::cout << "  UE[" << i + 1 << "] Position: (" << pos.x << ", " << pos.y << ", " << pos.z << ")" << std::endl;
//         }
//     }

//     Simulator::Schedule(Seconds(1.0), &PrintUePositions, ueNodes);
// }

// static std::ofstream&
// GetRsrpSinrLog(const std::string& filename = "rsrp-sinr-log.csv",
//                bool               append    = false)
// {
//     static std::ofstream logFile;
//     static bool          init = false;

//     if (!init)
//     {
//         std::ios_base::openmode mode = std::ios::out |
//                                        (append ? std::ios::app
//                                                : std::ios::trunc);

//         logFile.open(filename, mode);
//         if (!logFile.is_open())
//         {
//             throw std::runtime_error("Cannot open log file: " + filename);
//         }

//         if (!append)      // header once
//             logFile << "time_s,nodeId,rnti,imsi,rsrp_w,sinr_lin,ccId\n";

//         init = true;
//     }
//     return logFile;
// }

// void PrintSbsUeCounts()
// {
//     std::cout << "=== SBS UE Count Summary at " << Simulator::Now().GetSeconds() << "s ===" << std::endl;

//     for (const auto& entry : sbsToUeMap)
//     {
//         uint32_t sbsNodeId = entry.first;
//         size_t ueCount = entry.second.size();

//         std::cout << "SBS NodeId " << sbsNodeId << " has " << ueCount << " UE(s) connected." << std::endl;
//     }

//     // Schedule this function to run again in 1 second
//     Simulator::Schedule(Seconds(1.0), &PrintSbsUeCounts);
// }


// void UpdateSbsSleepStates()
// {
//     // uint32_t sleepThreshold = 1; // Minimum active UEs to stay awake

//     for (const auto& entry : sbsToUeMap)
//     {
//         uint32_t sbsNodeId = entry.first;
//         const auto& ueList = entry.second;
        
//         if (sbsEnergyModels.find(sbsNodeId) != sbsEnergyModels.end())
//         {
//             Ptr<SmallCellEnergyModel> scEnergyModel = sbsEnergyModels[sbsNodeId];
            
//             // Count ACTIVE UEs connected to this SBS
//             uint32_t activeUeCount = 0;
//             for (uint64_t imsi : ueList) {
//                 if (ueActivityMap[imsi]) {
//                     activeUeCount++;
//                 }
//             }

//             // Set state based on active UEs
//             // if (activeUeCount == 0) {
//             //     scEnergyModel->SetState(SmallCellEnergyModel::SM3); // Deep sleep
//             // } 
//             // else if (activeUeCount <= sleepThreshold) {
//             //     scEnergyModel->SetState(SmallCellEnergyModel::SM1); // Light sleep
//             // }
//             // else {
//                 scEnergyModel->SetState(SmallCellEnergyModel::ACTIVE);
//             // }

//             // Update SBS transmission power
//             Ptr<Node> sbsNode = NodeList::GetNode(sbsNodeId);
//             Ptr<LteEnbNetDevice> enbDev = sbsNode->GetDevice(0)->GetObject<LteEnbNetDevice>();
//             Ptr<LteEnbPhy> enbPhy = enbDev->GetPhy();
            
//             // Get power from energy model (convert W to dBm if needed)
//             double txPowerDbm = scEnergyModel->GetTransmissionPower();
//             // double txPowerDbm = 10 * log10(txPowerW) + 30; // Convert W to dBm
//             enbPhy->SetTxPower(txPowerDbm);

//             std::cout << Simulator::Now().GetSeconds()
//                       << "s: SBS " << sbsNodeId << " has " << activeUeCount << " active UEs"
//                       << ", State: " << StateToString(scEnergyModel->GetState())
//                       << ", Tx Power: " << txPowerDbm << " dBm" << std::endl;
//         }
//     }

//     Simulator::Schedule(Seconds(0.5), &UpdateSbsSleepStates);
// }


// void SinrTraceCallback(std::string context, double rsrp, double sinr)
// {
//     // Optional: parse IMSI from context string if needed
//     std::cout << "RSRP: " << rsrp << " dBm, SINR: " << sinr << " dB" << std::endl;
// }


// void ReportTraceCallback(uint16_t cellId,
//                     uint16_t rnti,
//                     double   rsrp,   // linear W
//                     double   sinr,   // linear ratio
//                     uint8_t  componentCarrierId)
// {
//     uint32_t nodeId = cellIdToNodeId[cellId];
//     double   t      = Simulator::Now().GetSeconds();

//     auto key     = std::make_pair(cellId, rnti);
//     auto it      = rntiToImsiMap.find(key);
//     uint64_t imsi = (it != rntiToImsiMap.end()) ? it->second : 0;

//     GetRsrpSinrLog()               // <-- log entry
//         << t        << ','
//         << nodeId   << ','
//         << rnti     << ','
//         << imsi     << ','
//         << rsrp     << ','
//         << sinr     << ','
//         << uint32_t(componentCarrierId)
//         << '\n';
// }



// void PrintUeToSbsDistances(NodeContainer ueNodes,
//                            std::map<uint16_t, uint32_t> cellIdToNodeId)
// {
//     for (uint32_t i = 0; i < ueNodes.GetN(); ++i)
//     {
//         Ptr<Node> ue = ueNodes.Get(i);
//         Ptr<MobilityModel> ueMobility = ue->GetObject<MobilityModel>();
//         Vector uePos = ueMobility->GetPosition();

//         Ptr<NetDevice> ueDevice = ue->GetDevice(0);
//         Ptr<LteUeNetDevice> lteUeDev = DynamicCast<LteUeNetDevice>(ueDevice);

//         uint64_t imsi = lteUeDev->GetImsi(); 
//         uint16_t cellId = lteUeDev->GetRrc()->GetCellId();

//         if (cellIdToNodeId.find(cellId) != cellIdToNodeId.end())
//         {
//             uint32_t sbsNodeId = cellIdToNodeId[cellId];
//             Ptr<Node> sbsNode = NodeList::GetNode(sbsNodeId);
//             Ptr<MobilityModel> sbsMobility = sbsNode->GetObject<MobilityModel>();
//             Vector sbsPos = sbsMobility->GetPosition();

//             double distance = std::sqrt(std::pow(uePos.x - sbsPos.x, 2) +
//                                         std::pow(uePos.y - sbsPos.y, 2) +
//                                         std::pow(uePos.z - sbsPos.z, 2));

//             std::cout << "UE IMSI " << imsi
//                       << " connected to SBS NodeId " << sbsNodeId
//                       << " (CellId " << cellId << ")"
//                       << " at distance: " << distance << " m" << std::endl;
//         }
//         else
//         {
//             std::cout << "UE IMSI " << imsi
//                       << " is connected to unknown CellId " << cellId << std::endl;
//         }
//     }
// }

// std::vector<Vector> GetSbsPositions(NodeContainer smallCellEnbs)
// {
//     std::vector<Vector> sbsPositions;

//     for (uint32_t i = 0; i < smallCellEnbs.GetN(); ++i)
//     {
//         Ptr<Node> sbs = smallCellEnbs.Get(i);
//         Ptr<MobilityModel> sbsMobility = sbs->GetObject<MobilityModel>();

//         if (sbsMobility)
//         {
//             Vector pos = sbsMobility->GetPosition();
//             sbsPositions.push_back(pos);
//         }
//     }

//     return sbsPositions;
// }

// uint32_t CountActiveUesNearSbs(Ptr<Node> sbsNode, NodeContainer ueNodes, double radius)
// {
//     Ptr<MobilityModel> sbsMobility = sbsNode->GetObject<MobilityModel>();
//     uint32_t count = 0;

//     for (uint32_t i = 0; i < ueNodes.GetN(); ++i)
//     {
//         Ptr<Node> ueNode = ueNodes.Get(i);
//         Ptr<MobilityModel> ueMobility = ueNode->GetObject<MobilityModel>();

//         Ptr<NetDevice> ueDev = ueNode->GetDevice(0);
//         Ptr<LteUeNetDevice> ueLteDev = ueDev->GetObject<LteUeNetDevice>();
//         uint64_t imsi = ueLteDev->GetImsi();

//         if (ueActivityMap[imsi]) // active
//         {
//             double dist = ueMobility->GetDistanceFrom(sbsMobility);
//             if (dist <= radius)
//                 count++;
//         }
//     }

//     return count;
// }

// Add this helper function to get active UE count per SBS
// uint32_t GetActiveUeCount(uint32_t sbsNodeId)
// {
//     uint32_t count = 0;
//     if (sbsToUeMap.find(sbsNodeId) != sbsToUeMap.end()) {
//         for (uint64_t imsi : sbsToUeMap[sbsNodeId]) {
//             if (ueActivityMap[imsi]) {
//                 count++;
//             }
//         }
//     }
//     return count;
// }

// Add this to print active UE counts
// void PrintActiveUeCounts()
// {
//     std::cout << "=== Active UE Counts at " << Simulator::Now().GetSeconds() << "s ===" << std::endl;
//     for (const auto& entry : sbsToUeMap) {
//         uint32_t sbsNodeId = entry.first;
//         std::cout << "SBS " << sbsNodeId << ": " 
//                   << GetActiveUeCount(sbsNodeId) << " active UEs" << std::endl;
//     }
//     Simulator::Schedule(Seconds(1.0), &PrintActiveUeCounts);
// }


// static std::ofstream&
// GetSbsSinrLog(const std::string& filename = "sbs_avg_sinr.csv", bool append = false)
// {
//     static std::ofstream logFile;
//     static bool initialized = false;

//     if (!initialized) {
//         std::ios_base::openmode mode = std::ios::out | (append ? std::ios::app : std::ios::trunc);
//         logFile.open(filename, mode);
//         if (!logFile.is_open()) {
//             throw std::runtime_error("Cannot open SINR log file: " + filename);
//         }

//         // Write CSV header
//         logFile << "time_s,sbs_id,avg_sinr_db,valid_ue_count\n";
//         initialized = true;
//     }

//     return logFile;
// }
