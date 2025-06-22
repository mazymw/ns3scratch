// Helper conversions
double dBmToWatt(double dBm) { return pow(10.0, (dBm - 30) / 10.0); }
double WattTodBm(double watt) { return 10.0 * log10(watt) + 30.0; }

// Path loss model (simplified log-distance)
double GetReceivedPowerDbm(double txPowerDbm, double distance, double freqMHz = 2600.0, double pathLossExp = 3.5)
{
    if (distance < 1.0) distance = 1.0; // Avoid log(0)
    double lambda = 3e8 / (freqMHz * 1e6);
    double fspl = 20 * log10(4 * M_PI * distance / lambda);
    double pathLoss = fspl + 10 * pathLossExp * log10(distance);
    return txPowerDbm - pathLoss;
}

// Custom handover function called periodically
void CustomSinrHandover(
    Ptr<LteHelper> lteHelper,
    NodeContainer ueNodes,
    NetDeviceContainer ueLteDevs,
    NodeContainer smallCellEnbs,
    NetDeviceContainer smallCellLteDevs,
    double interval = 0.2,
    double sinrHysteresis = 2.0)
{
    // Use global cellIdToNodeId and sbsEnergyModels

    for (uint32_t i = 0; i < ueNodes.GetN(); ++i) {
        Ptr<Node> ueNode = ueNodes.Get(i);
        Ptr<LteUeNetDevice> ueDev = ueLteDevs.Get(i)->GetObject<LteUeNetDevice>();
        Ptr<LteUeRrc> rrc = ueDev->GetRrc();
        Ptr<MobilityModel> ueMobility = ueNode->GetObject<MobilityModel>();

        // Identify serving SBS
        uint16_t servingCellId = rrc->GetCellId();
        if (cellIdToNodeId.find(servingCellId) == cellIdToNodeId.end()) continue;
        uint32_t servingNodeId = cellIdToNodeId.at(servingCellId);


        // Find servingEnbDev and its index
        Ptr<NetDevice> servingEnbDev = nullptr;
        int servingIdx = -1;
        for (uint32_t idx = 0; idx < smallCellEnbs.GetN(); ++idx) {
            if (smallCellEnbs.Get(idx)->GetId() == servingNodeId) {
                servingEnbDev = smallCellLteDevs.Get(idx);
                servingIdx = idx;
                break;
            }
        }
        if (servingEnbDev == nullptr) continue;

        Ptr<Node> servingSbsNode = NodeList::GetNode(servingNodeId);
        Ptr<MobilityModel> servingSbsMobility = servingSbsNode->GetObject<MobilityModel>();
        double servingDist = ueMobility->GetDistanceFrom(servingSbsMobility);
        Ptr<SmallCellEnergyModel> servingModel = sbsEnergyModels.at(servingNodeId);
        double servingTxDbm = servingModel->GetTransmissionPower();
        if (servingTxDbm < 0) servingTxDbm = 44.77;
        double servingRxDbm = GetReceivedPowerDbm(servingTxDbm, servingDist);

        // Calculate interference+noise for serving
        double interferenceWatt = 0.0;
        for (uint32_t j = 0; j < smallCellEnbs.GetN(); ++j) {
            Ptr<Node> otherSbsNode = smallCellEnbs.Get(j);
            if (otherSbsNode->GetId() == servingNodeId) continue;
            Ptr<MobilityModel> otherMobility = otherSbsNode->GetObject<MobilityModel>();
            double otherDist = ueMobility->GetDistanceFrom(otherMobility);
            Ptr<SmallCellEnergyModel> otherModel = sbsEnergyModels.at(otherSbsNode->GetId());
            double otherTxDbm = otherModel->GetTransmissionPower();
            if (otherTxDbm < 0) otherTxDbm = 44.77;
            double interfRxDbm = GetReceivedPowerDbm(otherTxDbm, otherDist);
            interferenceWatt += dBmToWatt(interfRxDbm);
        }
        double noiseWatt = 1.38e-23 * 290 * 20e6;
        double sinrServing = dBmToWatt(servingRxDbm) / (interferenceWatt + noiseWatt);
        double sinrServingDb = 10 * log10(sinrServing);

        // Find best SBS for handover
        double bestSinrDb = sinrServingDb;
        int bestSbsIdx = servingIdx;
        Ptr<NetDevice> bestEnbDev = servingEnbDev;
        uint32_t bestSbsId = servingNodeId;

        for (uint32_t j = 0; j < smallCellEnbs.GetN(); ++j) {
            Ptr<Node> candidateSbsNode = smallCellEnbs.Get(j);
            uint32_t candidateSbsId = candidateSbsNode->GetId();
            if (candidateSbsId == servingNodeId) continue;

            Ptr<MobilityModel> candidateMob = candidateSbsNode->GetObject<MobilityModel>();
            double dist = ueMobility->GetDistanceFrom(candidateMob);
            Ptr<SmallCellEnergyModel> candidateModel = sbsEnergyModels.at(candidateSbsId);
            double candTxDbm = candidateModel->GetTransmissionPower();
            if (candTxDbm < 0) candTxDbm = 44.77;
            double rxDbm = GetReceivedPowerDbm(candTxDbm, dist);

            // Interference: all other SBS except candidate
            double interfW = 0.0;
            for (uint32_t k = 0; k < smallCellEnbs.GetN(); ++k) {
                Ptr<Node> otherNode = smallCellEnbs.Get(k);
                if (otherNode->GetId() == candidateSbsId) continue;
                Ptr<MobilityModel> otherMob = otherNode->GetObject<MobilityModel>();
                double otherDist = ueMobility->GetDistanceFrom(otherMob);
                Ptr<SmallCellEnergyModel> otherModel = sbsEnergyModels.at(otherNode->GetId());
                double otherTxDbm = otherModel->GetTransmissionPower();
                if (otherTxDbm < 0) otherTxDbm = 44.77;
                double interfDbm = GetReceivedPowerDbm(otherTxDbm, otherDist);
                interfW += dBmToWatt(interfDbm);
            }
            double candSinr = dBmToWatt(rxDbm) / (interfW + noiseWatt);
            double candSinrDb = 10 * log10(candSinr);

            if (candSinrDb > bestSinrDb) {
                bestSinrDb = candSinrDb;
                bestSbsIdx = j;
                bestEnbDev = smallCellLteDevs.Get(j);
                bestSbsId = candidateSbsId;
            }
        }

        // Handover if best is not serving and SINR gain is enough
        if (bestSbsId != servingNodeId && (bestSinrDb - sinrServingDb) > sinrHysteresis) {
            std::cout << Simulator::Now().GetSeconds() << "s: [CUSTOM HO] UE " << ueDev->GetImsi()
                      << " SINR " << sinrServingDb << " -> " << bestSinrDb << " dB. Handover to SBS " << bestSbsId << std::endl;
            // 4-parameter HandoverRequest: time, ueDev, sourceEnbDev, targetEnbDev
            lteHelper->HandoverRequest(Seconds(0), ueLteDevs.Get(i), servingEnbDev, bestEnbDev);
        }
    }

    // Reschedule
    Simulator::Schedule(Seconds(interval), &CustomSinrHandover,
            lteHelper, ueNodes, ueLteDevs, smallCellEnbs, smallCellLteDevs, interval, sinrHysteresis);
}

bool IsSetupComplete(
    NodeContainer ueNodes,
    NetDeviceContainer ueLteDevs)
{
    // Use global cellIdToNodeId and sbsEnergyModels
    for (uint32_t i = 0; i < ueNodes.GetN(); ++i) {
        Ptr<LteUeNetDevice> ueDev = ueLteDevs.Get(i)->GetObject<LteUeNetDevice>();
        Ptr<LteUeRrc> rrc = ueDev->GetRrc();
        uint16_t servingCellId = rrc->GetCellId();
        if (cellIdToNodeId.find(servingCellId) == cellIdToNodeId.end())
            return false;
    }
    for (const auto& entry : cellIdToNodeId) {
        uint32_t nodeId = entry.second;
        if (sbsEnergyModels.find(nodeId) == sbsEnergyModels.end())
            return false;
    }
    return true;
}

void WaitForSetupAndRunHandover(
    Ptr<LteHelper> lteHelper,
    NodeContainer ueNodes,
    NetDeviceContainer ueLteDevs,
    NodeContainer smallCellEnbs,
    NetDeviceContainer smallCellLteDevs,
    double interval = 0.2,
    double sinrHysteresis = 2.0)
{
    if (IsSetupComplete(ueNodes, ueLteDevs)) {
        std::cout << Simulator::Now().GetSeconds()
                  << "s: [INFO] Setup complete, starting custom SINR-based handover." << std::endl;
        Simulator::Schedule(Seconds(0.01), &CustomSinrHandover,
            lteHelper, ueNodes, ueLteDevs, smallCellEnbs, smallCellLteDevs, interval, sinrHysteresis);
    } else {
        Simulator::Schedule(Seconds(0.05), &WaitForSetupAndRunHandover,
            lteHelper, ueNodes, ueLteDevs, smallCellEnbs, smallCellLteDevs, interval, sinrHysteresis);
    }
}

    Simulator::Schedule(Seconds(0.1), &WaitForSetupAndRunHandover,
        lteHelper, ueNodes, ueLteDevs, smallCellEnbs, smallCellLteDevs, 0.2, 2.0);