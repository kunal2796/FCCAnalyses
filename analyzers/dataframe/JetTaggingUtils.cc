#include "JetTaggingUtils.h"
using namespace JetTaggingUtils;

ROOT::VecOps::RVec<int> JetTaggingUtils::get_flavour(ROOT::VecOps::RVec<fastjet::PseudoJet> in, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin){
  ROOT::VecOps::RVec<int> result(in.size(),0);

  int loopcount =0;
  for (size_t i = 0; i < MCin.size(); ++i) {
    auto & parton = MCin[i];
    //Select partons only (for pythia 71-79):
    //if (parton.generatorStatus>80 || parton.generatorStatus<70) continue;
    //Select outgoing partons from the hardest process - use ONLY(?) for Zqq events
    //(For Zqq, ideally, don't even need anything else, just identify the partons)
    if (parton.generatorStatus!=23) continue;
    if (parton.PDG > 5) continue;
    ROOT::Math::PxPyPzMVector lv(parton.momentum.x, parton.momentum.y, parton.momentum.z, parton.mass);

    for (size_t j = 0; j < in.size(); ++j) {
      auto & p = in[j];
      //float dEta = lv.Eta() - p.eta();
      //float dPhi = lv.Phi() - p.phi();
      //float deltaR = sqrt(dEta*dEta+dPhi*dPhi);
      //if (deltaR <= 0.5 && gRandom->Uniform() <= efficiency) result[j] = true;
      
      Float_t dot = p.px()*parton.momentum.x+p.py()*parton.momentum.y+p.pz()*parton.momentum.z;
      Float_t lenSq1 = p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz();
      Float_t lenSq2 = parton.momentum.x*parton.momentum.x+parton.momentum.y*parton.momentum.y+parton.momentum.z*parton.momentum.z;
      Float_t norm = sqrt(lenSq1*lenSq2);
      Float_t angle = acos(dot/norm);
      if (angle <= 0.3) result[j] = std::max(result[j], std::abs ( parton.PDG ));
    }
  }

  return result;
}

ROOT::VecOps::RVec<int> JetTaggingUtils::get_flavour_qqbar(ROOT::VecOps::RVec<fastjet::PseudoJet> in, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin){
  ROOT::VecOps::RVec<int> result(in.size(),0);

  // ONLY FOR DIJETS FOR NOW
  
  int loopcount =0;
  for (size_t i = 0; i < MCin.size(); ++i) {
    auto & parton = MCin[i];
    //Select partons only (for pythia 71-79):
    //if (parton.generatorStatus>80 || parton.generatorStatus<70) continue;
    //Select outgoing partons from the hardest process - use ONLY(?) for Zqq events
    //(For Zqq, ideally, don't even need anything else, just identify the partons)
    if (parton.generatorStatus!=23) continue;
    if (parton.PDG > 5) continue;
    ROOT::Math::PxPyPzMVector lv(parton.momentum.x, parton.momentum.y, parton.momentum.z, parton.mass);

    for (size_t j = 0; j < in.size(); ++j) {
      auto & p = in[j];
      //float dEta = lv.Eta() - p.eta();
      //float dPhi = lv.Phi() - p.phi();
      //float deltaR = sqrt(dEta*dEta+dPhi*dPhi);
      //if (deltaR <= 0.5 && gRandom->Uniform() <= efficiency) result[j] = true;
      
      Float_t dot = p.px()*parton.momentum.x+p.py()*parton.momentum.y+p.pz()*parton.momentum.z;
      Float_t lenSq1 = p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz();
      Float_t lenSq2 = parton.momentum.x*parton.momentum.x+parton.momentum.y*parton.momentum.y+parton.momentum.z*parton.momentum.z;
      Float_t norm = sqrt(lenSq1*lenSq2);
      Float_t angle = acos(dot/norm);
      if (angle <= 0.3) result[j] = parton.PDG; // unlikely that both q & qbar would be within the limit
    }
  }

  return result;
}

ROOT::VecOps::RVec<int> JetTaggingUtils::get_flavour_gm(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<float> pdg_gm){
  // can later change the argument from the jet vectors to the number of jets
  // (get_njets)

  ROOT::VecOps::RVec<int> result(in.size(),0);
  if(in.size() == 0) return result;

  for (size_t i = 0; i < in.size(); i++) {
    auto & p = in[i];

    for (int ele : inJC.at(i)) {
      if (pdg_gm.at(ele) == 0) continue;
      if (abs(pdg_gm.at(ele)) > abs(result[i])) result[i] = pdg_gm.at(ele);
    }
  }
  return result;
}

ROOT::VecOps::RVec<int> JetTaggingUtils::get_flavour_gm_auto(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin, std::vector<fastjet::PseudoJet> PJin){

  // pseudoJet vector is the one before addign ghosts
  // get the pdg_gm vector here (instead of in a separate fn)
  // if pseudoJet doesn't work, replace it with (not as general) reco particles

  // push back zeros for all the particles other except for ghosts
  ROOT::VecOps::RVec<float> pdg_gm(PJin.size(),0);
  // push back the MC pdg ID for all the ghosts
  for (size_t j = 0; j < MCin.size(); ++j) {
    auto & parton = MCin[i];
    // CAUTION: use the SAME selection here as in addGhosts_pseudoJets
    //if (parton.generatorStatus!=23) continue;
    if (parton.generatorStatus>80 || parton.generatorStatus<70) continue;
    if (parton.PDG > 5) continue;   
    //if (parton.PDG > 5 && parton.PDG != 23) continue; // while adding gluons  
    pdg_gm.push_back(parton.PDG);
  }
  
  ROOT::VecOps::RVec<int> result(in.size(),0);
  if(in.size() == 0) return result;

  for (size_t i = 0; i < in.size(); i++) {
    auto & p = in[i];

    for (int ele : inJC.at(i)) {
      if (pdg_gm.at(ele) == 0) continue;
      if (abs(pdg_gm.at(ele)) > abs(result[i])) result[i] = pdg_gm.at(ele);
    }
  }
  return result;
}

ROOT::VecOps::RVec<int> JetTaggingUtils::get_btag(ROOT::VecOps::RVec<int> in, float efficiency) {
  ROOT::VecOps::RVec<int> result(in.size(),0);

  for (size_t j = 0; j < in.size(); ++j) {
    if (in.at(j)==5 && gRandom->Uniform() <= efficiency) result[j] = 1;
  }
  return result;
}

ROOT::VecOps::RVec<int>	JetTaggingUtils::get_ctag(ROOT::VecOps::RVec<int> in, float efficiency) {
  ROOT::VecOps::RVec<int> result(in.size(),0);

  for (size_t j = 0; j < in.size(); ++j) {
    if (in.at(j)==4 && gRandom->Uniform() <= efficiency) result[j] = 1;
  }
  return result;
}
