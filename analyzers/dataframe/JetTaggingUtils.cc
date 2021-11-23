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

/* //get parton indices before this
ROOT::VecOps::RVec<int> JetTaggingUtils::get_flavour_gm(ROOT::VecOps::RVec<fastjet::PseudoJet> in, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin){
  ROOT::VecOps::RVec<int> result(in.size(),0);

  int loopcount =0;
  for (size_t i = 0; i < MCin.size(); ++i) {
    auto & parton = MCin[i];
    //Select partons only (for pythia 71-79):
    //if (parton.generatorStatus>80 || parton.generatorStatus<70) continue;
    //Select outgoing particles only - use ONLY for Zqq events
    if (parton.generatorStatus!=23) continue;
    if (parton.PDG > 5) continue;
    
    for (size_t j = 0; j < in.size(); ++j) {
      auto & p = in[j];

      for (size_t k = 0; k < p.constituents.size(); ++k) {
	
      }

      if (angle <= 0.3) result[j] = std::max(result[j], std::abs ( parton.PDG ));
    }
  }

  return result;
}
*/
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
