#include "JetTaggingUtils.h"
using namespace JetTaggingUtils;

ROOT::VecOps::RVec<int>
JetTaggingUtils::get_flavour(ROOT::VecOps::RVec<fastjet::PseudoJet> in,
                             ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin)
{
  ROOT::VecOps::RVec<int> result(in.size(),0);

  int loopcount =0;
  for (size_t i = 0; i < MCin.size(); ++i) {
    auto & parton = MCin[i];
    //Select partons only (for pythia8 71-79, for pythia6 2):
    if ((parton.generatorStatus>80 ||
         parton.generatorStatus<70) &&
        parton.generatorStatus != 2 ) continue;
    if (std::abs(parton.PDG) > 5 && parton.PDG!=21) continue;
    ROOT::Math::PxPyPzMVector lv(parton.momentum.x, parton.momentum.y,
                                 parton.momentum.z, parton.mass);

    for (size_t j = 0; j < in.size(); ++j) {
      auto & p = in[j];
      //float dEta = lv.Eta() - p.eta();
      //float dPhi = lv.Phi() - p.phi();
      //float deltaR = sqrt(dEta*dEta+dPhi*dPhi);
      //if (deltaR <= 0.5 && gRandom->Uniform() <= efficiency) result[j] = true;
      
      Float_t dot = p.px() * parton.momentum.x
                  + p.py() * parton.momentum.y
                  + p.pz() * parton.momentum.z;
      Float_t lenSq1 = p.px() * p.px()
                     + p.py() * p.py()
                     + p.pz() * p.pz();
      Float_t lenSq2 = parton.momentum.x * parton.momentum.x
                     + parton.momentum.y * parton.momentum.y
                     + parton.momentum.z * parton.momentum.z;
      Float_t norm = sqrt(lenSq1*lenSq2);
      Float_t angle = acos(dot/norm);

      if (angle <= 0.3) {
        if (result[j]==21 or result[j]==0) {
          // if no match before, or matched to gluon, match to
          // this particle (favour quarks over gluons)
          result[j] = std::abs ( parton.PDG );
        }
        else if (parton.PDG!=21) {
          // if matched to quark, and this is a quark, favour
          // heavier flavours
          result[j] = std::max(result[j], std::abs ( parton.PDG ));
        } else {
          // if matched to quark, and this is a gluon, keep
          // previous result (favour quark)
           ;
        }   
       }


    }
  }

  return result;
}

ROOT::VecOps::RVec<int>
JetTaggingUtils::get_btag(ROOT::VecOps::RVec<int> in,
                          float efficiency, float mistag_c,
                          float mistag_l, float mistag_g) {

  ROOT::VecOps::RVec<int> result(in.size(),0);

  for (size_t j = 0; j < in.size(); ++j) {
    if (in.at(j) ==  5 && gRandom->Uniform() <= efficiency) result[j] = 1;
    if (in.at(j) ==  4 && gRandom->Uniform() <= mistag_c) result[j] = 1;
    if (in.at(j)  <  4 && gRandom->Uniform() <= mistag_l) result[j] = 1;
    if (in.at(j) == 21 && gRandom->Uniform() <= mistag_g) result[j] = 1;
  }
  return result;
}

ROOT::VecOps::RVec<int>
JetTaggingUtils::get_ctag(ROOT::VecOps::RVec<int> in,
                          float efficiency, float mistag_b,
                          float mistag_l, float mistag_g) {

  ROOT::VecOps::RVec<int> result(in.size(),0);

  for (size_t j = 0; j < in.size(); ++j) {
    if (in.at(j) ==  4 && gRandom->Uniform() <= efficiency) result[j] = 1;
    if (in.at(j) ==  5 && gRandom->Uniform() <= mistag_b) result[j] = 1;
    if (in.at(j)  <  4 && gRandom->Uniform() <= mistag_l) result[j] = 1;
    if (in.at(j) == 21 && gRandom->Uniform() <= mistag_g) result[j] = 1;
  }
  return result;
}

ROOT::VecOps::RVec<int>
JetTaggingUtils::get_ltag(ROOT::VecOps::RVec<int> in,
                          float efficiency, float mistag_b,
                          float mistag_c, float mistag_g) {

  ROOT::VecOps::RVec<int> result(in.size(),0);

  for (size_t j = 0; j < in.size(); ++j) {
    if (in.at(j) <  4  && gRandom->Uniform() <= efficiency) result[j] = 1;
    if (in.at(j) ==  5 && gRandom->Uniform() <= mistag_b) result[j] = 1;
    if (in.at(j) ==  4 && gRandom->Uniform() <= mistag_c) result[j] = 1;
    if (in.at(j) == 21 && gRandom->Uniform() <= mistag_g) result[j] = 1;
  }
  return result;
}

ROOT::VecOps::RVec<int>
JetTaggingUtils::get_gtag(ROOT::VecOps::RVec<int> in,
                          float efficiency, float mistag_b,
                          float mistag_c, float mistag_l) {

  ROOT::VecOps::RVec<int> result(in.size(),0);

  for (size_t j = 0; j < in.size(); ++j) {
    if (in.at(j) == 21 && gRandom->Uniform() <= efficiency) result[j] = 1;
    if (in.at(j) ==  5 && gRandom->Uniform() <= mistag_b) result[j] = 1;
    if (in.at(j) ==  4 && gRandom->Uniform() <= mistag_c) result[j] = 1;
    if (in.at(j)  <  4 && gRandom->Uniform() <= mistag_l) result[j] = 1;
  }
  return result;
}

JetTaggingUtils::sel_tag::sel_tag(bool arg_pass): m_pass(arg_pass) {};
ROOT::VecOps::RVec<fastjet::PseudoJet>
JetTaggingUtils::sel_tag::operator()(ROOT::VecOps::RVec<bool> tags,
                                     ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<fastjet::PseudoJet> result;
  for (size_t i = 0; i < in.size(); ++i) {
    if (m_pass) {
      if (tags.at(i)) result.push_back(in.at(i));
    }
    else {
      if (!tags.at(i)) result.push_back(in.at(i));
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
    if (parton.generatorStatus>80 || parton.generatorStatus<70) continue;
    //Select outgoing partons from the hardest process - use ONLY(?) for Zqq events
    //(For Zqq, ideally, don't even need anything else, just identify the partons)
    //if (parton.generatorStatus!=23) continue;
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

ROOT::VecOps::RVec<int> JetTaggingUtils::get_flavour_gm_manual(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<float> pdg_gm){
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

  // push back zeros for all the particles except for ghosts
  ROOT::VecOps::RVec<float> pdg_gm(PJin.size(),0);
  // push back the MC pdg ID for all the ghosts
  for (size_t j = 0; j < MCin.size(); ++j) {
    auto & parton = MCin[j];
    // CAUTION: use the SAME selection here as in addGhosts_pseudoJets
    if (parton.generatorStatus!=23) continue;
    //if (parton.generatorStatus>80 || parton.generatorStatus<70) continue;
    if (parton.PDG > 5) continue;                     // only partons
    //if (parton.PDG > 5 && parton.PDG != 21) continue; // partons + gluons  
    pdg_gm.push_back(parton.PDG);
  }
  
  ROOT::VecOps::RVec<int> result(in.size(),0);
  if(in.size() == 0) return result;

  for (size_t i = 0; i < in.size(); i++) {
    auto & p = in[i];

    for (int ele : inJC.at(i)) {
      if (pdg_gm.at(ele) == 0) continue;
      if (abs(pdg_gm.at(ele)) > abs(result[i])) result[i] = pdg_gm.at(ele);

      // adding gluons
      //if(find(pdg_gm.begin(), pdg_gm.end(), 5) != pdg_gm.end())
    }
  }
  return result;
}

ROOT::VecOps::RVec<int> JetTaggingUtils::get_flavour_gm7x_auto(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin, std::vector<fastjet::PseudoJet> PJin){

  // pseudoJet vector is the one before addign ghosts
  // get the pdg_gm vector here (instead of in a separate fn)
  // if pseudoJet doesn't work, replace it with (not as general) reco particles

  // push back zeros for all the particles except for ghosts
  ROOT::VecOps::RVec<float> pdg_gm(PJin.size(),0);
  // push back the MC pdg ID for all the ghosts
  for (size_t j = 0; j < MCin.size(); ++j) {
    auto & parton = MCin[j];
    // CAUTION: use the SAME selection here as in addGhosts_pseudoJets
    //if (parton.generatorStatus!=23) continue;
    if (parton.generatorStatus>80 || parton.generatorStatus<70) continue;
    if (parton.PDG > 5) continue;                     // only partons
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

ROOT::VecOps::RVec<int> JetTaggingUtils::get_flavour_gm_pcut(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin, std::vector<fastjet::PseudoJet> PJin, float p_cut) {

  // push back zeros for all the particles except for ghosts
  ROOT::VecOps::RVec<float> pdg_gm(PJin.size(),0);
  TLorentzVector zero(0,0,0,0);
  ROOT::VecOps::RVec<TLorentzVector> p4_gm(PJin.size(),zero);
  // push back the MC pdg ID for all the ghosts
  for (size_t j = 0; j < MCin.size(); ++j) {
    auto & parton = MCin[j];
    // CAUTION: use the SAME selection here as in addGhosts_pseudoJets
    //if (parton.generatorStatus!=23) continue;
    if (parton.generatorStatus>80 || parton.generatorStatus<70) continue;
    if (parton.PDG > 5) continue;                     // only partons
    pdg_gm.push_back(parton.PDG);
    //
    TLorentzVector p4;
    p4.SetXYZM(parton.momentum.x, parton.momentum.y, parton.momentum.z, parton.mass);
    p4_gm.push_back(p4);
  }
  
  ROOT::VecOps::RVec<int> result(in.size(),0);
  if(in.size() == 0) return result;

  for (size_t i = 0; i < in.size(); i++) {
    auto & p = in[i];

    for (int ele : inJC.at(i)) {
      if (pdg_gm.at(ele) == 0) continue;
      if (p4_gm.at(ele).P() < p_cut) continue; // apply momentum cut
      if (abs(pdg_gm.at(ele)) > abs(result[i])) result[i] = pdg_gm.at(ele);
    }
  }
  return result;
}

std::vector<std::vector<int>> JetTaggingUtils::get_flavour_gm(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin, std::vector<fastjet::PseudoJet> PJin, int statCode, float p_cut) {

  // CAUTION: use the SAME statCode here as in addGhosts_pseudoJets
    
  // push back zeros (later update to NaN) for all the particles except for ghosts
  ROOT::VecOps::RVec<float> pdg_gm(PJin.size(),0);
  TLorentzVector zero(0,0,0,0);
  ROOT::VecOps::RVec<TLorentzVector> p4_gm(PJin.size(),zero);
  // push back the MC pdg ID for all the ghosts
  for (size_t j = 0; j < MCin.size(); ++j) {
    auto & parton = MCin[j];

    TLorentzVector p4;
    
    // statCode==0 : select outgoing particles from hardest reaction
    if (statCode==0 && parton.generatorStatus==23) {
      //if (parton.PDG > 5) continue;                     // only partons
      pdg_gm.push_back(parton.PDG);
      //
      p4.SetXYZM(parton.momentum.x, parton.momentum.y, parton.momentum.z, parton.mass);
      p4_gm.push_back(p4);
    }
    // statCode==1 : select partons just before hadronisation
    else if (statCode==1 && parton.generatorStatus<80 && parton.generatorStatus>70) {
      //if (parton.PDG > 5) continue;                     // only partons
      pdg_gm.push_back(parton.PDG);
      //
      p4.SetXYZM(parton.momentum.x, parton.momentum.y, parton.momentum.z, parton.mass);
      p4_gm.push_back(p4);
    }
    
    // select primary hadrons after hadronisation
    if (parton.generatorStatus>80 && parton.generatorStatus<90)
      {
	pdg_gm.push_back(parton.PDG);
	//
	p4.SetXYZM(parton.momentum.x, parton.momentum.y, parton.momentum.z, parton.mass);
	p4_gm.push_back(p4);
      }
  }

  // initiate with 0 because unassigned jets are given flvour=0 for now
  std::vector<int> partonFlv(in.size(),0);
  std::vector<int> hadronFlv(in.size(),0); // currently only absolute flavour
  std::vector<std::vector<int>> result;
  if(in.size() == 0) return result;

  for (size_t i = 0; i < in.size(); i++) {
    auto & p = in[i];

    for (int ele : inJC.at(i)) {
      if (pdg_gm.at(ele) == 0) continue;
      if (p4_gm.at(ele).P() < p_cut) continue; // apply momentum cut (p_cut=0 for no cut)

      // hadron flavour
      if (abs(int(pdg_gm.at(ele)/1000)%10) == 5 || abs(int(pdg_gm.at(ele)/100)%10) == 5) hadronFlv[i] = 5;
      else if (abs(int(pdg_gm.at(ele)/1000)%10) == 4 || abs(int(pdg_gm.at(ele)/100)%10) == 4) hadronFlv[i] = 4;
      else if (abs(int(pdg_gm.at(ele)/1000)%10) == 3 || abs(int(pdg_gm.at(ele)/100)%10) == 3 || abs(pdg_gm.at(ele)) == 130) hadronFlv[i] = 3;
      // can be updated to only 1 flavour vector as an output instead of separate hadron and parton flavour vectors
      // in that case, if no b- or c-hadron is found, start searching for partons

      // parton flavour
      float parton_mom = 0;
      if (abs(pdg_gm.at(ele)) == 5 && p4_gm.at(ele).P() > parton_mom) {
	partonFlv[i] = pdg_gm.at(ele);
	parton_mom = p4_gm.at(ele).P(); }
      else if (abs(pdg_gm.at(ele)) == 4 && p4_gm.at(ele).P() > parton_mom) {
	partonFlv[i] = pdg_gm.at(ele);
	parton_mom = p4_gm.at(ele).P(); }
      else if (p4_gm.at(ele).P() > parton_mom) {
	partonFlv[i] = pdg_gm.at(ele);
	parton_mom = p4_gm.at(ele).P(); }
    }
  }
  //
  result.push_back(hadronFlv);
  result.push_back(partonFlv);
  return result;
}
=======
