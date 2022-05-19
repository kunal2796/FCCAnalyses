#include "JetClusteringUtils.h"
using namespace JetClusteringUtils;


ROOT::VecOps::RVec<fastjet::PseudoJet> JetClusteringUtils::get_pseudoJets(FCCAnalysesJet jets){
  return jets.jets;
}

std::vector<std::vector<int>> JetClusteringUtils::get_constituents(FCCAnalysesJet jets){
  return jets.constituents;
}

std::vector<fastjet::PseudoJet> JetClusteringUtils::set_pseudoJets(ROOT::VecOps::RVec<float> px, 
								   ROOT::VecOps::RVec<float> py, 
								   ROOT::VecOps::RVec<float> pz, 
								   ROOT::VecOps::RVec<float> e) {
  std::vector<fastjet::PseudoJet> result;
  unsigned index = 0;
  for (size_t i = 0; i < px.size(); ++i) {
    result.emplace_back(px.at(i), py.at(i), pz.at(i), e.at(i));
    result.back().set_user_index(index);
    ++index;
  }
  return result;
}

std::vector<fastjet::PseudoJet> JetClusteringUtils::addMore_pseudoJets(std::vector<fastjet::PseudoJet> pseudoJ,
								       ROOT::VecOps::RVec<float> px, 
								       ROOT::VecOps::RVec<float> py, 
								       ROOT::VecOps::RVec<float> pz, 
								       ROOT::VecOps::RVec<float> e) {
  unsigned index = pseudoJ.size();
  for (size_t i = 0; i < px.size(); ++i) {
    pseudoJ.emplace_back(px.at(i), py.at(i), pz.at(i), e.at(i));
    pseudoJ.back().set_user_index(index);
    ++index;
  }
  return pseudoJ;
}

std::vector<fastjet::PseudoJet> JetClusteringUtils::addGhosts_pseudoJets_old(std::vector<fastjet::PseudoJet> pseudoJ,
									     ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin) {
  unsigned index = pseudoJ.size();
  for (size_t i = 0; i < MCin.size(); ++i) {
    auto & parton = MCin[i];
    // select outgoing partons from the hardest reaction
    // later introduce a way for user to choose from the 2 status code options
    if (parton.generatorStatus!=23) continue;
    //if (parton.generatorStatus>80 || parton.generatorStatus<70) continue;
    if (parton.PDG > 5) continue;                     // only partons
    //if (parton.PDG > 5 && parton.PDG != 21) continue; // partons + gluons

    TLorentzVector tlv;
    tlv.SetXYZM(parton.momentum.x, parton.momentum.y, parton.momentum.z, parton.mass);
    pseudoJ.emplace_back(parton.momentum.x * 1.e-18, parton.momentum.y * 1.e-18, parton.momentum.z * 1.e-18, tlv.E() * 1.e-18);
    pseudoJ.back().set_user_index(index);
    ++index;
  }
  return pseudoJ;
}

std::vector<fastjet::PseudoJet> JetClusteringUtils::addGhosts_pseudoJets(std::vector<fastjet::PseudoJet> pseudoJ,
									 ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin,
									 int statCode) {
  unsigned index = pseudoJ.size();
  for (size_t i = 0; i < MCin.size(); ++i) {
    auto & parton = MCin[i];

    TLorentzVector tlv;
    tlv.SetXYZM(parton.momentum.x, parton.momentum.y, parton.momentum.z, parton.mass);

    // statCode==0 : select outgoing particles from hardest reaction
    if (statCode==0 && parton.generatorStatus==23) {
      //if (parton.PDG > 5) continue;                     // only partons
      pseudoJ.emplace_back(parton.momentum.x * 1.e-18, parton.momentum.y * 1.e-18, parton.momentum.z * 1.e-18, tlv.E() * 1.e-18);
      pseudoJ.back().set_user_index(index);
      ++index;
    }
    // statCode==1 : select partons just before hadronisation
    else if (statCode==1 && parton.generatorStatus<80 && parton.generatorStatus>70) {
      //if (parton.PDG > 5) continue;                     // only partons
      pseudoJ.emplace_back(parton.momentum.x * 1.e-18, parton.momentum.y * 1.e-18, parton.momentum.z * 1.e-18, tlv.E() * 1.e-18);
      pseudoJ.back().set_user_index(index);
      ++index;
    }
    
    // select primary hadrons after hadronisation
    if (parton.generatorStatus>80 && parton.generatorStatus<90)
      {
	pseudoJ.emplace_back(parton.momentum.x * 1.e-18, parton.momentum.y * 1.e-18, parton.momentum.z * 1.e-18, tlv.E() * 1.e-18);
	pseudoJ.back().set_user_index(index);
	++index;
      }
  }
  return pseudoJ;
}

std::vector<fastjet::PseudoJet> JetClusteringUtils::addGhosts7x_pseudoJets(std::vector<fastjet::PseudoJet> pseudoJ,
									   ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin) {
  unsigned index = pseudoJ.size();
  for (size_t i = 0; i < MCin.size(); ++i) {
    auto & parton = MCin[i];
    // select outgoing partons from the hardest reaction
    // later introduce a way for user to choose from the 2 status code options
    //if (parton.generatorStatus!=23) continue;
    if (parton.generatorStatus>80 || parton.generatorStatus<70) continue;
    if (parton.PDG > 5) continue;                     // only partons
    
    TLorentzVector tlv;
    tlv.SetXYZM(parton.momentum.x, parton.momentum.y, parton.momentum.z, parton.mass);
    pseudoJ.emplace_back(parton.momentum.x * 1.e-18, parton.momentum.y * 1.e-18, parton.momentum.z * 1.e-18, tlv.E() * 1.e-18);
    pseudoJ.back().set_user_index(index);
    ++index;
  }
  return pseudoJ;
}

ROOT::VecOps::RVec<float> JetClusteringUtils::get_gmPDG(ROOT::VecOps::RVec<float> pdg, 
							ROOT::VecOps::RVec<float> px, 
							ROOT::VecOps::RVec<float> px_g) {
  // the arguments can later be changed to accept size (get_n) instead of the momentum vectors
  // push back zeros for all the reco particles
  ROOT::VecOps::RVec<float> result(px.size(),0);
  // push back the MC pdg ID for all the ghosts
  for (size_t j = 0; j < px_g.size(); ++j) {
    result.push_back(pdg.at(j));
  }
  return result;
}

std::vector<fastjet::PseudoJet> JetClusteringUtils::set_pseudoJets_xyzm(ROOT::VecOps::RVec<float> px, 
								   ROOT::VecOps::RVec<float> py, 
								   ROOT::VecOps::RVec<float> pz, 
								   ROOT::VecOps::RVec<float> m) {
  std::vector<fastjet::PseudoJet> result;
  unsigned index = 0;
  for (size_t i = 0; i < px.size(); ++i) {
    double px_d = px.at(i);
    double py_d = py.at(i);
    double pz_d = pz.at(i);
    double  m_d =  m.at(i);
    double  E_d = sqrt(px_d*px_d + py_d*py_d + pz_d*pz_d + m_d*m_d);
    result.emplace_back(px_d, py_d, pz_d, E_d);
    result.back().set_user_index(index);
    ++index;
  }
  return result;
}


ROOT::VecOps::RVec<float> JetClusteringUtils::get_px(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.px());
  }
  return result;
}


ROOT::VecOps::RVec<float> JetClusteringUtils::get_py(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.py());
  }
  return result;
}

ROOT::VecOps::RVec<float> JetClusteringUtils::get_pz(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.pz());
  }
  return result;
}

ROOT::VecOps::RVec<float> JetClusteringUtils::get_e(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.E());
  }
  return result;
}

ROOT::VecOps::RVec<float> JetClusteringUtils::get_p(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    float mom = sqrt(p.px()*p.px() + p.py()*p.py() + p.pz()*p.pz());
    result.push_back(mom);
  }
  return result;
}

ROOT::VecOps::RVec<float> JetClusteringUtils::get_pt(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.pt());
  }
  return result;
}

ROOT::VecOps::RVec<float> JetClusteringUtils::get_m(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.m());
  }
  return result;
}

ROOT::VecOps::RVec<float> JetClusteringUtils::get_eta(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.eta());
  }
  return result;
}

ROOT::VecOps::RVec<float> JetClusteringUtils::get_phi(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.phi());
  }
  return result;
}

ROOT::VecOps::RVec<float> JetClusteringUtils::get_theta(ROOT::VecOps::RVec<fastjet::PseudoJet> in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.theta());
  }
  return result;
}


/// ------ ///
// From Edi //

ROOT::VecOps::RVec<float> JetClusteringUtils::get_nConstituents(std::vector<std::vector<int>> constituents){
  ROOT::VecOps::RVec<float> result;
  for (auto& constis : constituents){
    result.push_back(constis.size());
  }
  return result;
}

//Note that the below doesn't make sense atm, there is no consti vec of pseudojets (to my knowledge, maybe this refers to an old implementation)...
///ROOT::VecOps::RVec<float> JetClusteringUtils::get_dTheta(ROOT::VecOps::RVec<fastjet::PseudoJet> in, ROOT::VecOps::RVec<fastjet::PseudoJet> constituents, std::vector<std::vector<int>> indices){
///  ROOT::VecOps::RVec<float> result(constituents.size(), 0);
///  //for (auto& ind : indices){
///  for (int j=0; j<indices.size(); j++){
///    for (auto& i : indices[j]){
///    //for (int i=0; i<contituents.size())
///      result[i] = in[j].theta()-constituents[i].theta();
///    }
///  }
///  return result;
///}

std::vector<std::vector<float>> JetClusteringUtils::get_dTheta(ROOT::VecOps::RVec<float> jet_theta, std::vector<std::vector<float>> constituents_theta){
  std::vector<std::vector<float>> result;
  //for (auto& ind : indices){
  for (int j=0; j<jet_theta.size(); j++){
    std::vector<float> tmp_res;
    for (auto& consti_theta : constituents_theta[j]){
      tmp_res.push_back(jet_theta[j]-consti_theta);
    }
    result.push_back(tmp_res);
  }
  return result;
}

std::vector<std::vector<float>> JetClusteringUtils::get_dPhi(ROOT::VecOps::RVec<float> jet_phi, std::vector<std::vector<float>> constituents_phi){
  //const float  PI = 3.14159265358979f;
  float dphi; 
  std::vector<std::vector<float>> result;
  for (int j=0; j<jet_phi.size(); j++){
    std::vector<float> tmp_res;
    for (auto& consti_phi : constituents_phi[j]){
      dphi = jet_phi[j]-consti_phi;
      if(std::abs(dphi)<=(fastjet::pi)){
        tmp_res.push_back(dphi);
      }
      else if(dphi>fastjet::pi){
        tmp_res.push_back(dphi-2*fastjet::pi);
      }
      else if(dphi<-fastjet::pi){
        tmp_res.push_back(dphi+2*fastjet::pi);
      }
    }
    result.push_back(tmp_res);
  }
  return result;
}

std::vector<std::vector<float>> JetClusteringUtils::get_pRel(ROOT::VecOps::RVec<float> jet_p, std::vector<std::vector<float>> constituents_p){
  std::vector<std::vector<float>> result;
  //for (auto& ind : indices){
  for (int j=0; j<jet_p.size(); j++){
    std::vector<float> tmp_res;
    for (auto& consti_p : constituents_p[j]){
      tmp_res.push_back(consti_p/jet_p[j]);
    }
    result.push_back(tmp_res);
  }
  return result;
}

std::vector<std::vector<float>> JetClusteringUtils::reshape2jet(ROOT::VecOps::RVec<float> var, std::vector<std::vector<int>> constituents){
  std::vector<std::vector<float>> result;
  for (auto& ind : constituents){
    //ROOT::VecOps::RVec<float> tmp_res(ind.size(), 0);
    std::vector<float> tmp_res;
    for (auto& i : ind){
      tmp_res.push_back(var.at(i));
    }
    result.push_back(tmp_res);
  }
  return result;
}

/// ------ ///


FCCAnalysesJet JetClusteringUtils::initialise_FCCAnalysesJet(){
  
  JetClusteringUtils::FCCAnalysesJet result;
  ROOT::VecOps::RVec<fastjet::PseudoJet> jets;
  std::vector<std::vector<int>> constituents;

  result.jets = jets;
  result.constituents = constituents;

  return result;
};

FCCAnalysesJet JetClusteringUtils::build_FCCAnalysesJet(std::vector<fastjet::PseudoJet> in){
  JetClusteringUtils::FCCAnalysesJet result = JetClusteringUtils::initialise_FCCAnalysesJet();
  for (const auto& pjet : in) {
    result.jets.push_back(pjet);
    
    std::vector<fastjet::PseudoJet> consts = pjet.constituents();
    std::vector<int> tmpvec;
    for (const auto& constituent : consts){
      tmpvec.push_back(constituent.user_index());  
    }
    result.constituents.push_back(tmpvec);
  }
  return result;
}



std::vector<fastjet::PseudoJet> JetClusteringUtils::build_jets(fastjet::ClusterSequence & cs, int exclusive, float cut, int sorted){
  std::vector<fastjet::PseudoJet> pjets;

  if (sorted == 0){
    if(exclusive ==  0 )       pjets = fastjet::sorted_by_pt(cs.inclusive_jets(cut));
    else if( exclusive ==  1)  pjets = fastjet::sorted_by_pt(cs.exclusive_jets(cut));
    else if( exclusive ==  2)  pjets = fastjet::sorted_by_pt(cs.exclusive_jets(int(cut)));
    else if( exclusive ==  3)  pjets = fastjet::sorted_by_pt(cs.exclusive_jets_up_to(int(cut)));
    else if( exclusive ==  4)  pjets = fastjet::sorted_by_pt(cs.exclusive_jets_ycut(cut));
  }
  else if (sorted == 1){
    if(exclusive ==  0 )       pjets = fastjet::sorted_by_E(cs.inclusive_jets(cut));
    else if( exclusive ==  1)  pjets = fastjet::sorted_by_E(cs.exclusive_jets(cut));
    else if( exclusive ==  2)  pjets = fastjet::sorted_by_E(cs.exclusive_jets(int(cut)));
    else if( exclusive ==  3)  pjets = fastjet::sorted_by_E(cs.exclusive_jets_up_to(int(cut)));
    else if( exclusive ==  4)  pjets = fastjet::sorted_by_E(cs.exclusive_jets_ycut(cut));
  }
  return pjets;
}


bool JetClusteringUtils::check(unsigned int n, int exclusive, float cut){
  //if (exclusive>0 && n<=int(cut)) return false;
  if (exclusive>0 && n<int(cut)) return false; //atleast for exclusive 2 & 3
  return true;
}

fastjet::RecombinationScheme JetClusteringUtils::recomb_scheme(int recombination){
  fastjet::RecombinationScheme recomb_scheme;

  if(recombination == 0) recomb_scheme = fastjet::RecombinationScheme::E_scheme;
  else if (recombination == 1) recomb_scheme = fastjet::RecombinationScheme::pt_scheme;
  else if (recombination == 2) recomb_scheme = fastjet::RecombinationScheme::pt2_scheme;
  else if (recombination == 3) recomb_scheme = fastjet::RecombinationScheme::Et_scheme;
  else if (recombination == 4) recomb_scheme = fastjet::RecombinationScheme::Et2_scheme;
  else if (recombination == 5) recomb_scheme = fastjet::RecombinationScheme::BIpt_scheme;
  else if (recombination == 6) recomb_scheme = fastjet::RecombinationScheme::BIpt2_scheme;
  else recomb_scheme = fastjet::RecombinationScheme::external_scheme;

  return recomb_scheme;
}
