#include "FCCAnalyses/JetClusteringUtils.h"

namespace FCCAnalyses{

namespace JetClusteringUtils{


ROOT::VecOps::RVec<fastjet::PseudoJet> get_pseudoJets(const FCCAnalysesJet &jets){
  return jets.jets;
}

std::vector<std::vector<int>> get_constituents(const FCCAnalysesJet &jets){
  return jets.constituents;
}


float get_exclusive_dmerge(const FCCAnalysesJet &in,
                                               int n) {
  float d = -1;
  if ( n >= 1 &&  n <= Nmax_dmerge) d= in.exclusive_dmerge[n-1] ;
  return d;
}

float get_exclusive_dmerge_max(const FCCAnalysesJet &in,
                                                   int n) {
  float d = -1;
  if ( n >= 1 &&  n <= Nmax_dmerge) d= in.exclusive_dmerge_max[n-1] ;
  return d;
}

std::vector<fastjet::PseudoJet> set_pseudoJets(const ROOT::VecOps::RVec<float> &px,
					       const ROOT::VecOps::RVec<float> &py,
					       const ROOT::VecOps::RVec<float> &pz,
					       const ROOT::VecOps::RVec<float> &e){
  std::vector<fastjet::PseudoJet> result;
  unsigned index = 0;
  for (size_t i = 0; i < px.size(); ++i) {
    result.emplace_back(px[i], py[i], pz[i], e[i]);
    result.back().set_user_index(index);
    ++index;
  }
  return result;
}

std::vector<fastjet::PseudoJet> set_pseudoJets_xyzm(const ROOT::VecOps::RVec<float> &px,
                                                                        const ROOT::VecOps::RVec<float> &py,
                                                                        const ROOT::VecOps::RVec<float> &pz,
                                                                        const ROOT::VecOps::RVec<float> &m){
  std::vector<fastjet::PseudoJet> result;
  unsigned index = 0;
  for (size_t i = 0; i < px.size(); ++i) {
    double px_d = px[i];
    double py_d = py[i];
    double pz_d = pz[i];
    double  m_d =  m[i];
    double  E_d = sqrt(px_d*px_d + py_d*py_d + pz_d*pz_d + m_d*m_d);
    result.emplace_back(px_d, py_d, pz_d, E_d);
    result.back().set_user_index(index);
    ++index;
  }
  return result;
}


std::vector<fastjet::PseudoJet> set_ghostPseudoJets_xyzm_primitive(const ROOT::VecOps::RVec<float> px, 
								   const ROOT::VecOps::RVec<float> py, 
								   const ROOT::VecOps::RVec<float> pz, 
								   const ROOT::VecOps::RVec<float> m, 
								   const ROOT::VecOps::RVec<float> genStatus, 
								   const ROOT::VecOps::RVec<float> px_MC, 
								   const ROOT::VecOps::RVec<float> py_MC, 
								   const ROOT::VecOps::RVec<float> pz_MC, 
								   const ROOT::VecOps::RVec<float> m_MC) {
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
  for (size_t i = 0; i < px_MC.size(); ++i) {
    bool flag1 = false;
    bool flag2 = false;
    if (((genStatus[i]>29)) || (genStatus[i]<21)) continue;
    //if ((genStatus[i]>89) || (genStatus[i]<81)) flag1 = true;
    //if (((genStatus[i]>29)) || (genStatus[i]<21)) flag2 = true;
    //if (flag1 && flag2) continue;
    double px_d = px_MC.at(i);
    double py_d = py_MC.at(i);
    double pz_d = pz_MC.at(i);
    double  m_d =  m_MC.at(i);
    double  E_d = sqrt(px_d*px_d + py_d*py_d + pz_d*pz_d + m_d*m_d);
    result.emplace_back(px_d*pow(10, -18), py_d*pow(10, -18), pz_d*pow(10, -18), E_d*pow(10, -18));
    result.back().set_user_index(index);
    ++index;
  }
  return result;
}

std::vector<fastjet::PseudoJet> set_ghostPseudoJets_xyzm(std::vector<fastjet::PseudoJet> pseudoJets, ROOT::VecOps::RVec<edm4hep::MCParticleData> ghosts) {
  //std::vector<fastjet::PseudoJet> result;
  unsigned index = pseudoJets.size();

  std::vector<int> c_hadrons = {411, 421};
  std::vector<int> b_hadrons = {511, 521};


  for (auto& p : ghosts) {
    bool isGhost = false;
    //if ((p.generatorStatus<80) && (p.generatorStatus>70)) {
    if ((p.generatorStatus<30) && (p.generatorStatus>20)) {
        isGhost = true;
    }
    else if (p.generatorStatus == 1) {
        if (std::find(c_hadrons.begin(), c_hadrons.end(), std::abs(int(p.PDG))) != c_hadrons.end()) {
            isGhost = true;
        }
        else if (std::find(b_hadrons.begin(), b_hadrons.end(), std::abs(int(p.PDG))) != b_hadrons.end()) {
            isGhost = true;
        }
    }
    if (isGhost){
        double px = p.momentum.x;//*pow(10, -1);
        double py = p.momentum.y;
        double pz = p.momentum.z;
        double m = p.mass;
        double  E = sqrt(px*px + py*py + pz*pz + m*m);
        pseudoJets.emplace_back(px*pow(10, -18), py*pow(10, -18), pz*pow(10, -18), E*pow(10, -18));
        //std::cout << px << std::endl;
        //pseudoJets.emplace_back(px, py, pz, E);
        pseudoJets.back().set_user_index(index);
        ++index;
    }
  }
  return pseudoJets;
}

  

ROOT::VecOps::RVec<float> get_px(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.px());
  }
  return result;
}


ROOT::VecOps::RVec<float> get_py(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.py());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_pz(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.pz());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_e(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.E());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_pt(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.pt());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_m(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.m());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_eta(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.eta());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_phi(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.phi());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_theta(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in){
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.theta());
  }
  return result;
}

// Ghost Matching
std::vector<fastjet::PseudoJet> addMore_pseudoJets(std::vector<fastjet::PseudoJet> pseudoJ,
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

std::vector<fastjet::PseudoJet> addGhosts_pseudoJets_old(std::vector<fastjet::PseudoJet> pseudoJ,
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

std::vector<fastjet::PseudoJet> addGhosts_pseudoJets(std::vector<fastjet::PseudoJet> pseudoJ,
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

std::vector<fastjet::PseudoJet> addGhosts7x_pseudoJets(std::vector<fastjet::PseudoJet> pseudoJ,
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

ROOT::VecOps::RVec<float> get_gmPDG(ROOT::VecOps::RVec<float> pdg, 
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

/// ------ ///
// From Edi //

ROOT::VecOps::RVec<float> get_nConstituents(std::vector<std::vector<int>> constituents){
  ROOT::VecOps::RVec<float> result;
  for (auto& constis : constituents){
    result.push_back(constis.size());
  }
  return result;
}

std::vector<std::vector<float>> get_dTheta(ROOT::VecOps::RVec<float> jet_theta, std::vector<std::vector<float>> constituents_theta){
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

std::vector<std::vector<float>> get_dPhi(ROOT::VecOps::RVec<float> jet_phi, std::vector<std::vector<float>> constituents_phi){
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

std::vector<std::vector<float>> get_pRel(ROOT::VecOps::RVec<float> jet_p, std::vector<std::vector<float>> constituents_p){
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

std::vector<std::vector<float>> reshape2jet(ROOT::VecOps::RVec<float> var, std::vector<std::vector<int>> constituents){
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


FCCAnalysesJet initialise_FCCAnalysesJet(){

  FCCAnalysesJet result;
  ROOT::VecOps::RVec<fastjet::PseudoJet> jets;
  std::vector<std::vector<int>> constituents;

  result.jets = jets;
  result.constituents = constituents;

  std::vector<float> exclusive_dmerge;
  std::vector<float> exclusive_dmerge_max;
  exclusive_dmerge.reserve(Nmax_dmerge);
  exclusive_dmerge_max.reserve(Nmax_dmerge);

  result.exclusive_dmerge = exclusive_dmerge;
  result.exclusive_dmerge_max = exclusive_dmerge_max;

  return result;
};

FCCAnalysesJet build_FCCAnalysesJet(const std::vector<fastjet::PseudoJet> &in,
                                                        const std::vector<float> &dmerge,
                                                        const std::vector<float> &dmerge_max){

  FCCAnalysesJet result = initialise_FCCAnalysesJet();
  for (const auto& pjet : in) {
    result.jets.push_back(pjet);

    std::vector<fastjet::PseudoJet> consts = pjet.constituents();
    std::vector<int> tmpvec;
    for (const auto& constituent : consts){
      tmpvec.push_back(constituent.user_index());
    }
    result.constituents.push_back(tmpvec);
  }
  result.exclusive_dmerge = dmerge;
  result.exclusive_dmerge_max = dmerge_max;
  return result;
}



std::vector<fastjet::PseudoJet> build_jets(fastjet::ClusterSequence &cs,
                                                               int exclusive,
                                                               float cut,
                                                               int sorted){
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

std::vector<float> exclusive_dmerge(fastjet::ClusterSequence &cs,
                                                        int do_dmarge_max) {

  const int Nmax = Nmax_dmerge;
  std::vector<float>  result;
  for (int i=1; i <= Nmax; i++) {
     	float  d;
	const int j = i;
     	if ( do_dmarge_max == 0) d = cs.exclusive_dmerge( j );
	else d = cs.exclusive_dmerge_max( j ) ;
	result.push_back( d );
  }
  return result;
}


bool check(unsigned int n,
                               int exclusive,
                               float cut){
  if (exclusive>0 && n<=int(cut)) return false;
  return true;
}

fastjet::RecombinationScheme recomb_scheme(int recombination){
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

}//end NS JetClusteringUtils

}//end NS FCCAnalyses
