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

ROOT::VecOps::RVec<float> get_p(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in){
  ROOT::VecOps::RVec<float> result;
  for (auto & part: in) {
    double px = part.px();
    double py = part.py();
    double pz = part.pz();
    double p = sqrt(px*px + py*py + pz*pz);
    result.push_back(p);
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



//FCCAnalysesJet initialise_FCCAnalysesJet(){
//
//  FCCAnalysesJet result;
//  ROOT::VecOps::RVec<fastjet::PseudoJet> jets;
//  std::vector<std::vector<int>> constituents;
//
//  result.jets = jets;
//  result.constituents = constituents;
//
//  std::vector<float> exclusive_dmerge;
//  std::vector<float> exclusive_dmerge_max;
//  exclusive_dmerge.reserve(Nmax_dmerge);
//  exclusive_dmerge_max.reserve(Nmax_dmerge);
//
//  result.exclusive_dmerge = exclusive_dmerge;
//  result.exclusive_dmerge_max = exclusive_dmerge_max;
//
//  return result;
//};



FCCAnalysesJet initialise_FCCAnalysesJet(TString clustering_algo, ROOT::VecOps::RVec<float> clustering_params){

  FCCAnalysesJet result;
  ROOT::VecOps::RVec<fastjet::PseudoJet> jets;
  std::vector<std::vector<int>> constituents;

  result.jets = jets;
  result.constituents = constituents;
  result.clustering_algo = clustering_algo;
  result.clustering_params = clustering_params;

  std::vector<float> exclusive_dmerge;
  std::vector<float> exclusive_dmerge_max;
  exclusive_dmerge.reserve(Nmax_dmerge);
  exclusive_dmerge_max.reserve(Nmax_dmerge);

  result.exclusive_dmerge = exclusive_dmerge;
  result.exclusive_dmerge_max = exclusive_dmerge_max;

  return result;
};

//FCCAnalysesJet build_FCCAnalysesJet(const std::vector<fastjet::PseudoJet> &in,
//                                                        const std::vector<float> &dmerge,
//                                                        const std::vector<float> &dmerge_max){
//
//  FCCAnalysesJet result = initialise_FCCAnalysesJet();
//  for (const auto& pjet : in) {
//    result.jets.push_back(pjet);
//
//    std::vector<fastjet::PseudoJet> consts = pjet.constituents();
//    std::vector<int> tmpvec;
//    for (const auto& constituent : consts){
//      tmpvec.push_back(constituent.user_index());
//    }
//    result.constituents.push_back(tmpvec);
//  }
//  result.exclusive_dmerge = dmerge;
//  result.exclusive_dmerge_max = dmerge_max;
//  return result;
//}


FCCAnalysesJet build_FCCAnalysesJet(const std::vector<fastjet::PseudoJet> &in,
                                                        const std::vector<float> &dmerge,
                                                        const std::vector<float> &dmerge_max,
                                                        TString clustering_algo,
                                                        ROOT::VecOps::RVec<float> clustering_params){

  FCCAnalysesJet result = initialise_FCCAnalysesJet(clustering_algo, clustering_params);
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
