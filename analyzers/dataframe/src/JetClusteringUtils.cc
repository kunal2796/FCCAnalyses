#include "FCCAnalyses/JetClusteringUtils.h"
#include "TLorentzVector.h"

namespace FCCAnalyses {
namespace JetClusteringUtils {

ROOT::VecOps::RVec<fastjet::PseudoJet>
get_pseudoJets(const JetClustering::FCCAnalysesJet &jets) {
  return jets.jets;
}

std::vector<std::vector<int>>
get_constituents(const JetClustering::FCCAnalysesJet &jets) {
  return jets.constituents;
}

float get_exclusive_dmerge(const JetClustering::FCCAnalysesJet &in, int n) {
  float d = -1;
  if (n >= 1 && n <= Nmax_dmerge && in.exclusive_dmerge.size() > n - 1)
    d = in.exclusive_dmerge[n - 1];
  return d;
}

float get_exclusive_dmerge_max(const JetClustering::FCCAnalysesJet &in, int n) {
  float d = -1;
  if (n >= 1 && n <= Nmax_dmerge && in.exclusive_dmerge.size() > n - 1)
    d = in.exclusive_dmerge_max[n - 1];
  return d;
}

std::vector<fastjet::PseudoJet> set_pseudoJets(
    const ROOT::VecOps::RVec<float> &px, const ROOT::VecOps::RVec<float> &py,
    const ROOT::VecOps::RVec<float> &pz, const ROOT::VecOps::RVec<float> &e) {
  std::vector<fastjet::PseudoJet> result;
  unsigned index = 0;
  for (size_t i = 0; i < px.size(); ++i) {
    result.emplace_back(px[i], py[i], pz[i], e[i]);
    result.back().set_user_index(index);
    ++index;
  }
  return result;
}

std::vector<fastjet::PseudoJet> set_pseudoJets_xyzm(
    const ROOT::VecOps::RVec<float> &px, const ROOT::VecOps::RVec<float> &py,
    const ROOT::VecOps::RVec<float> &pz, const ROOT::VecOps::RVec<float> &m) {
  std::vector<fastjet::PseudoJet> result;
  unsigned index = 0;
  for (size_t i = 0; i < px.size(); ++i) {
    double px_d = px[i];
    double py_d = py[i];
    double pz_d = pz[i];
    double m_d = m[i];
    double E_d = sqrt(px_d * px_d + py_d * py_d + pz_d * pz_d + m_d * m_d);
    result.emplace_back(px_d, py_d, pz_d, E_d);
    result.back().set_user_index(index);
    ++index;
  }
  return result;
}

ROOT::VecOps::RVec<float>
get_px(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : in) {
    result.push_back(p.px());
  }
  return result;
}

ROOT::VecOps::RVec<float>
get_py(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : in) {
    result.push_back(p.py());
  }
  return result;
}

ROOT::VecOps::RVec<float>
get_pz(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : in) {
    result.push_back(p.pz());
  }
  return result;
}

ROOT::VecOps::RVec<float>
get_e(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : in) {
    result.push_back(p.E());
  }
  return result;
}

ROOT::VecOps::RVec<float>
get_pt(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : in) {
    result.push_back(p.pt());
  }
  return result;
}

ROOT::VecOps::RVec<float>
get_p(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : in) {
    result.push_back(sqrt(p.pt() * p.pt() + p.pz() * p.pz()));
  }
  return result;
}

ROOT::VecOps::RVec<float>
get_m(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : in) {
    result.push_back(p.m());
  }
  return result;
}

ROOT::VecOps::RVec<float>
get_eta(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : in) {
    result.push_back(p.eta());
  }
  return result;
}

ROOT::VecOps::RVec<float>
get_phi(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : in) {
    result.push_back(p.phi());
  }
  return result;
}

ROOT::VecOps::RVec<float>
get_phi_std(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : in) {
    result.push_back(p.phi_std());
  }
  return result;
}

ROOT::VecOps::RVec<float>
get_theta(const ROOT::VecOps::RVec<fastjet::PseudoJet> &in) {
  ROOT::VecOps::RVec<float> result;
  for (auto &p : in) {
    result.push_back(p.theta());
  }
  return result;
}

sel_pt::sel_pt(float arg_min_pt) : m_min_pt(arg_min_pt){};
ROOT::VecOps::RVec<fastjet::PseudoJet>
sel_pt::operator()(ROOT::VecOps::RVec<fastjet::PseudoJet> in) {
  ROOT::VecOps::RVec<fastjet::PseudoJet> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto &p = in[i];
    if (std::sqrt(std::pow(p.px(), 2) + std::pow(p.py(), 2)) > m_min_pt) {
      result.emplace_back(p);
    }
  }
  return result;
}

JetClustering::FCCAnalysesJet initialise_FCCAnalysesJet() {
  JetClustering::FCCAnalysesJet result;
  std::vector<fastjet::PseudoJet> jets;
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

JetClustering::FCCAnalysesJet
build_FCCAnalysesJet(const std::vector<fastjet::PseudoJet> &in,
                     const std::vector<float> &dmerge,
                     const std::vector<float> &dmerge_max) {
  JetClustering::FCCAnalysesJet result = initialise_FCCAnalysesJet();
  for (const auto &pjet : in) {
    result.jets.push_back(pjet);

    std::vector<fastjet::PseudoJet> consts = pjet.constituents();
    std::vector<int> tmpvec;
    for (const auto &constituent : consts) {
      tmpvec.push_back(constituent.user_index());
    }
    result.constituents.push_back(tmpvec);
  }
  result.exclusive_dmerge = dmerge;
  result.exclusive_dmerge_max = dmerge_max;
  return result;
}

std::vector<fastjet::PseudoJet>
build_jets(fastjet::ClusterSequence &cs, int exclusive, float cut, int sorted) {
  std::vector<fastjet::PseudoJet> pjets;

  if (sorted == 0) {
    if (exclusive == 0)
      pjets = fastjet::sorted_by_pt(cs.inclusive_jets(cut));
    else if (exclusive == 1)
      pjets = fastjet::sorted_by_pt(cs.exclusive_jets(cut));
    else if (exclusive == 2)
      pjets = fastjet::sorted_by_pt(cs.exclusive_jets(int(cut)));
    else if (exclusive == 3)
      pjets = fastjet::sorted_by_pt(cs.exclusive_jets_up_to(int(cut)));
    else if (exclusive == 4)
      pjets = fastjet::sorted_by_pt(cs.exclusive_jets_ycut(cut));
  } else if (sorted == 1) {
    if (exclusive == 0)
      pjets = fastjet::sorted_by_E(cs.inclusive_jets(cut));
    else if (exclusive == 1)
      pjets = fastjet::sorted_by_E(cs.exclusive_jets(cut));
    else if (exclusive == 2)
      pjets = fastjet::sorted_by_E(cs.exclusive_jets(int(cut)));
    else if (exclusive == 3)
      pjets = fastjet::sorted_by_E(cs.exclusive_jets_up_to(int(cut)));
    else if (exclusive == 4)
      pjets = fastjet::sorted_by_E(cs.exclusive_jets_ycut(cut));
  }
  return pjets;
}

std::vector<float> exclusive_dmerge(fastjet::ClusterSequence &cs,
                                    int do_dmarge_max) {
  const int Nmax = Nmax_dmerge;
  std::vector<float> result;
  for (int i = 1; i <= Nmax; i++) {
    float d;
    const int j = i;
    if (do_dmarge_max == 0)
      d = cs.exclusive_dmerge(j);
    else
      d = cs.exclusive_dmerge_max(j);
    result.push_back(d);
  }
  return result;
}

bool check(unsigned int n, int exclusive, float cut) {
  if (exclusive > 0 && n <= int(cut))
    return false;
  return true;
}

fastjet::RecombinationScheme recomb_scheme(int recombination) {
  fastjet::RecombinationScheme recomb_scheme;

  if (recombination == 0)
    recomb_scheme = fastjet::RecombinationScheme::E_scheme;
  else if (recombination == 1)
    recomb_scheme = fastjet::RecombinationScheme::pt_scheme;
  else if (recombination == 2)
    recomb_scheme = fastjet::RecombinationScheme::pt2_scheme;
  else if (recombination == 3)
    recomb_scheme = fastjet::RecombinationScheme::Et_scheme;
  else if (recombination == 4)
    recomb_scheme = fastjet::RecombinationScheme::Et2_scheme;
  else if (recombination == 5)
    recomb_scheme = fastjet::RecombinationScheme::BIpt_scheme;
  else if (recombination == 6)
    recomb_scheme = fastjet::RecombinationScheme::BIpt2_scheme;
  else
    recomb_scheme = fastjet::RecombinationScheme::external_scheme;

  return recomb_scheme;
}

resonanceBuilder::resonanceBuilder(float arg_resonance_mass) {
  m_resonance_mass = arg_resonance_mass;
}
ROOT::VecOps::RVec<fastjet::PseudoJet>
resonanceBuilder::operator()(ROOT::VecOps::RVec<fastjet::PseudoJet> legs) {
  ROOT::VecOps::RVec<fastjet::PseudoJet> result;
  int n = legs.size();
  if (n > 1) {
    ROOT::VecOps::RVec<bool> v(n);
    std::fill(v.end() - 2, v.end(), true);
    do {
      TLorentzVector reso;
      TLorentzVector reso_lv;
      for (int i = 0; i < n; ++i) {
        if (v[i]) {
          TLorentzVector leg_lv;
          leg_lv.SetXYZM(legs[i].px(), legs[i].py(), legs[i].pz(), legs[i].m());
          reso_lv += leg_lv;
        }
      }
      result.emplace_back(reso_lv);
    } while (std::next_permutation(v.begin(), v.end()));
  }
  if (result.size() > 1) {
    auto resonancesort = [&](fastjet::PseudoJet i, fastjet::PseudoJet j) {
      return (abs(m_resonance_mass - i.m()) < abs(m_resonance_mass - j.m()));
    };
    std::sort(result.begin(), result.end(), resonancesort);
    ROOT::VecOps::RVec<fastjet::PseudoJet>::const_iterator first =
        result.begin();
    ROOT::VecOps::RVec<fastjet::PseudoJet>::const_iterator last =
        result.begin() + 1;
    ROOT::VecOps::RVec<fastjet::PseudoJet> onlyBestReso(first, last);
    return onlyBestReso;
  } else {
    return result;
  }
}
recoilBuilder::recoilBuilder(float arg_sqrts) : m_sqrts(arg_sqrts){};
double recoilBuilder::operator()(ROOT::VecOps::RVec<fastjet::PseudoJet> in) {
  double result;
  auto recoil_p4 = TLorentzVector(0, 0, 0, m_sqrts);
  for (auto &v1 : in) {
    TLorentzVector tv1;
    tv1.SetPxPyPzE(v1.px(), v1.py(), v1.pz(), v1.e());
    recoil_p4 -= tv1;
  }
  result = recoil_p4.M();
  return result;
}

/// from Edi ///
  
std::vector<std::vector<int>> int_2d(){
  std::vector<std::vector<int>> result;
  for (int i=0; i<5; ++i){
    std::vector<int> result_tmp;
    for (int j=0; j<10; ++j){
      result_tmp.push_back(j);
    }
    result.push_back(result_tmp);
  }
  return result;
}

std::vector<std::vector<float>> float_2d(std::vector<std::vector<int>> constituents){
  std::vector<std::vector<float>> result;
  std::cout<<"-------EVENT---------"<<std::endl;
  for (int i=0; i<constituents.size(); ++i){
    std::cout<<"inner loop..."<<std::endl;
    std::vector<float> result_tmp;
    for (int j=0; j<constituents.at(i).size(); ++j){
      std::cout<<float(j)<<std::endl;
      result_tmp.push_back(float(j));
    }
    result.push_back(result_tmp);
  }
  return result;
}


std::vector<std::vector<float>> reshape2jet(ROOT::VecOps::RVec<float> var, std::vector<std::vector<int>> constituents){
  //for (auto& stuff : var) std::cout<<stuff<<std::endl;
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
  float dphi=0; 
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

std::vector<std::vector<float>> get_dR(std::vector<std::vector<float>> constituents_dTheta, std::vector<std::vector<float>> constituents_dPhi){
  std::vector<std::vector<float>> result;
  //for (auto& ind : indices){
  for (int j=0; j<constituents_dTheta.size(); j++){
    std::vector<float> tmp_res;
    for (int k=0; k<constituents_dTheta[j].size(); k++){
        tmp_res.push_back(std::sqrt(std::pow(constituents_dTheta[j][k], 2)+std::pow(constituents_dPhi[j][k], 2)));
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

//the same as above but I want the naming to not be confusing
std::vector<std::vector<float>> get_eRel(ROOT::VecOps::RVec<float> jet_e, std::vector<std::vector<float>> constituents_e){
  std::vector<std::vector<float>> result;
  //for (auto& ind : indices){
  for (int j=0; j<jet_e.size(); j++){
    std::vector<float> tmp_res;
    for (auto& consti_e : constituents_e[j]){
      tmp_res.push_back(consti_e/jet_e[j]);
    }
    result.push_back(tmp_res);
  }
  return result;
}

std::vector<std::vector<float>> get_dAngle(ROOT::VecOps::RVec<float> jet_px, 
                                           ROOT::VecOps::RVec<float> jet_py, 
                                           ROOT::VecOps::RVec<float> jet_pz, 
                                           std::vector<std::vector<float>> constituents_px,
                                           std::vector<std::vector<float>> constituents_py,
                                           std::vector<std::vector<float>> constituents_pz
                                           ){
    std::vector<std::vector<float>> result;
    for(int j=0; j<constituents_px.size(); ++j){ 
        std::vector<float> tmp_res; 
        for(int k=0; k<constituents_px[j].size(); ++k){  
            float dot = constituents_px[j][k] * jet_px[j]
                        + constituents_py[j][k] * jet_py[j]
                        + constituents_pz[j][k] * jet_pz[j];
            float lenSq1 = constituents_px[j][k] * constituents_px[j][k]
                           + constituents_py[j][k] * constituents_py[j][k]
                           + constituents_pz[j][k] * constituents_pz[j][k];
            float lenSq2 = jet_px[j] * jet_px[j]
                           + jet_py[j] * jet_py[j]
                           + jet_pz[j] * jet_pz[j];
            float norm = std::sqrt(lenSq1*lenSq2);
            float angle = std::acos(dot/norm);
            tmp_res.push_back(angle);
        }
        result.push_back(tmp_res);
    }
      return result;
}

std::vector<float> get_angularity(float kappa, float beta, ROOT::VecOps::RVec<float> jet_e, std::vector<std::vector<float>> constituents_e, std::vector<std::vector<float>> constituents_dAngle){
    
    std::vector<float> result;
    for(int j=0; j<constituents_e.size(); ++j){ 
        float tmp_res=0; 
        for(int k=0; k<constituents_e[j].size(); ++k){ 
            tmp_res+=std::pow((constituents_e[j][k]/jet_e[j]), kappa)*std::pow((constituents_dAngle[j][k]/TMath::Pi()), beta);
        }
        result.push_back(tmp_res);
    }
      return result;
}

  /// from Edi ///

} // namespace JetClusteringUtils

} // namespace FCCAnalyses
