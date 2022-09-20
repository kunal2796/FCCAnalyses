#include "FCCAnalyses/ReconstructedParticle.h"
#include <iostream>

namespace FCCAnalyses{

namespace ReconstructedParticle{

sel_pt::sel_pt(float arg_min_pt) : m_min_pt(arg_min_pt) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sel_pt::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if (std::sqrt(std::pow(p.momentum.x,2) + std::pow(p.momentum.y,2)) > m_min_pt) {
      result.emplace_back(p);
    }
  }
  return result;
}

sel_p::sel_p(float arg_min_p, float arg_max_p) : m_min_p(arg_min_p), m_max_p(arg_max_p)  {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sel_p::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    float momentum = std::sqrt(   std::pow(p.momentum.x,2)
                                + std::pow(p.momentum.y,2)
                                + std::pow(p.momentum.z,2) );
    if ( momentum > m_min_p && momentum < m_max_p ) {
      result.emplace_back(p);
    }
  }
  return result;
}

sel_charge::sel_charge(int arg_charge, bool arg_abs){m_charge = arg_charge; m_abs = arg_abs;};

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sel_charge::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if ((m_abs && abs(in[i].charge)==m_charge) || (m_charge==in[i].charge) ) {
      result.emplace_back(p);
    }
  }
  return result;
}



resonanceBuilder::resonanceBuilder(float arg_resonance_mass) {m_resonance_mass = arg_resonance_mass;}
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> resonanceBuilder::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  int n = legs.size();
  if (n >1) {
    ROOT::VecOps::RVec<bool> v(n);
    std::fill(v.end() - 2, v.end(), true);
    do {
      edm4hep::ReconstructedParticleData reso;
      TLorentzVector reso_lv;
      for (int i = 0; i < n; ++i) {
          if (v[i]) {
            reso.charge += legs[i].charge;
            TLorentzVector leg_lv;
            leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
            reso_lv += leg_lv;
          }
      }
      reso.momentum.x = reso_lv.Px();
      reso.momentum.y = reso_lv.Py();
      reso.momentum.z = reso_lv.Pz();
      reso.mass = reso_lv.M();
      result.emplace_back(reso);
    } while (std::next_permutation(v.begin(), v.end()));
  }
  if (result.size() > 1) {
    auto resonancesort = [&] (edm4hep::ReconstructedParticleData i ,edm4hep::ReconstructedParticleData j) { return (abs( m_resonance_mass -i.mass)<abs(m_resonance_mass-j.mass)); };
    std::sort(result.begin(), result.end(), resonancesort);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator first = result.begin();
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator last = result.begin() + 1;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> onlyBestReso(first, last);
    return onlyBestReso;
  } else {
    return result;
  }
}


recoilBuilder::recoilBuilder(float arg_sqrts) : m_sqrts(arg_sqrts) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  recoilBuilder::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  auto recoil_p4 = TLorentzVector(0, 0, 0, m_sqrts);
  for (auto & v1: in) {
    TLorentzVector tv1;
    tv1.SetXYZM(v1.momentum.x, v1.momentum.y, v1.momentum.z, v1.mass);
    recoil_p4 -= tv1;
  }
  auto recoil_fcc = edm4hep::ReconstructedParticleData();
  recoil_fcc.momentum.x = recoil_p4.Px();
  recoil_fcc.momentum.y = recoil_p4.Py();
  recoil_fcc.momentum.z = recoil_p4.Pz();
  recoil_fcc.mass = recoil_p4.M();
  result.push_back(recoil_fcc);
  return result;
};


sel_axis::sel_axis(bool arg_pos): m_pos(arg_pos) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> sel_axis::operator()(ROOT::VecOps::RVec<float> angle, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  for (size_t i = 0; i < angle.size(); ++i) {
    if (m_pos==1 && angle.at(i)>0.) result.push_back(in.at(i));
    if (m_pos==0 && angle.at(i)<0.) result.push_back(in.at(i));;
  }
  return result;
}


sel_tag::sel_tag(bool arg_pass): m_pass(arg_pass) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> sel_tag::operator()(ROOT::VecOps::RVec<bool> tags, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
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



// Angular separation between the particles of a collection:
//   arg_delta = 0 / 1 / 2 :   return delta_max, delta_min, delta_average

angular_separationBuilder::angular_separationBuilder( int  arg_delta) : m_delta(arg_delta) {};
float angular_separationBuilder::operator() ( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {

 float result = -9999;

 float dmax = -999;
 float dmin = 999;
 float sum = 0;
 float npairs = 0;
 for (int i=0; i < in.size(); i++) {
  if ( in.at(i).energy < 0) continue;    // "dummy" particle - cf selRP_matched_to_list
  TVector3 p1( in.at(i).momentum.x, in.at(i).momentum.y, in.at(i).momentum.z );
  for (int j=i+1; j < in.size(); j++) {
    if ( in.at(j).energy < 0) continue;   // "dummy" particle
    TVector3 p2( in.at(j).momentum.x, in.at(j).momentum.y, in.at(j).momentum.z );
    float delta_ij = fabs( p1.Angle( p2 ) );
    if ( delta_ij > dmax) dmax = delta_ij;
    if ( delta_ij < dmin) dmin = delta_ij;
    sum = sum + delta_ij;
    npairs ++;
  }
 }
 float delta_max = dmax;
 float delta_min = dmin;
 float delta_ave = sum / npairs;

 if (m_delta == 0 ) result = delta_max;
 if (m_delta == 1 ) result = delta_min;
 if (m_delta == 2 ) result = delta_ave;

 return result;
}


ROOT::VecOps::RVec<float> get_pt(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
 ROOT::VecOps::RVec<float> result;
 for (size_t i = 0; i < in.size(); ++i) {
   result.push_back(sqrt(in[i].momentum.x * in[i].momentum.x + in[i].momentum.y * in[i].momentum.y));
 }
 return result;
}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> merge(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> y) {
  //to be keept as ROOT::VecOps::RVec
  std::vector<edm4hep::ReconstructedParticleData> result;
  result.reserve(x.size() + y.size());
  result.insert( result.end(), x.begin(), x.end() );
  result.insert( result.end(), y.begin(), y.end() );
  return ROOT::VecOps::RVec(result);
}


ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> remove(
  		ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x,
  		ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> y) {
  //to be kept as ROOT::VecOps::RVec
  std::vector<edm4hep::ReconstructedParticleData> result;
  result.reserve( x.size() );
  result.insert( result.end(), x.begin(), x.end() );
  float epsilon = 1e-8;
  for (size_t i = 0; i < y.size(); ++i) {
    float mass1 = y.at(i).mass;
    float px1 = y.at(i).momentum.x;
    float py1 = y.at(i).momentum.y;
    float pz1 = y.at(i).momentum.z;
    for(std::vector<edm4hep::ReconstructedParticleData>::iterator
          it = std::begin(result); it != std::end(result); ++it) {
      float mass2 = it->mass;
      float px2 = it->momentum.x;
      float py2 = it->momentum.y;
      float pz2 = it->momentum.z;
      if ( abs(mass1-mass2) < epsilon &&
	   abs(px1-px2) < epsilon &&
	   abs(py1-py2) < epsilon &&
	   abs(pz1-pz2) < epsilon ) {
        result.erase(it);
        break;
      }
    }
  }
  return ROOT::VecOps::RVec(result);
}




ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> get(ROOT::VecOps::RVec<int> index, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  for (size_t i = 0; i < index.size(); ++i) {
    if (index[i]>-1)
      result.push_back(in.at(index[i]));
    //else
    //  std::cout << "electron index negative " << index[i]<<std::endl;
  }
  return result;
}


ROOT::VecOps::RVec<float> get_mass(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.mass);
  }
  return result;
}

ROOT::VecOps::RVec<float> get_eta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Eta());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_phi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Phi());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_e(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.energy);
  }
  return result;
}

ROOT::VecOps::RVec<float> get_p(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.P());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_px(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.x);
  }
  return result;
}


ROOT::VecOps::RVec<float> get_py(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.y);
  }
  return result;
}

ROOT::VecOps::RVec<float> get_pz(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.momentum.z);
  }
  return result;
}

ROOT::VecOps::RVec<float> get_charge(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    result.push_back(p.charge);
  }
  return result;
}

ROOT::VecOps::RVec<float> get_y(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Rapidity());
  }
  return result;
}

ROOT::VecOps::RVec<float> get_theta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<float> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv.Theta());
  }
  return result;
}

ROOT::VecOps::RVec<TLorentzVector> get_tlv(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<TLorentzVector> result;
  for (auto & p: in) {
    TLorentzVector tlv;
    tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
    result.push_back(tlv);
  }
  return result;
}

TLorentzVector get_tlv(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, int index) {
  TLorentzVector result;
  auto & p = in[index];
  result.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
  return result;
}

TLorentzVector get_tlv(edm4hep::ReconstructedParticleData in) {
  TLorentzVector result;
  result.SetXYZM(in.momentum.x, in.momentum.y, in.momentum.z, in.mass);
  return result;
}

ROOT::VecOps::RVec<int>
get_type(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
  ROOT::VecOps::RVec<int> result;
  for (auto & p: in) {
    result.push_back(p.type);
  }
  return result;
}


int get_n(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> x) {
  int result =  x.size();
  return result;
}


ROOT::VecOps::RVec<bool> getJet_btag(ROOT::VecOps::RVec<int> index, ROOT::VecOps::RVec<edm4hep::ParticleIDData> pid, ROOT::VecOps::RVec<float> values){
  ROOT::VecOps::RVec<bool> result;
  //std::cout << "========================new event=======================" <<std::endl;
  for (size_t i = 0; i < index.size(); ++i) {
    result.push_back(values.at(pid.at(index.at(i)).parameters_begin));

    //std::cout << pid.at(index.at(i)).parameters_begin << "  ==  " << pid.at(index.at(i)).parameters_end << std::endl;
    //for (unsigned j = pid.at(index.at(i)).parameters_begin; j != pid.at(index.at(i)).parameters_end; ++j) {
    //  std::cout << " values : " << values.at(j) << std::endl;
    //}
  }
  return result;
}

int getJet_ntags(ROOT::VecOps::RVec<bool> in) {
  int result =  0;
  for (size_t i = 0; i < in.size(); ++i)
    if (in.at(i))result+=1;
  return result;
}


/// ------ ///
// From Edi //

sel_template::sel_template(float arg_pass){m_pass = arg_pass;}

template<class G>
ROOT::VecOps::RVec<G> sel_template::operator()(ROOT::VecOps::RVec<float> tags, ROOT::VecOps::RVec<G> in){
//ROOT::VecOps::RVec<G> sel_template::operator()(ROOT::VecOps::RVec<float> tags, ROOT::VecOps::RVec<G> in){
//////////ROOT::VecOps::RVec<int> sel_template::operator()(ROOT::VecOps::RVec<float> tags, ROOT::VecOps::RVec<int> in){
//ROOT::VecOps::RVec<G> sel_template(ROOT::VecOps::RVec<G> in){
  ////////ROOT::VecOps::RVec<G> result;
  ROOT::VecOps::RVec<int> result;
  //ROOT::VecOps::RVec<G> result = in;
  bool m_pass=true;
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
std::vector<std::vector<int>> sel_template::operator()(std::vector<std::vector<float>> tag_vector, std::vector<std::vector<int>> in){
//ROOT::VecOps::RVec<G> sel_template::operator()(ROOT::VecOps::RVec<float> tags, ROOT::VecOps::RVec<G> in){
//////////ROOT::VecOps::RVec<int> sel_template::operator()(ROOT::VecOps::RVec<float> tags, ROOT::VecOps::RVec<int> in){
//ROOT::VecOps::RVec<G> sel_template(ROOT::VecOps::RVec<G> in){
  ////////ROOT::VecOps::RVec<G> result;
  std::vector<std::vector<int>> result;
  //ROOT::VecOps::RVec<G> result = in;
  //////bool m_pass=true;
  for (size_t i = 0; i < in.size(); ++i) {
    std::vector<int> tmp_res;
    for (size_t j = 0; j < in[i].size(); ++j) {
      if (m_pass) {
        if (tag_vector.at(i).at(j)) tmp_res.push_back(in.at(i).at(j));
      }
      else {
        if (!tag_vector.at(i).at(j)) tmp_res.push_back(in.at(i).at(j));
      }
    }
    result.push_back(tmp_res);
  }
  return result;
}

std::vector<std::vector<float>> sel_template::operator()(std::vector<std::vector<float>> tag_vector, std::vector<std::vector<float>> in){
  std::vector<std::vector<float>> result;
  //ROOT::VecOps::RVec<G> result = in;
  //////bool m_pass=true;
  for (size_t i = 0; i < in.size(); ++i) {
    std::vector<float> tmp_res;
    for (size_t j = 0; j < in[i].size(); ++j) {
      if (m_pass) {
        if (tag_vector.at(i).at(j)) tmp_res.push_back(in.at(i).at(j));
      }
      else {
        if (!tag_vector.at(i).at(j)) tmp_res.push_back(in.at(i).at(j));
      }
    }
    result.push_back(tmp_res);
  }
  return result;
}

//template ROOT::VecOps::RVec<float> sel_template<float>(ROOT::VecOps::RVec<float>);
///template ROOT::VecOps::RVec<float> sel_template<float>(ROOT::VecOps::RVec<bool>, ROOT::VecOps::RVec<float>);
///template ROOT::VecOps::RVec<int> sel_template<int>(ROOT::VecOps::RVec<bool>, ROOT::VecOps::RVec<int>);
///template ROOT::VecOps::RVec<double> sel_template<double>(ROOT::VecOps::RVec<bool>, ROOT::VecOps::RVec<double>);
////////////////////template ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> sel_template::operator()(ROOT::VecOps::RVec<bool>, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>);
//template ROOT::VecOps::RVec<edm4hep::MCParticleData> sel_template<edm4hep::MCParticleData>(ROOT::VecOps::RVec<bool>, ROOT::VecOps::RVec<edm4hep::MCParticleData>);
//template ROOT::VecOps::RVec<edm4hep::TrackState> sel_template<edm4hep::TrackState>(ROOT::VecOps::RVec<bool>, ROOT::VecOps::RVec<edm4hep::TrackState>);
//template ROOT::VecOps::RVec<JetClusteringUtils::FCCAnalysesJet> sel_template<JetClusteringUtils::FCCAnalysesJet>(ROOT::VecOps::RVec<bool>, ROOT::VecOps::RVec<JetClusteringUtils::FCCAnalysesJet>);
//think about adding additional instantiations as the need arises...

ROOT::VecOps::RVec<int> index_splitter(ROOT::VecOps::RVec<int> ind){
//ROOT::VecOps::RVec<int> index_splitter(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<int> ind){
  ROOT::VecOps::RVec<int> result;
  int charged_counter, neutral_counter=0;
  for(size_t i = 0; i < ind.size(); ++i){
    if(ind[i]==1){
      result.push_back(charged_counter);
      charged_counter+=1;
    }
    else if(ind[i]==0){
      result.push_back(neutral_counter); 
      neutral_counter+=1;
    }
  }
  return result;
}

std::vector<std::vector<int>> index_converter(std::vector<std::vector<int>> RP_ind, ROOT::VecOps::RVec<int> split_ind){
  std::vector<std::vector<int>> result;
  //for(size_t i = 0; i < in.size(); ++i){
  for(auto& indices : RP_ind){
    std::vector<int> tmp_res;
    for(auto& i : indices){
      tmp_res.push_back(split_ind.at(i));
    }
    result.push_back(tmp_res);
  }
  return result;
}

ROOT::VecOps::RVec<float> is_particle(ROOT::VecOps::RVec<int> index, const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> & in){
  ROOT::VecOps::RVec<float> result(in.size(), 0);
  for(auto& i : index){
    result[i] = 1;
  }
  return result;
}

std::vector<std::vector<float>> one_hot_encode(ROOT::VecOps::RVec<float> flavour){
  std::vector<std::vector<float>> result;
  int min = -5; //std::min_element(flavour.begin(),flavour.end());
  int max = 5; //std::max_element(flavour.begin(),flavour.end());
  for(int i = min; i<=max; ++i){
    std::vector<float> zeros(flavour.size(),0);
    result.push_back(zeros);
  }
  for(size_t k = 0; k<flavour.size(); ++k){
    for(int j = min; j<=max; ++j){
      if(j==flavour[k]){
        result.at(j+5).at(k)=1;
      }
    }
  }
  return result;
}

ROOT::VecOps::RVec<float> get_PID(const ROOT::VecOps::RVec<int> & recind, 
				  const ROOT::VecOps::RVec<int> & mcind, 
				  const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> & reco,
				  const ROOT::VecOps::RVec<edm4hep::MCParticleData> & mc){
  ROOT::VecOps::RVec<float> result;
  result = ReconstructedParticle2MC::getRP2MC_pdg(recind, mcind, reco, mc);
  return result;
}

ROOT::VecOps::RVec<int> is_S(ROOT::VecOps::RVec<float> PID){
  
  ROOT::VecOps::RVec<int> result;
  for(auto & id : PID){
    //the below conditions checks on strange mesons (and K0L separately), and strange baryons
    if((id==130)||(std::abs(int((id/100))%10)==3)||(std::abs(int((id/1000))%10)==3)){
      result.push_back(1);
    }
    else result.push_back(0);
  }
  return result;
  
  
}
  
ROOT::VecOps::RVec<int> is_Kaon(ROOT::VecOps::RVec<float> PID){
  
  ROOT::VecOps::RVec<int> result;
  for(auto & id : PID){
    //the below conditions checks only on Kaons
    if((id==130)||(std::abs(id)==321)||(std::abs(id)==311)||(std::abs(id)==310)){
      result.push_back(1);
    }
    else result.push_back(0);
    }
  return result;
  
  
}
  

ROOT::VecOps::RVec<int> is_Kaon_smearedUniform010(ROOT::VecOps::RVec<float> PID){
  
  ROOT::VecOps::RVec<int> result;
  int tmp_res;
  for(auto & id : PID){
    //the below conditions checks only on Kaons
    if((id==130)||(std::abs(id)==321)||(std::abs(id)==311)||(std::abs(id)==310)){
      tmp_res=1;
      }
    else tmp_res=0;
    if(gRandom->Uniform()<=0.1) tmp_res = (tmp_res+1)%2;
    result.push_back(tmp_res);
  } 
  return result;
}
  

ROOT::VecOps::RVec<int> is_Kaon_smearedUniform005(ROOT::VecOps::RVec<float> PID){

  ROOT::VecOps::RVec<int> result;
  int tmp_res;
  for(auto & id : PID){
    //the below conditions checks only on Kaons
    if((id==130)||(std::abs(id)==321)||(std::abs(id)==311)||(std::abs(id)==310)){
      tmp_res=1;
    }
    else tmp_res=0;
      if(gRandom->Uniform()<=0.05) tmp_res = (tmp_res+1)%2;
      result.push_back(tmp_res);
  } 
  return result;
}
  

ROOT::VecOps::RVec<int> is_Kaon_smearedUniform001(ROOT::VecOps::RVec<float> PID){
  
  ROOT::VecOps::RVec<int> result;
  int tmp_res;
  for(auto & id : PID){
      //the below conditions checks only on Kaons
    if((id==130)||(std::abs(id)==321)||(std::abs(id)==311)||(std::abs(id)==310)){
      tmp_res=1;
    }
    else tmp_res=0;
    if(gRandom->Uniform()<=0.01) tmp_res = (tmp_res+1)%2;
    result.push_back(tmp_res);
  } 
  return result;
}
  
  
/// ------ ///

}//end NS ReconstructedParticle

}//end NS FCCAnalyses
