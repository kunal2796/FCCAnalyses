#include "VertexFinderLCFIPlus.h"
#include <iostream>

#include "TFile.h"
#include "TString.h"

using namespace VertexFinderLCFIPlus;

bool debug_me = false;


VertexingUtils::FCCAnalysesSV VertexFinderLCFIPlus::get_SV_jets(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
								ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
								VertexingUtils::FCCAnalysesVertex PV,
								ROOT::VecOps::RVec<bool> isInPrimary,
								ROOT::VecOps::RVec<fastjet::PseudoJet> jets,
								std::vector<std::vector<int>> jet_consti,
								bool V0_rej,
								double chi2_cut, double invM_cut, double chi2Tr_cut) {

  // find SVs using LCFI+ (clustering first)
  // change to vec of vec (RVec of RVec breaking) to separate SV from diff jets, currently don't separate SVs by jet
  // added a vector containing number of SVs per jet
  
  VertexingUtils::FCCAnalysesSV SV;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;
  ROOT::VecOps::RVec<int> nSV_jet;
  SV.vtx = result;
  SV.nSV_jet = nSV_jet;

  int nJet = jets.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks;

  // retrieve tracks from reco particles & get a vector with their indices in the reco collection
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks   = ReconstructedParticle2Track::getRP2TRK( recoparticles, thetracks );
  ROOT::VecOps::RVec<int> reco_ind_tracks     = ReconstructedParticle2Track::get_recoindTRK( recoparticles, thetracks );
  if(tracks.size() != reco_ind_tracks.size()) std::cout<<"ERROR: reco index vector not the same size as no of tracks"<<std::endl;

  if(tracks.size() != isInPrimary.size()) std::cout<<"ERROR: isInPrimary vector size not the same as no of tracks"<<std::endl;

  if(debug_me) std::cout<<"tracks extracted from the reco particles"<<std::endl;

  // find SVs inside jet loop
  //
  for (unsigned int j=0; j<nJet; j++) {

    // remove primary tracks
    // separate non-primary tracks by jet
    std::vector<int> i_jetconsti = jet_consti[j];
    for (int ctr=0; ctr<tracks.size(); ctr++) {
      if(isInPrimary[ctr]) continue; // remove primary tracks
      if(std::find(i_jetconsti.begin(), i_jetconsti.end(), reco_ind_tracks[ctr]) == i_jetconsti.end()) {
	np_tracks.push_back(tracks[ctr]); // separate tracks by jet
      }
    }
    
    if(debug_me) std::cout<<"primary tracks removed; there are "<<np_tracks.size()<<" non-primary tracks in jet#"<<j+1<<std::endl;
    
    // V0 rejection (tight)
    // perform V0 rejection with tight constraints if user chooses
    ROOT::VecOps::RVec<edm4hep::TrackState> tracks_fin = V0rejection_tight(np_tracks, PV, V0_rej);
    
    if(debug_me) {
      std::cout<<np_tracks.size()-tracks_fin.size()<<" V0 tracks removed"<<std::endl;
      std::cout<<"now starting to find secondary vertices..."<<std::endl;
    }
    
    // start finding SVs
    ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> i_result = findSVfromTracks(tracks_fin, PV, chi2_cut, invM_cut, chi2Tr_cut);

    int i_nSV = i_result.size();
    nSV_jet.push_back(i_nSV);

    //result.insert(result.end(), i_result.begin(), i_result.end()); // compilation error
    for(VertexingUtils::FCCAnalysesVertex i_sv : i_result) result.push_back(i_sv);
    
    // clean-up
    i_result.clear();
    np_tracks.clear();
    tracks_fin.clear();
  }

  if(debug_me) std::cout<<"no more SVs can be reconstructed"<<std::endl;
  
  // nSV_jet gives a handle on deciding which SVs are in which jet
  SV.vtx = result;
  SV.nSV_jet = nSV_jet;
  //
  return SV;
}

//ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> VertexFinderLCFIPlus::get_SV_event(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
VertexingUtils::FCCAnalysesSV VertexFinderLCFIPlus::get_SV_event(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
								 ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
								 VertexingUtils::FCCAnalysesVertex PV,
								 ROOT::VecOps::RVec<bool> isInPrimary,
								 bool V0_rej,
								 double chi2_cut, double invM_cut, double chi2Tr_cut) {
  
  
  if(debug_me) std::cout << "Starting SV finding!" << std::endl;

  // find SVs using LCFI+ (w/o clustering)
  // still need to think a little about jet clustering using SVs & pseudo-vertices as seeds
  
  VertexingUtils::FCCAnalysesSV SV;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;
  SV.vtx = result;

  // retrieve the tracks associated to the recoparticles
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks = ReconstructedParticle2Track::getRP2TRK( recoparticles, thetracks );

  if(debug_me) std::cout<<"tracks extracted from the reco particles"<<std::endl;

  if(tracks.size() != isInPrimary.size()) std::cout<<"ISSUE: track vector and primary-nonprimary vector of diff sizes"<<std::endl;

  // remove primary tracks
  ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks;
  for(unsigned int i=0; i<tracks.size(); i++) {
    if (!isInPrimary[i]) np_tracks.push_back(tracks[i]);
  }

  if(debug_me) std::cout<<"primary tracks removed; there are "<<np_tracks.size()<<" non-primary tracks in the event"<<std::endl;

  // V0 rejection (tight)
  // perform V0 rejection with tight constraints if user chooses
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks_fin = V0rejection_tight(np_tracks, PV, V0_rej);

  if(debug_me) {
    std::cout<<np_tracks.size()-tracks_fin.size()<<" V0 tracks removed"<<std::endl;
    std::cout<<"now starting to find secondary vertices..."<<std::endl;
  }
  
  //if(debug_me) std::cout << "tracks_fin.size() = " << tracks_fin.size() << std::endl;

  // start finding SVs (only if there are 2 or more tracks)
  result = findSVfromTracks(tracks_fin, PV, chi2_cut, invM_cut, chi2Tr_cut);

  //if(debug_me) std::cout<<"no more SVs can be reconstructed"<<std::endl;
  
  SV.vtx = result;
  //
  return SV;
  //return result;
}

//ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> VertexFinderLCFIPlus::get_SV_event(ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
VertexingUtils::FCCAnalysesSV VertexFinderLCFIPlus::get_SV_event(ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
								 VertexingUtils::FCCAnalysesVertex PV,
								 bool V0_rej,
								 double chi2_cut, double invM_cut, double chi2Tr_cut) {
  
  // find SVs from non-primary tracks using LCFI+ (w/o clustering)
  // still need to think a little about jet clustering using SVs & pseudo-vertices as seeds
  // primary - non-primary separation done externally
  
  VertexingUtils::FCCAnalysesSV SV;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;
  SV.vtx = result;

  // V0 rejection (tight)
  // perform V0 rejection with tight constraints if user chooses
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks_fin = V0rejection_tight(np_tracks, PV, V0_rej);

  if(debug_me) {
    std::cout<<np_tracks.size()-tracks_fin.size()<<" V0 tracks removed"<<std::endl;
    std::cout<<"now starting to find secondary vertices..."<<std::endl;
  }

  // start finding SVs (only if there are 2 or more tracks)
  result = findSVfromTracks(tracks_fin, PV, chi2_cut, invM_cut, chi2Tr_cut);
  
  SV.vtx = result;
  //
  return SV;
  //return result;
}


ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> VertexFinderLCFIPlus::get_all_vertices(VertexingUtils::FCCAnalysesVertex PV,
											     VertexingUtils::FCCAnalysesSV SV) {
  // Returns a vector of all vertices (PV and SVs)
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;
  result.push_back(PV);
  for (auto &p:SV.vtx){
    result.push_back(p);
  }
  return result;  
}


//
ROOT::VecOps::RVec<int> VertexFinderLCFIPlus::VertexSeed_best(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							      VertexingUtils::FCCAnalysesVertex PV,
							      double chi2_cut, double invM_cut) {

  // gives indices of the best pair of tracks
  // maybe also update and write one to get first pair to pass constraints

  ROOT::VecOps::RVec<int> result;
  int isel = 0;
  int jsel = 1;
  
  int nTr = tracks.size();
  // push empty tracks to make a size=2 vector
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  tr_pair.push_back(tracks[0]);
  tr_pair.push_back(tracks[1]);
  VertexingUtils::FCCAnalysesVertex vtx_seed;
  double chi2_min = 99;
  
  for(unsigned int i=0; i<nTr-1; i++) {
    tr_pair[0] = tracks[i];
    
    for(unsigned int j=i+1; j<nTr; j++) {
      tr_pair[1] = tracks[j];
      
      // V0 rejection (loose)
      ROOT::VecOps::RVec<bool> isInV0 = isV0(tr_pair, PV, false);
      if(isInV0[0] && isInV0[1]) continue;
      
      vtx_seed = VertexFitterSimple::VertexFitter_Tk(2, tr_pair);
      
      // Constraints
      bool pass = check_constraints(vtx_seed, tr_pair, PV, true, chi2_cut, invM_cut);
      if(!pass) continue;
      
      // if a pair passes all constraints compare chi2, store lowest chi2
      double chi2_seed = vtx_seed.vertex.chi2; // normalised but nDOF=1 for nTr=2      
      if(chi2_seed < chi2_min) {
	isel = i; jsel =j;
	chi2_min = chi2_seed;
      }
    }
  }

  if(chi2_min != 99){
    result.push_back(isel); 
    result.push_back(jsel);
  }
  return result;
}

ROOT::VecOps::RVec<int> VertexFinderLCFIPlus::addTrack_best(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							    ROOT::VecOps::RVec<int> vtx_tr,
							    VertexingUtils::FCCAnalysesVertex PV,
							    double chi2_cut, double invM_cut, double chi2Tr_cut) {
  // adds index of the best track to the (seed) vtx
  
  ROOT::VecOps::RVec<int> result = vtx_tr;
  if(tracks.size() == vtx_tr.size()) return result;
  
  int isel = -1;

  int nTr = tracks.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_vtx;
  VertexingUtils::FCCAnalysesVertex vtx;
  double chi2_min = 99;

  // add tracks of the previously formed vtx to a vector
  for(int tr : vtx_tr) {
    if(debug_me) std::cout << "Track integer: " << tr << std::endl;
    if(debug_me) std::cout <<  "Track value: " << tracks[tr] << std::endl;
    tr_vtx.push_back(tracks[tr]);
  }
  int iTr = tr_vtx.size();
  // add an empty track to increase vector size by 1
  tr_vtx.push_back(tracks[0]);

  // find best track to add to the vtx
  for(unsigned int i=0; i<nTr; i++) {
    if(std::find(vtx_tr.begin(), vtx_tr.end(), i) != vtx_tr.end()) continue;
    tr_vtx[iTr] = tracks[i];
    
    vtx = VertexFitterSimple::VertexFitter_Tk(2, tr_vtx);

    // Constraints
    bool pass = check_constraints(vtx, tr_vtx, PV, false, chi2_cut, invM_cut, chi2Tr_cut);
    if(!pass) continue;
    
    // if a track passes all constraints compare chi2, store lowest chi2
    double chi2_vtx = vtx.vertex.chi2; // normalised
    double nDOF = 2*(iTr+1) - 3;       // nDOF = 2*nTr - 3
    chi2_vtx = chi2_vtx * nDOF;
    if(chi2_vtx < chi2_min) {
      isel = i;
      chi2_min = chi2_vtx;
    }    
  }

  if(isel>=0) result.push_back(isel);
  return result;
}

ROOT::VecOps::RVec<edm4hep::TrackState> VertexFinderLCFIPlus::V0rejection_tight(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
										VertexingUtils::FCCAnalysesVertex PV,
										bool V0_rej) {
  // perform V0 rejection with tight constraints if user chooses
  ROOT::VecOps::RVec<edm4hep::TrackState> result;
  if(V0_rej) {
    bool tight = true;
    ROOT::VecOps::RVec<bool> isInV0 = isV0(tracks, PV, tight);
    for(unsigned int i=0; i<isInV0.size(); i++) {
      if (!isInV0[i]) result.push_back(tracks[i]);
    }
  }
  else result = tracks;
  //
  return result;
}

ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> VertexFinderLCFIPlus::findSVfromTracks(ROOT::VecOps::RVec<edm4hep::TrackState> tracks_fin,
											     VertexingUtils::FCCAnalysesVertex PV,
											     double chi2_cut, double invM_cut, double chi2Tr_cut) {

  // find SVs (only if there are 2 or more tracks)
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;

  while(tracks_fin.size() > 1) {
    // find vertex seed
    ROOT::VecOps::RVec<int> vtx_seed = VertexSeed_best(tracks_fin, PV, chi2_cut, invM_cut);
    
    if(debug_me){
      std::cout << "tracks_fin.size(): " << tracks_fin.size() << std::endl;
      for(int i=0; i<vtx_seed.size();i++)
	std::cout << "vtx_seed: " << vtx_seed[i] << std::endl;
    }
    if(vtx_seed.size() == 0) break;
    
    // add tracks to the seed
    // check if a track is added; if not break loop
    ROOT::VecOps::RVec<int> vtx_fin = vtx_seed;
    int vtx_fin_size = 0; // to start the loop
    while(vtx_fin_size != vtx_fin.size()) {
      vtx_fin_size = vtx_fin.size();
      vtx_fin = addTrack_best(tracks_fin, vtx_fin, PV, chi2_cut, invM_cut, chi2Tr_cut);
    }
    
    // fit tracks to SV and remove from tracks_fin
    ROOT::VecOps::RVec<edm4hep::TrackState> tr_vtx_fin;
    for(int i_tr : vtx_fin){
      tr_vtx_fin.push_back(tracks_fin[i_tr]);
      if(debug_me) std::cout << "Pushing back tracks_fin[i_tr]" << std::endl;
    }
    VertexingUtils::FCCAnalysesVertex sec_vtx = VertexFitterSimple::VertexFitter_Tk(2, tr_vtx_fin); // flag 2 for SVs

    // see if we can also get indices in the reco collection (for tracks forming an SV)
    //sec_vtx.reco_ind = VertexFitterSimple::get_reco_ind(recoparticles,thetracks); // incorrect

    result.push_back(sec_vtx);
    //
    ROOT::VecOps::RVec<edm4hep::TrackState> temp = tracks_fin;
    tracks_fin.clear();
    for(unsigned int t=0; t<temp.size(); t++) {
      if(std::find(vtx_fin.begin(), vtx_fin.end(), t) == vtx_fin.end()) tracks_fin.push_back(temp[t]);
    }
    // all this cause don't know how to remove multiple elements at once

    if(debug_me) std::cout<<result.size()<<" SV found"<<std::endl;
  }

  //
  return result;
}

bool VertexFinderLCFIPlus::check_constraints(VertexingUtils::FCCAnalysesVertex vtx,
					     ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
					     VertexingUtils::FCCAnalysesVertex PV,
					     bool seed,
					     double chi2_cut, double invM_cut, double chi2Tr_cut) {
  // if all constraints pass -> true
  // if any constraint fails -> false

  bool result = true;

  int iTr = chi2_tr.size() - 1;  // final track index
  
  // Constraints
  // chi2 < cut (9)
  double chi2 = vtx.vertex.chi2; // normalised
  double nDOF = 2*(iTr+1) - 3;   // nDOF = 2*nTr - 3
  chi2 = chi2 * nDOF;
  if(chi2 >= chi2_cut) result = false;
  //
  // invM < cut (10GeV)
  double invM = VertexingUtils::get_invM(vtx);
  if(invM >= invM_cut) result = false;
  //
  // invM < sum of energy
  double E_tracks = 0.;
  for(edm4hep::TrackState tr_e : tracks) E_tracks += VertexingUtils::get_trackE(tr_e);
  if(invM >= E_tracks) result = false;
  //
  // momenta sum & vtx r on same side
  double angle = VertexingUtils::get_PV2vtx_angle(tracks, vtx, PV);
  if(angle<0) result = false;
  //
  if(!seed) {
    // chi2_contribution(track) < threshold
    ROOT::VecOps::RVec<float> chi2_tr = vtx.reco_chi2;
    if(chi2_tr[iTr] >= chi2Tr_cut) result = false;    // threshold = 5 ok?
  }
  //
  return result;
}

ROOT::VecOps::RVec<bool> VertexFinderLCFIPlus::isV0(ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
						    VertexingUtils::FCCAnalysesVertex PV,
						    bool tight) {
  // V0 rejection
  //
  // take all non-primary tracks & assign "true" to pairs that form V0
  // if(tight)  -> tight constraints
  // if(!tight) -> loose constraints

  int nTr = np_tracks.size();

  ROOT::VecOps::RVec<bool> result(nTr, false);
  // true -> forms a V0, false -> doesn't form a V0
  if(nTr<2) return result;
  
  edm4hep::Vector3f r_PV = PV.vertex.position; // in mm  
  
  ROOT::VecOps::RVec<edm4hep::TrackState> t_pair;
  // push empty tracks to make a size=2 vector
  edm4hep::TrackState tr_i;
  edm4hep::TrackState tr_j;
  t_pair.push_back(tr_i);
  t_pair.push_back(tr_j);
  VertexingUtils::FCCAnalysesVertex V0;
  //
  const double m_pi = 0.13957039; // pi+- mass [GeV]
  const double m_p  = 0.93827208; // p+- mass
  const double m_e  = 0.00051099; // e+- mass
  //
  for(unsigned int i=0; i<nTr-1; i++) {
    if(result[i] == true) continue;
    t_pair[0] = np_tracks[i];

    for(unsigned int j=i+1; j<nTr; j++) {
      if(result[j] == true) continue;
      t_pair[1] = np_tracks[j];

      V0 = VertexFitterSimple::VertexFitter_Tk(2, t_pair);

      // invariant masses for V0 candidates
      double invM_Ks      = VertexingUtils::get_invM_pairs(V0, m_pi, m_pi);
      double invM_Lambda1 = VertexingUtils::get_invM_pairs(V0, m_pi, m_p);
      double invM_Lambda2 = VertexingUtils::get_invM_pairs(V0, m_p, m_pi);
      double invM_Gamma   = VertexingUtils::get_invM_pairs(V0, m_e, m_e);

      // V0 candidate distance from PV
      edm4hep::Vector3f r_V0 = V0.vertex.position; // in mm
      // does Vector3f class has similar functions as root vectors?
      TVector3 r_V0_PV(r_V0[0] - r_PV[0], r_V0[1] - r_PV[1], r_V0[2] - r_PV[2]);
      double r = r_V0_PV.Mag(); // in mm

      // angle b/n V0 candidate momentum & PV-V0 displacement vector
      double p_r = VertexingUtils::get_PV2V0angle(V0, PV);

      if(tight) {
	// Ks
	if(invM_Ks>0.493 && invM_Ks<0.503 && r>0.5 && p_r>0.999) {
	  result[i] = true;
	  result[j] = true;
	  break;
	}

	// Lambda0
	else if(invM_Lambda1>1.111 && invM_Lambda1<1.121 && r>0.5 && p_r>0.99995) {
	  result[i] = true;
	  result[j] = true;
	  break;
	}
	else if(invM_Lambda2>1.111 && invM_Lambda2<1.121 && r>0.5 && p_r>0.99995) {
	  result[i] = true;
	  result[j] = true;
	  break;
	}

	// photon conversion
	else if(invM_Gamma<0.005 && r>9 && p_r>0.99995) {
	  result[i] = true;
	  result[j] = true;
	  break;
	}	
      }

      else {
	// Ks
	if(invM_Ks>0.488 && invM_Ks<0.508 && r>0.3 && p_r>0.999) {
	  result[i] = true;
	  result[j] = true;
	}
	
	// Lambda0
	else if(invM_Lambda1>1.106 && invM_Lambda1<1.126 && r>0.3 && p_r>0.999) {
	  result[i] = true;
	  result[j] = true;
	}
	else if(invM_Lambda2>1.106 && invM_Lambda2<1.126 && r>0.3 && p_r>0.999) {
	  result[i] = true;
	  result[j] = true;
	}
	
	// photon conversion
	else if(invM_Gamma<0.01 && r>9 && p_r>0.999) {
	  result[i] = true;
	  result[j] = true;
	}	
      }
      
    }
  }

  return result;
}


///////////////////////////
//** V0 Reconstruction **//
///////////////////////////

VertexingUtils::FCCAnalysesV0 VertexFinderLCFIPlus::get_V0s(ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
							    VertexingUtils::FCCAnalysesVertex PV,
							    bool tight,
							    double chi2_cut) {
  // V0 reconstruction
  // if(tight)  -> tight constraints
  // if(!tight) -> loose constraints

  // also look into how to reconstruct pi0 soon

  // can make it stand-alone (removing primary tracks etc)

  VertexingUtils::FCCAnalysesV0 result;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vtx; // FCCAnalyses vertex object
  ROOT::VecOps::RVec<int> pdgAbs;                            // absolute PDG ID
  ROOT::VecOps::RVec<double> invM;                           // invariant mass
  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;

  VertexingUtils::FCCAnalysesVertex V0_vtx;
  
  int nTr = np_tracks.size();
  if(nTr<2) return result;
  ROOT::VecOps::RVec<bool> isInV0(nTr, false);

  edm4hep::Vector3f r_PV = PV.vertex.position; // in mm  
  
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  // push empty tracks to make a size=2 vector
  edm4hep::TrackState tr_i;
  edm4hep::TrackState tr_j;
  tr_pair.push_back(tr_i);
  tr_pair.push_back(tr_j);
  //
  const double m_pi = 0.13957039; // pi+- mass [GeV]
  const double m_p  = 0.93827208; // p+- mass
  const double m_e  = 0.00051099; // e+- mass
  //
  for(unsigned int i=0; i<nTr-1; i++) {
    if(isInV0[i] == true) continue; // don't pair a track if it already forms a V0
    tr_pair[0] = np_tracks[i];

    for(unsigned int j=i+1; j<nTr; j++) {
      if(isInV0[j] == true) continue; // don't pair a track if it already forms a V0
      tr_pair[1] = np_tracks[j];

      V0_vtx = VertexFitterSimple::VertexFitter_Tk(2, tr_pair);

      // constraint on chi2: chi2 < cut (9)
      double chi2_V0 = V0_vtx.vertex.chi2; // normalised but DOF=1
      if(chi2_V0 >= chi2_cut) continue;

      // invariant masses for V0 candidates
      double invM_Ks      = VertexingUtils::get_invM_pairs(V0_vtx, m_pi, m_pi);
      double invM_Lambda1 = VertexingUtils::get_invM_pairs(V0_vtx, m_pi, m_p);
      double invM_Lambda2 = VertexingUtils::get_invM_pairs(V0_vtx, m_p, m_pi);
      double invM_Gamma   = VertexingUtils::get_invM_pairs(V0_vtx, m_e, m_e);

      // V0 candidate distance from PV
      edm4hep::Vector3f r_V0 = V0_vtx.vertex.position; // in mm
      // does Vector3f class has similar functions as root vectors?
      TVector3 r_V0_PV(r_V0[0] - r_PV[0], r_V0[1] - r_PV[1], r_V0[2] - r_PV[2]);
      double r = r_V0_PV.Mag(); // in mm

      // angle b/n V0 candidate momentum & PV-V0 displacement vector
      double p_r = VertexingUtils::get_PV2V0angle(V0_vtx, PV);

      if(tight) {
	// Ks
	if(invM_Ks>0.493 && invM_Ks<0.503 && r>0.5 && p_r>0.999) {
	  if(debug_me) std::cout<<"Found a Ks"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(310);
	  invM.push_back(invM_Ks);
	  break;
	}
      
	// Lambda0
	else if(invM_Lambda1>1.111 && invM_Lambda1<1.121 && r>0.5 && p_r>0.99995) {
	  if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(3122);
	  invM.push_back(invM_Lambda1);
	  break;
	}
	else if(invM_Lambda2>1.111 && invM_Lambda2<1.121 && r>0.5 && p_r>0.99995) {
	  if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(3122);
	  invM.push_back(invM_Lambda2);
	  break;
	}
	
	// photon conversion
	else if(invM_Gamma<0.005 && r>9 && p_r>0.99995) {
	  if(debug_me) std::cout<<"Found a Photon coversion"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(22);
	  invM.push_back(invM_Gamma);
	  break;
	}
      }

      else {
	// Ks
	if(invM_Ks>0.488 && invM_Ks<0.508 && r>0.3 && p_r>0.999) {
	  if(debug_me) std::cout<<"Found a Ks"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(310);
	  invM.push_back(invM_Ks);
	  break;
	}
      
	// Lambda0
	else if(invM_Lambda1>1.106 && invM_Lambda1<1.126 && r>0.3 && p_r>0.999) {
	  if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(3122);
	  invM.push_back(invM_Lambda1);
	  break;
	}
	else if(invM_Lambda2>1.106 && invM_Lambda2<1.126 && r>0.3 && p_r>0.999) {
	  if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(3122);
	  invM.push_back(invM_Lambda2);
	  break;
	}
	
	// photon conversion
	else if(invM_Gamma<0.01 && r>9 && p_r>0.999) {
	  if(debug_me) std::cout<<"Found a Photon coversion"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(22);
	  invM.push_back(invM_Gamma);
	  break;
	}
      }
      
    }
  }

  //std::cout<<"Found "<<vtx.size()<<" V0s"<<std::endl;
  
  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;
  //
  return result;
}

VertexingUtils::FCCAnalysesV0 VertexFinderLCFIPlus::get_V0s_jet(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
								ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
								ROOT::VecOps::RVec<bool> isInPrimary,
								ROOT::VecOps::RVec<fastjet::PseudoJet> jets,
								std::vector<std::vector<int>> jet_consti,
								VertexingUtils::FCCAnalysesVertex PV,
								bool tight,
								double chi2_cut) {
  // V0 reconstruction after jet clustering
  // if(tight)  -> tight constraints
  // if(!tight) -> loose constraints

  // write a separate fn to get non-primary tracks separated by jet so as not to replicate the calculation

  VertexingUtils::FCCAnalysesV0 result;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vtx; // FCCAnalyses vertex object
  ROOT::VecOps::RVec<int> pdgAbs;                            // absolute PDG ID
  ROOT::VecOps::RVec<double> invM;                           // invariant mass
  ROOT::VecOps::RVec<int> nSV_jet;
  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;
  result.nSV_jet = nSV_jet;

  VertexingUtils::FCCAnalysesVertex V0_vtx;

  int n_par = recoparticles.size();
  if(n_par<2) return result;
  
  edm4hep::Vector3f r_PV = PV.vertex.position; // in mm

  const double m_pi = 0.13957039; // pi+- mass [GeV]
  const double m_p  = 0.93827208; // p+- mass
  const double m_e  = 0.00051099; // e+- mass

  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  // push empty tracks to make a size=2 vector
  edm4hep::TrackState tr_i;
  edm4hep::TrackState tr_j;
  tr_pair.push_back(tr_i);
  tr_pair.push_back(tr_j);

  // find V0s inside the jet loop (only from non-primary tracks)
  // first separate reco particles by jet then get the associated tracks
  int nJet = jets.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks;
  
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks   = ReconstructedParticle2Track::getRP2TRK( recoparticles, thetracks );
  ROOT::VecOps::RVec<int> reco_ind_tracks     = ReconstructedParticle2Track::get_recoindTRK( recoparticles, thetracks );
  if(tracks.size() != reco_ind_tracks.size()) std::cout<<"ERROR: reco index vector not the same size as no of tracks"<<std::endl;

  if(tracks.size() != isInPrimary.size()) std::cout<<"ERROR: isInPrimary vector size not the same as no of tracks"<<std::endl;
  
  if(debug_me) std::cout<<"tracks extracted from the reco particles"<<std::endl;
  
  //
  for (unsigned int j=0; j<nJet; j++) {

    int i_nSV = 0;
    
    std::vector<int> i_jetconsti = jet_consti[j];
    for (int ctr=0; ctr<tracks.size(); ctr++) {
      if(isInPrimary[ctr]) continue; // remove primary tracks
      if(std::find(i_jetconsti.begin(), i_jetconsti.end(), reco_ind_tracks[ctr]) == i_jetconsti.end()) {
	np_tracks.push_back(tracks[ctr]); // separate tracks by jet
      }
    }
    
    if(debug_me) std::cout<<"primary tracks removed; there are "<<np_tracks.size()<<" non-primary tracks in jet#"<<j+1<<std::endl;

    int nTr = np_tracks.size();
    if(nTr<2) continue;    
    ROOT::VecOps::RVec<bool> isInV0(nTr, false);
    
    //
    for(unsigned int i=0; i<nTr-1; i++) {
      if(isInV0[i] == true) continue; // don't pair a track if it already forms a V0
      tr_pair[0] = np_tracks[i];
      
      for(unsigned int j=i+1; j<nTr; j++) {
	if(isInV0[j] == true) continue; // don't pair a track if it already forms a V0
	tr_pair[1] = np_tracks[j];
	
	V0_vtx = VertexFitterSimple::VertexFitter_Tk(2, tr_pair);
	
	// constraint on chi2: chi2 < cut (9)
	double chi2_V0 = V0_vtx.vertex.chi2; // normalised but DOF=1
	if(chi2_V0 >= chi2_cut) continue;
	
	// invariant masses for V0 candidates
	double invM_Ks      = VertexingUtils::get_invM_pairs(V0_vtx, m_pi, m_pi);
	double invM_Lambda1 = VertexingUtils::get_invM_pairs(V0_vtx, m_pi, m_p);
	double invM_Lambda2 = VertexingUtils::get_invM_pairs(V0_vtx, m_p, m_pi);
	double invM_Gamma   = VertexingUtils::get_invM_pairs(V0_vtx, m_e, m_e);
	
	// V0 candidate distance from PV
	edm4hep::Vector3f r_V0 = V0_vtx.vertex.position; // in mm
	// does Vector3f class has similar functions as root vectors?
	TVector3 r_V0_PV(r_V0[0] - r_PV[0], r_V0[1] - r_PV[1], r_V0[2] - r_PV[2]);
	double r = r_V0_PV.Mag(); // in mm
	
	// angle b/n V0 candidate momentum & PV-V0 displacement vector
	double p_r = VertexingUtils::get_PV2V0angle(V0_vtx, PV);
	
	if(tight) {
	  // Ks
	  if(invM_Ks>0.493 && invM_Ks<0.503 && r>0.5 && p_r>0.999) {
	    if(debug_me) std::cout<<"Found a Ks"<<std::endl;
	    isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(310);
	    invM.push_back(invM_Ks);
	    i_nSV++;
	    break;
	  }
	  
	  // Lambda0
	  else if(invM_Lambda1>1.111 && invM_Lambda1<1.121 && r>0.5 && p_r>0.99995) {
	    if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	    isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(3122);
	    invM.push_back(invM_Lambda1);
	    i_nSV++;
	    break;
	  }
	  else if(invM_Lambda2>1.111 && invM_Lambda2<1.121 && r>0.5 && p_r>0.99995) {
	    if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	    isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(3122);
	    invM.push_back(invM_Lambda2);
	    i_nSV++;
	    break;
	  }
	
	  // photon conversion
	  else if(invM_Gamma<0.005 && r>9 && p_r>0.99995) {
	    if(debug_me) std::cout<<"Found a Photon coversion"<<std::endl;
	    isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(22);
	    invM.push_back(invM_Gamma);
	    i_nSV++;
	    break;
	  }
	}

	else {
	  // Ks
	  if(invM_Ks>0.488 && invM_Ks<0.508 && r>0.3 && p_r>0.999) {
	    if(debug_me) std::cout<<"Found a Ks"<<std::endl;
	    isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(310);
	    invM.push_back(invM_Ks);
	    i_nSV++;
	    break;
	  }
      
	  // Lambda0
	  else if(invM_Lambda1>1.106 && invM_Lambda1<1.126 && r>0.3 && p_r>0.999) {
	    if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	    isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(3122);
	    invM.push_back(invM_Lambda1);
	    i_nSV++;
	    break;
	  }
	  else if(invM_Lambda2>1.106 && invM_Lambda2<1.126 && r>0.3 && p_r>0.999) {
	    if(debug_me) std::cout<<"Found a Lambda0"<<std::endl;
	    isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(3122);
	    invM.push_back(invM_Lambda2);
	    i_nSV++;
	    break;
	  }
	
	  // photon conversion
	  else if(invM_Gamma<0.01 && r>9 && p_r>0.999) {
	    if(debug_me) std::cout<<"Found a Photon coversion"<<std::endl;
	    isInV0[i] = true;
	    isInV0[j] = true;
	    vtx.push_back(V0_vtx);
	    pdgAbs.push_back(22);
	    invM.push_back(invM_Gamma);
	    i_nSV++;
	    break;
	  }
	}
      
      }
    }

    nSV_jet.push_back(i_nSV);
    // clean-up
    np_tracks.clear();
  } // jet loop ends

  //std::cout<<"Found "<<vtx.size()<<" V0s"<<std::endl;
  
  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;
  result.nSV_jet = nSV_jet;
  //
  return result;
}


// ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> VertexFinderLCFIPlus::VertexSeed_all(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
// 					  				            VertexingUtils::FCCAnalysesVertex PV,
// 									            double chi2_cut, double invM_cut) {

  // // gives indices of the all pairs of tracks which pass the constraints

  // ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> result;
  // ROOT::VecOps::RVec<int> ij_sel;
  
  // int nTr = tracks.size();
  // ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  // // push empty tracks to make a size=2 vector
  // edm4hep::TrackState tr_i;
  // edm4hep::TrackState tr_j;
  // tr_pair.push_back(tr_i);
  // tr_pair.push_back(tr_j);
  // VertexingUtils::FCCAnalysesVertex vtx_seed;
  
//   for(unsigned int i=0; i<nTr-1; i++) {
//     if(i!=0) tr_pair[0] = tracks[i];

//     for(unsigned int j=i+1; j<nTr; j++) {
//       if(j!=1) tr_pair[1] = tracks[j];

//       // V0 rejection (loose)
//       ROOT::VecOps::RVec<bool> isInV0 = isV0(tr_pair, PV, false);
//       if(isInV0[0] && isInV0[1]) continue;
      
//       vtx_seed = VertexFitterSimple::VertexFitter_Tk(2, tr_pair);

//       // Constraints
//       bool pass = check_constraints(vtx_seed, tr_pair, PV, true, chi2_cut, invM_cut, chi2Tr_cut);
//       if(!pass) continue;

//       // if a pair passes all constraints, store indices
//       ij_sel.push_back(i); ij_sel.push_back(j);
//       result.push_back(ij_sel);
//       ij_sel.clear();
//     }
//   }

//   return result;
// }


// ROOT::VecOps::RVec<int> VertexFinderLCFIPlus::addTrack_multi(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
// 							        ROOT::VecOps::RVec<int> vtx_tr,
// 							        VertexingUtils::FCCAnalysesVertex PV,
// 							        double chi2_cut, double invM_cut, double chi2Tr_cut) {
//   // adds indices of all tracks passing constraints to the (seed) vtx
  
//   ROOT::VecOps::RVec<int> result = vtx_tr;
//   if(tracks.size() == vtx_tr.size()) return result;

//   int nTr = tracks.size();
//   ROOT::VecOps::RVec<edm4hep::TrackState> tr_vtx;
//   VertexingUtils::FCCAnalysesVertex vtx;

//   // tracks from the previously formed vtx
//   for(int tr : vtx_tr) {
//     tr_vtx.push_back(tracks[tr]);
//   }
//   int iTr = tr_vtx.size();
  
//   // find best track to add to the vtx
//   for(unsigned int i=0; i<nTr; i++) {
//     if(std::find(vtx_tr.begin(), vtx_tr.end(), i) != vtx_tr.end()) continue;

//     if(iTr != tr_vtx.size()) tr_vtx[iTr] = tracks[i];
//     else tr_vtx.push_back(tracks[i]);
    
//     vtx = VertexFitterSimple::VertexFitter_Tk(2, tr_vtx);

//     // Constraints
//     bool pass = check_constraints(vtx_seed, tr_pair, PV, false, chi2_cut, invM_cut, chi2Tr_cut);
//     if(!pass) continue;
    
//     // if a pair passes all constraints add to the vtx
//     result.push_back(i);
//     iTr++;
//   }

//   return result;
// }
