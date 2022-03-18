// LCFI+ (SV Finding - Clustering First)
// 1. select all non-primary tracks in a jet
// 2. remove V0s (or store them separately)
// 3. find best possible pair of tracks (seed vertex) using vertex fitter
// 4. keep adding the best possible track to the seed (non V0) until no track passes the constraints
// 5. remove tracks forming this vertex from the set of non-primary tracks
// 6. repeat 2, 3, 4 until no seed more seed vertices can be found

// LCFI+ (Vertex Merging - Clustering First)
// 1. merge best possible pair of vertices as long as resulting vertex passes the constraints
// 2. repeat until there are at-most 2 vertices in the jet

// LCFI+ (Pseudo-vertices)
// 1. find pseudo-vertices from the tracks that don't form any vertex

// one function to output a vector of SVs
// another to output a vector of psedo-vertices
// yet another to get V0s (isV0 tells which tracks form a V0)

// change to vec of vec to separate SV from diff jets
ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> VertexFitterSimple::get_SV(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
										 ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
										 VertexingUtils::FCCAnalysesVertex PV,
										 ROOT::VecOps::RVec<bool> isInPrimary,
										 std::vector<fastjet::PseudoJet> jets,
										 std::vector<std::vector<int>> jet_consti) {

  // find SVs using LCFI+ (clustering first)
  // still need to think a little about jet clustering using SVs & pseudo-vertices as seeds
  // write one for finding SVs in jets, another for finding SVs in whole evt
  
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;
  
  // retrieve the tracks associated to the recoparticles
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks = ReconstructedParticle2Track::getRP2TRK( recoparticles, thetracks );

  if(tracks.size() != isInPrimary) cout<<"ISSUE: track vector and primary-nonprimary vector of diff sizes"<<endl;

  /*
  // separate the tracks from different jets (only non-primary)
  int nJet = jets.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks_j[nJet]; // access track#3 in jet#2 from this collection as tracks_j[1].at(2)
  //
  for (unsigned int j=0; j<nJet; j++) {
    for (unsigned int ele : jet_consti.at(j)) {
      if (!isInPrimary.at(ele)) tracks_j[j].push_back(tracks.at(ele));
    }
  }
  */

  // V0 rejection
  // should all this be inside the jet loop
  // that way don't need to define the track collection as an array

  // find SV inside the jet loop (only from non-primary tracks)
  int nJet = jets.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks_j;
  //
  for (unsigned int j=0; j<nJet; j++) {
    for (unsigned int ele : jet_consti.at(j)) {
      if (!isInPrimary.at(ele)) tracks_j.push_back(tracks.at(ele));
    }

    // V0 rejection (tight)
    ROOT::VecOps::RVec<edm4hep::TrackState> tracks_fin;
    bool tight = true;
    ROOT::VecOps::RVec<bool> isInV0 = isV0(tracks_j, PV, tight);
    for(unsigned int i=0; i<isInV0.size(); i++) {
      if (!isInV0[i]) tracks_fin.push_back(tracks_j[i]);
    }

    while(tracks_fin.size > 1) {
      // find vertex seed
      ROOT::VecOps::RVec<int> vtx_seed = VertexSeed_best(tracks_fin, PV);
      // constraints can be chosen by user, here using default cuts
      if(vtx_seed.size() == 0) break;
      
      // add tracks to the seed
      // check if a track is added; if not break loop
      ROOT::VecOps::RVec<int> vtx_fin = vtx_seed;
      int vtx_fin_size = 0; // to start the loop
      while(vtx_fin_size != vtx_fin.size()) {
	vtx_fin_size = vtx_fin.size();
	vtx_fin = addTrack_best(tracks_fin, vtx_fin, PV);
      // constraints can be chosen by user, here using default cuts
      }
      
      // fit tracks to SV and remove from tracks_fin
      VertexingUtils::FCCAnalysesVertex sec_vtx = VertexFitter_Tk(0, vtx_fin);
      result.push_back(sec_vtx);
      //
      ROOT::VecOps::RVec<edm4hep::TrackState> temp = tracks_fin;
      tracks_fin.clear();
      for(unsigned int t=0; t<temp.size(); t++) {
	if(std::find(vtx_fin.begin(), vtx_fin.end(), t) == vtx_fin.end()) tracks_fin.push_back(temp[t]);
      }
      // all this cause don't know how to remove multiple elements at once
    }

    tracks_j.clear();
  }

  // currently don't know which SV is from which jet (FIX SOON)
  return result;
}

ROOT::VecOps::RVec<int> VertexFitterSimple::VertexSeed_best(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							    VertexingUtils::FCCAnalysesVertex PV,
							    double chi2_cut, double invM_cut) {

  // gives indices of the best pair of tracks
  // maybe also update and write one to get first pair to pass constraints

  ROOT::VecOps::RVec<int> result;
  int isel, jsel;
  
  int nTr = tracks.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  tr_pair.push_back(tracks[0]);
  tr_pair.push_back(tracks[1]);
  VertexingUtils::FCCAnalysesVertex vtx_seed;
  double chi2_min = 99;
  
  for(unsigned int i=0; i<nTr-1; i++) {
    if(i!=0) tr_pair[0] = tracks[i];

    for(unsigned int j=i+1; j<nTr; j++) {
      if(j!=1) tr_pair[1] = tracks[j];

      // V0 rejection (loose)
      ROOT::VecOps::RVec<bool> isInV0 = isV0(tr_pair, PV, false);
      if(isInV0[0] && isInV0[1]) continue;
      
      vtx_seed = VertexFitter_Tk(0, tr_pair);

      // Constraints
      // chi2 < cut (9)
      double chi2_seed = vtx_seed.vertex.chi2; // normalised
      if(chi2_seed >= chi2_cut) continue; // nDOF for 2 track vtx = 1
      //
      // invM < cut (10GeV)
      double invM_seed = get_invM(vtx_seed);
      if(invM_seed >= invM_cut) continue;
      //
      // invM < sum of energy
      double E_pair = 0.;
      for(edm4hep::TrackState tr_e : tr_pair) E_pair += get_trackE(tr_e);
      if(invM_seed >= E_pair) continue;
      //
      // momenta sum & vtx r on same side
      double angle = get_PV2vtx_angle(tr_pair, vtx_seed, PV);
      if(angle<0) continue;

      // if a pair passes all constraints compare chi2, store lowest chi2
      if(chi2_seed < chi2_min) {
	isel = i; jsel =j;
	chi2_min = chi2_seed;
      }
    }
  }

  result.push_back(isel); result.push_back(jsel);
  return result;
}

std::vector<std::vector<int>> VertexFitterSimple::VertexSeed_all(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
								 VertexingUtils::FCCAnalysesVertex PV,
								 double chi2_cut, double invM_cut) {

  // gives indices of the all pairs of tracks which pass the constraints

  std::vector<std::vector<int>> result;
  std::vector<int> ij_sel;
  
  int nTr = tracks.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  tr_pair.push_back(tracks[0]);
  tr_pair.push_back(tracks[1]);
  VertexingUtils::FCCAnalysesVertex vtx_seed;
  
  for(unsigned int i=0; i<nTr-1; i++) {
    if(i!=0) tr_pair[0] = tracks[i];

    for(unsigned int j=i+1; j<nTr; j++) {
      if(j!=1) tr_pair[1] = tracks[j];

      // V0 rejection (loose)
      ROOT::VecOps::RVec<bool> isInV0 = isV0(tr_pair, PV, false);
      if(isInV0[0] && isInV0[1]) continue;
      
      vtx_seed = VertexFitter_Tk(0, tr_pair);

      // Constraints
      // chi2 < cut (9)
      double chi2_seed = vtx_seed.vertex.chi2; // normalised
      if(chi2_seed >= chi2_cut) continue; // nDOF for 2 track vtx = 1
      //
      // invM < cut (10GeV)
      double invM_seed = get_invM(vtx_seed);
      if(invM_seed >= invM_cut) continue;
      //
      // invM < sum of energy
      double E_pair = 0.;
      for(edm4hep::TrackState tr_e : tr_pair) E_pair += get_trackE(tr_e);
      if(invM_seed >= E_pair) continue;
      //
      // momenta sum & vtx r on same side
      double angle = get_PV2vtx_angle(tr_pair, vtx_seed, PV);
      if(angle<0) continue;

      // if a pair passes all constraints, store indices
      ij_sel.push_back(i); ij_sel.push_back(j);
      result.push_back(ij_sel);
      ij_sel.clear();
    }
  }

  return result;
}

ROOT::VecOps::RVec<int> VertexFitterSimple::addTrack_best(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
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
    tr_vtx.push_back(tracks[tr]);
  }
  int iTr = tr_vtx.size();
  tr_vtx.push_back(tracks[0]);

  // find best track to add to the vtx
  for(unsigned int i=0; i<nTr; i++) {
    if(std::find(vtx_tr.begin(), vtx_tr.end(), i) != vtx_tr.end()) continue;
    if(i!=0) tr_vtx[iTr] = tracks[i];
    
    vtx = VertexFitter_Tk(0, tr_vtx);

    // Constraints
    // chi2_contribution(track) < threshold
    ROOT::VecOps::RVec<float> chi2_tr = vtx.reco_chi2;
    if(chi2_tr[iTr] >= chi2Tr_cut) continue; // threshold = 5 ok?
    //
    // chi2 < cut (9)
    double chi2_vtx = vtx.vertex.chi2; // normalised
    double nDOF = 2*(iTr+1) - 3; // nDOF = 2*nTr - 3
    chi2_vtx = chi2_vtx * nDOF;
    if(chi2_vtx >= chi2_cut) continue;
    //
    // invM < cut (10GeV)
    double invM_vtx = get_invM(vtx);
    if(invM_vtx >= invM_cut) continue;
    //
    // invM < sum of energy (should it be or not?)
    //
    // momenta sum & vtx r on same side
    double angle = get_PV2vtx_angle(tr_vtx, vtx, PV);
    if(angle<0) continue;
    
    // if a track passes all constraints compare chi2, store lowest chi2
    if(chi2_vtx < chi2_min) {
      isel = i;
      chi2_min = chi2_vtx;
    }    
  }

  if(isel>=0) result.push_back(isel);
  return result;
}

ROOT::VecOps::RVec<int> VertexFitterSimple::addTrack_multi(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							   ROOT::VecOps::RVec<int> vtx_tr,
							   VertexingUtils::FCCAnalysesVertex PV,
							   double chi2_cut, double invM_cut, double chi2Tr_cut) {
  // adds indices of all tracks passing constraints to the (seed) vtx
  
  ROOT::VecOps::RVec<int> result = vtx_tr;
  if(tracks.size() == vtx_tr.size()) return result;

  int nTr = tracks.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_vtx;
  VertexingUtils::FCCAnalysesVertex vtx;

  // tracks from the previously formed vtx
  for(int tr : vtx_tr) {
    tr_vtx.push_back(tracks[tr]);
  }
  int iTr = tr_vtx.size();
  
  // find best track to add to the vtx
  for(unsigned int i=0; i<nTr; i++) {
    if(std::find(vtx_tr.begin(), vtx_tr.end(), i) != vtx_tr.end()) continue;

    if(iTr != tr_vtx.size()) tr_vtx[iTr] = tracks[i];
    else tr_vtx.push_back(tracks[i]);
    
    vtx = VertexFitter_Tk(0, tr_vtx);

    // Constraints
    // chi2_contribution < threshold
    ROOT::VecOps::RVec<float> chi2_tr = vtx.reco_chi2;
    if(chi2_tr[iTr] >= chi2Tr_cut) continue; // threshold = 5 ok?
    //
    // chi2 < cut (9)
    double chi2_vtx = vtx.vertex.chi2; // normalised
    double nDOF = 2*(iTr+1) - 3; // nDOF = 2*nTr - 3
    chi2_vtx = chi2_vtx * nDOF;
    if(chi2_vtx >= chi2_cut) continue;
    //
    // invM < cut (10GeV)
    double invM_vtx = get_invM(vtx);
    if(invM_vtx >= invM_cut) continue;
    //
    // invM < sum of energy (should it be or not?)
    //
    // momenta sum & vtx r on same side
    double angle = get_PV2vtx_angle(tr_vtx, vtx, PV);
    if(angle<0) continue;
    
    // if a pair passes all constraints add to the vtx
    result.push_back(i);
    iTr++;
  }

  return result;
}

ROOT::VecOps::RVec<bool> VertexFitterSimple::isV0(ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
						  VertexingUtils::FCCAnalysesVertex PV,
						  bool tight) {
  // V0 rejection
  // fn can be updated to reconstruct V0 and output its momentum, PID, etc instead

  // take all non-primary tracks & assign "true" to pairs that form V0
  // if(tight)  -> tight constraints
  // if(!tight) -> loose constraints

  int nTr = np_tracks.size();

  ROOT::VecOps::RVec<bool> result(nTr, false);
  // true -> forms a V0, false -> doesn't form a V0
  if(nTracks<2) return result;
  
  edm4hep::Vector3f r_PV = PV.vertex.position; // in mm  
  
  ROOT::VecOps::RVec<edm4hep::TrackState> t_pair;
  t_pair.push_back(np_tracks[0]);
  t_pair.push_back(np_tracks[1]);
  VertexingUtils::FCCAnalysesVertex V0;
  //
  const double m_pi = 0.13957039; // pi+- mass [GeV]
  const double m_p  = 0.93827208; // p+- mass
  const double m_e  = 0.00051099; // e+- mass
  //
  for(unsigned int i=0; i<nTr-1; i++) {
    if(result[i] == true) continue;
    if(i!=0) t_pair[0] = np_tracks[i];

    for(unsigned int j=i+1; j<nTr; j++) {
      if(result[j] == true) continue;
      if(j!=1) t_pair[1] = np_tracks[j];

      V0 = VertexFitter_Tk(0, t_pair);

      // invariant masses for V0 candidates
      double invM_Ks      = get_invM_pairs(V0, m_pi, m_pi);
      double invM_Lambda1 = get_invM_pairs(V0, m_pi, m_p);
      double invM_Lambda2 = get_invM_pairs(V0, m_p, m_pi);
      double invM_Gamma   = get_invM_pairs(V0, m_e, m_e);

      // V0 candidate distance from PV
      edm4hep::Vector3f r_V0 = V0.vertex.position; // in mm
      // does Vector3f class has similar functions as root vectors?
      TVector3 r_V0_PV = (r_V0[0] - r_PV[0], r_V0[1] - r_PV[1], r_V0[2] - r_PV[2]);
      double r = r_V0_PV.Mag(); // in mm

      // angle b/n V0 candidate momentum & PV-V0 displacement vector
      double p_r = get_PV2V0angle(V0, PV);

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

// invariant mass of a two track vertex
double VertexFitterSimple::get_invM_pairs(VertexingUtils::FCCAnalysesVertex vertex,
					  double m1,
					  double m2) {
  // CAUTION: m1 -> first track; m2 -> second track

  double result;
  
  ROOT::VecOps::RVec<TVector3> p_tracks = vertex.updated_track_momentum_at_vertex;

  TLorentzVector p4_vtx;
  double m[2] = {m1, m2};
  int nTr = p_tracks.size();
  
  for(unsigned int i=0; i<nTr; i++) {
    TLorentzVector p4_tr;
    p4_tr.SetXYZM(p_tracks[i].X(), p_tracks[i].Y(), p_tracks[i].Z(), m[i]);
    p4_vtx += p4_tr;
  }

  result = p4_vtx.M();
  return result;
}

// invariant mass of a vertex (assuming all tracks to be pions)
double VertexFitterSimple::get_invM(VertexingUtils::FCCAnalysesVertex vertex) {

  double result;
  
  ROOT::VecOps::RVec<TVector3> p_tracks = vertex.updated_track_momentum_at_vertex;

  TLorentzVector p4_vtx;
  const double m = 0.13957039; // pion mass

  for(TVector3 p_tr : p_tracks) {
    TLorentzVector p4_tr;
    p4_tr.SetXYZM(p_tr.X(), p_tr.Y(), p_tr.Z(), m);
    p4_vtx += p4_tr;
  }

  result = p4_vtx.M();
  return result;
}

// cos(angle) b/n V0 candidate's momentum & PV to V0 displacement vector
double VertexFitterSimple::get_PV2V0angle(VertexingUtils::FCCAnalysesVertex V0,
					  VertexingUtils::FCCAnalysesVertex PV) {
  double result;

  ROOT::VecOps::RVec<TVector3> p_tracks = V0.updated_track_momentum_at_vertex;

  TVector3 p_vtx;
  for(TVector3 p_tr : p_tracks) p_vtx += p_tr;

  edm4hep::Vector3f r_V0 = V0.vertex.position; // in mm
  edm4hep::Vector3f r_PV = PV.vertex.position; // in mm

  TVector3 r_V0_PV = (r_V0[0] - r_PV[0], r_V0[1] - r_PV[1], r_V0[2] - r_PV[2]);
  
  double pDOTr = p_vtx.Dot(r_V0_PV);
  double p_mag = p_vtx.Mag();
  double r_mag = r_V0_PV.Mag();

  result = pDOTr / (p_mag * r_mag);
  return result;
}

// cos(angle) b/n track momentum sum & PV to vtx displacement vector
double VertexFitterSimple::get_PV2vtx_angle(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
					    VertexingUtils::FCCAnalysesVertex vtx,
					    VertexingUtils::FCCAnalysesVertex PV) {
  double result;

  TVector3 p_sum;
  for(edm4hep::TrackState tr : tracks) {
    TVectorD ipar = VertexingUtils::get_trackParam(tr);
    TVector3 ip   = ParToP(ipar);
    p_sum += ip;
  }
  
  edm4hep::Vector3f r_vtx = vtx.vertex.position; // in mm
  edm4hep::Vector3f r_PV  = PV.vertex.position;  // in mm

  TVector3 r_vtx_PV = (r_vtx[0] - r_PV[0], r_vtx[1] - r_PV[1], r_vtx[2] - r_PV[2]);
  
  double pDOTr = p_sum.Dot(r_vtx_PV);
  double p_mag = p_sum.Mag();
  double r_mag = r_vtx_PV.Mag();

  result = pDOTr / (p_mag * r_mag);
  return result;
}

double VertexFitterSimple::get_trackE(edm4hep::TrackState track) {

  // get track's energy assuming it to be a pion

  double result;

  const double m_pi = 0.13957039;
  
  TVectorD par = VertexingUtils::get_trackParam(track);
  TVector3 p   = ParToP(par);
  TLorentzVector p4;
  p4.SetXYZM(p[0], p[1], p[2], m_pi);

  result = p4.E();
  return result;
}
