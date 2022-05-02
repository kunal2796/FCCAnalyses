#ifndef  VERTEXFITTERSIMPLE_ANALYZERS_H
#define  VERTEXFITTERSIMPLE_ANALYZERS_H

#include <cmath>
#include <vector>

#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/TrackState.h"

#include "TVectorD.h"
#include "TVector3.h"
#include "TMatrixDSym.h"
#include "TMath.h"
#include "TDecompChol.h"
#include "TMatrixD.h"

#include "ReconstructedParticle2Track.h"
#include "ReconstructedParticle2MC.h"
#include "VertexingUtils.h"

#include "edm4hep/VertexData.h"
#include "edm4hep/Vertex.h"

#include "fastjet/JetDefinition.hh"

/** Vertex interface using Franco Bedeshi's code. 
This represents a set functions and utilities to perfom vertexing from a list of tracks.  
*/

namespace VertexFitterSimple{

  /// Vertex (code from Franco Bedeschi): passing the recoparticles. Units for the beamspot constraint: mum
  VertexingUtils::FCCAnalysesVertex  VertexFitter( int Primary, 
						   ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
						   ROOT::VecOps::RVec<edm4hep::TrackState> alltracks,
						   bool BeamSpotConstraint = false,
						   double sigmax=0., double sigmay=0., double sigmaz=0.,
                                                   double bsc_x=0., double bsc_y=0., double bsc_z=0. )  ;


  /// Vertex (code from Franco Bedeschi): passing the tracks. Units for the beamspot constraint: mum
  VertexingUtils::FCCAnalysesVertex  VertexFitter_Tk( int Primary, ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
						      bool BeamSpotConstraint = false,		
						      double sigmax=0., double sigmay=0., double sigmaz=0., 
                                                      double bsc_x=0., double bsc_y=0., double bsc_z=0. )  ;

/// Return the tracks that are flagged as coming from the primary vertex
  ROOT::VecOps::RVec<edm4hep::TrackState> get_PrimaryTracks( VertexingUtils::FCCAnalysesVertex  initialVertex,
							     ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							     bool BeamSpotConstraint,
							     double bsc_sigmax, double bsc_sigmay, double bsc_sigmaz,
							     double bsc_x, double bsc_y, double bsc_z,
							     int ipass = 0 ) ;


/// Return the tracks that are NOT flagged as coming from the primary vertex
  ROOT::VecOps::RVec<edm4hep::TrackState>  get_NonPrimaryTracks( ROOT::VecOps::RVec<edm4hep::TrackState> allTracks,
                                                                 ROOT::VecOps::RVec<edm4hep::TrackState> primaryTracks ) ;

/// for an input vector of tracks, return a  vector of bools that tell if the track  was identified as a primary track
   ROOT::VecOps::RVec<bool> IsPrimary_forTracks( ROOT::VecOps::RVec<edm4hep::TrackState> allTracks,
						 ROOT::VecOps::RVec<edm4hep::TrackState> primaryTracks ) ;


  ///////////////////////////
  //** SV Finder (LCFI+) **//
  ///////////////////////////

  /** returns SVs reconstructed from non-primary tracks of jets
   *  non-primary separated from all tracks using isInPrimary (bool) vector
   *  currently not separating SVs by jet
   */

  VertexingUtils::FCCAnalysesSV get_SV_jets( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
					     ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
					     VertexingUtils::FCCAnalysesVertex PV,
					     ROOT::VecOps::RVec<bool> isInPrimary,
					     ROOT::VecOps::RVec<fastjet::PseudoJet> jets,
					     std::vector<std::vector<int>> jet_consti,
					     bool V0_rej=true,
					     double chi2_cut=9., double invM_cut=10., double chi2Tr_cut=5. ) ;
  
  /** returns SVs reconstructed from non-primary tracks of the event
   *  SV finding done before jet clustering
   *  non-primary separated from all tracks using isInPrimary (bool) vector
   */
  VertexingUtils::FCCAnalysesSV get_SV_event( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
					      ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
					      VertexingUtils::FCCAnalysesVertex PV,
					      ROOT::VecOps::RVec<bool> isInPrimary,
					      bool V0_rej=true,
					      double chi2_cut=9., double invM_cut=10., double chi2Tr_cut=5. ) ;

  /** returns SVs reconstructed from non-primary tracks of the event
   *  SV finding done before jet clustering
   */
  VertexingUtils::FCCAnalysesSV get_SV_event( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
					      ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
					      ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
					      VertexingUtils::FCCAnalysesVertex PV,
					      bool V0_rej=true,
					      double chi2_cut=9., double invM_cut=10., double chi2Tr_cut=5. ) ;

  /** returns a vector of all vertices (PV and SVs), e.g to use in myUtils::get_Vertex_d2PV
   *  first entry: PV, all subsequent entries: SVs
  */
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> get_all_vertices( VertexingUtils::FCCAnalysesVertex PV,
									  VertexingUtils::FCCAnalysesSV SV ); 

  /** returns indices of the best pair of tracks from a vector of (non-primary) tracks 
   *  default chi2 threshold is 9 and default invariant mass threshold is 10GeV
   */
  ROOT::VecOps::RVec<int> VertexSeed_best( ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
					   VertexingUtils::FCCAnalysesVertex PV,
					   double chi2_cut=9., double invM_cut=10.) ;

  /** returns indices of the all pairs of tracks that pass a set of constraints from a vector of (non-primary) tracks
   *  default chi2 threshold is 9 and default invariant mass threshold is 10GeV
   */
  ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> VertexSeed_all( ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							      VertexingUtils::FCCAnalysesVertex PV,
							      double chi2_cut=9., double invM_cut=10.) ;

  /** adds index of the best track (from the remaining tracks) to the (seed) vtx 
   *  default chi2 threshold is 9 and default invariant mass threshold is 10GeV
   *  default threshold for track's chi2 contribution is 5 (?)
   */
  ROOT::VecOps::RVec<int> addTrack_best( ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
					 ROOT::VecOps::RVec<int> vtx_tr,
					 VertexingUtils::FCCAnalysesVertex PV,
					 double chi2_cut=9., double invM_cut=10., double chi2Tr_cut=5.) ;

  /** adds indices of tracks (from the remaining tracks) that pass a set of constraints to the (seed) vtx 
   *  default chi2 threshold is 9 and default invariant mass threshold is 10GeV
   *  default threshold for track's chi2 contribution is 5 (?)
   */
  ROOT::VecOps::RVec<int> addTrack_multi( ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
					  ROOT::VecOps::RVec<int> vtx_tr,
					  VertexingUtils::FCCAnalysesVertex PV,
					  double chi2_cut=9., double invM_cut=10., double chi2Tr_cut=5.) ;

  /** Get the reco indices of all tracks
   */
  ROOT::VecOps::RVec<int> get_reco_ind( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
					ROOT::VecOps::RVec<edm4hep::TrackState> tracks ) ;

  /** V0 rejection/identification
   *  takes all (non-primary) tracks & assigns "true" to pairs that form a V0
   *  if(tight)  -> tight constraints
   *  if(!tight) -> loose constraints
   *  by default loose constraints
   */
  ROOT::VecOps::RVec<bool> isV0( ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
				 VertexingUtils::FCCAnalysesVertex PV,
				 bool tight = false ) ;

  ///
  
  /** returns V0s reconstructed from a set of tracks (as an FCCAnalysesV0 object)
   */
  VertexingUtils::FCCAnalysesV0 get_V0s( ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
					 VertexingUtils::FCCAnalysesVertex PV,
					 bool tight = true,
					 double chi2_cut=9.) ;



  Double_t FastRv(TVectorD p1, TVectorD p2) ;
  TMatrixDSym RegInv3(TMatrixDSym &Smat0) ;
  TMatrixD Fill_A(TVectorD par, Double_t phi) ;
  TVectorD Fill_a(TVectorD par, Double_t phi) ;
  TVectorD Fill_x0(TVectorD par) ;
  TVectorD Fill_x(TVectorD par, Double_t phi) ;

  TVectorD XPtoPar(TVector3 x, TVector3 p, Double_t Q);
  TVector3 ParToP(TVectorD Par);
}

#endif

