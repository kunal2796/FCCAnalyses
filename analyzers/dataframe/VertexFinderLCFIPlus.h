#ifndef  VERTEXFINDERLCFIPLUS_ANALYZERS_H
#define  VERTEXFINDERLCFIPLUS_ANALYZERS_H

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
#include "VertexFitterSimple.h"
#include "VertexingUtils.h"

#include "edm4hep/VertexData.h"
#include "edm4hep/Vertex.h"

#include "fastjet/JetDefinition.hh"

/** Primary and Seconday Vertex Finder interface using vertex fitter from VertexFitterSimple. 
This represents a set functions and utilities to find vertices from a list of tracks following the algorithm from LCFIPlus framework.  
*/

namespace VertexFinderLCFIPlus{

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
  //ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> get_SV_event( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
  VertexingUtils::FCCAnalysesSV get_SV_event( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
					      ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
					      VertexingUtils::FCCAnalysesVertex PV,
					      ROOT::VecOps::RVec<bool> isInPrimary,
					      bool V0_rej=true,
					      double chi2_cut=9., double invM_cut=10., double chi2Tr_cut=5. ) ;

  /** returns SVs reconstructed from non-primary tracks of the event
   *  SV finding done before jet clustering
   */
  //ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> get_SV_event( ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
  VertexingUtils::FCCAnalysesSV get_SV_event( ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
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

  /** adds index of the best track (from the remaining tracks) to the (seed) vtx 
   *  default chi2 threshold is 9 and default invariant mass threshold is 10GeV
   *  default threshold for track's chi2 contribution is 5 (?)
   */
  ROOT::VecOps::RVec<int> addTrack_best( ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
  					 ROOT::VecOps::RVec<int> vtx_tr,
  					 VertexingUtils::FCCAnalysesVertex PV,
  					 double chi2_cut=9., double invM_cut=10., double chi2Tr_cut=5.) ;

  /** V0 rejection (tight)
   *  takes all (non-primary tracks) & removes tracks coming from V0s if user chooses
   *  by default V0 rejection is done
   */
  ROOT::VecOps::RVec<edm4hep::TrackState> V0rejection_tight( ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							     VertexingUtils::FCCAnalysesVertex PV,
							     bool V0_rej=true ) ;

  /** find SVs from a set of tracks
   *  default values of thresholds for the constraints are set
   */
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> findSVfromTracks( ROOT::VecOps::RVec<edm4hep::TrackState> tracks_fin,
									  VertexingUtils::FCCAnalysesVertex PV,
									  double chi2_cut=9., double invM_cut=10., double chi2Tr_cut=5.) ;

  /** check constraints of vertex candidates
   *  default values of thresholds for the constraints are set
   *  default constraint check is that for finding vertex seed
   *  seed=true -> constraints for seed; seed=false -> constraints for adding tracks
   */
  bool check_constraints( VertexingUtils::FCCAnalysesVertex vtx,
			  ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
			  VertexingUtils::FCCAnalysesVertex PV,
			  bool seed=true,
			  double chi2_cut=9., double invM_cut=10., double chi2Tr_cut=5.) ;
  
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
					 double chi2_cut=9. ) ;

  VertexingUtils::FCCAnalysesV0 get_V0s_jet( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
					     ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
					     ROOT::VecOps::RVec<bool> isInPrimary,
					     ROOT::VecOps::RVec<fastjet::PseudoJet> jets,
					     std::vector<std::vector<int>> jet_consti,
					     VertexingUtils::FCCAnalysesVertex PV,
					     bool tight = true,
					     double chi2_cut=9. );

  /** returns indices of the all pairs of tracks that pass a set of constraints from a vector of (non-primary) tracks
   *  default chi2 threshold is 9 and default invariant mass threshold is 10GeV
   */
  //ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> VertexSeed_all( ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
  //							        VertexingUtils::FCCAnalysesVertex PV,
  //							        double chi2_cut=9., double invM_cut=10.) ;

  /** adds indices of tracks (from the remaining tracks) that pass a set of constraints to the (seed) vtx 
   *  default chi2 threshold is 9 and default invariant mass threshold is 10GeV
   *  default threshold for track's chi2 contribution is 5 (?)
   */
  //ROOT::VecOps::RVec<int> addTrack_multi( ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
  //					    ROOT::VecOps::RVec<int> vtx_tr,
  //					    VertexingUtils::FCCAnalysesVertex PV,
  //					    double chi2_cut=9., double invM_cut=10., double chi2Tr_cut=5.) ;


}

#endif
