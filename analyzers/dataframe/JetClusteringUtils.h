#ifndef  JETCLUSTERINGUTILS_ANALYZERS_H
#define  JETCLUSTERINGUTILS_ANALYZERS_H

#include "TLorentzVector.h"
#include <cmath>
#include <vector>
#include "Math/Vector4D.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleData.h"
#include "fastjet/JetDefinition.hh"
#include "TRandom3.h"

/** Jet clustering utilities interface. 
This represents a set functions and utilities to perfom jet clustering from a list of.  
*/

namespace JetClusteringUtils{

  /** @name JetClusteringUtils
   *  Jet clustering interface utilities. 
  */
  ///@{
  

  /** Structure to keep useful informations for the jets*/
  struct FCCAnalysesJet{
    ROOT::VecOps::RVec<fastjet::PseudoJet> jets;
    std::vector<std::vector<int>> constituents;
  };

  /** Set fastjet pseudoJet for later reconstruction*/
  std::vector<fastjet::PseudoJet> set_pseudoJets(ROOT::VecOps::RVec<float> px, 
						 ROOT::VecOps::RVec<float> py, 
						 ROOT::VecOps::RVec<float> pz, 
						 ROOT::VecOps::RVec<float> e);

  /** Add more particles to fastjet pseudoJet for later reconstruction
   *(ghostmatching trial)
  */
  std::vector<fastjet::PseudoJet> addMore_pseudoJets(std::vector<fastjet::PseudoJet> pseudoJ,
						     ROOT::VecOps::RVec<float> px, 
						     ROOT::VecOps::RVec<float> py, 
						     ROOT::VecOps::RVec<float> pz, 
						     ROOT::VecOps::RVec<float> e);

  /** Add partons with scaled momenta to fastjet pseudoJet for later reconstruction
   *(ghostmatching trial)
  */
  std::vector<fastjet::PseudoJet> addGhosts_pseudoJets_old(std::vector<fastjet::PseudoJet> pseudoJ,
							   ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin);

  /** Add partons & hadrons with scaled momenta to fastjet pseudoJet for later reconstruction
   *(ghostmatching trial)
  */
  std::vector<fastjet::PseudoJet> addGhosts_pseudoJets(std::vector<fastjet::PseudoJet> pseudoJ,
						       ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin,
						       int statCode);
  /** statCode: 0->23, 1->71-79 (choice of parton status codes) */
  
  /** Add partons with scaled momenta to fastjet pseudoJet for later reconstruction
   *(ghostmatching trial with status 71-79)
  */
  std::vector<fastjet::PseudoJet> addGhosts7x_pseudoJets(std::vector<fastjet::PseudoJet> pseudoJ,
							 ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin);

  /** Get PDG IDs of reco+ghosts
   * PDG_reco = 0
   * PDG_ghosts = PDG
  */
  ROOT::VecOps::RVec<float> get_gmPDG(ROOT::VecOps::RVec<float> pdg, 
				      ROOT::VecOps::RVec<float> px, 
				      ROOT::VecOps::RVec<float> px_g);
  
  /** Set fastjet pseudoJet for later reconstruction using px, py, pz and m
   * 
   * This version is to be preferred over the px,py,pz,E version when m is known 
   * accurately, because it uses double precision to reconstruct the energy, 
   * reducing the size of rounding errors on FastJet calculations (e.g. of
   * PseudoJet masses)
   * 
  */
  std::vector<fastjet::PseudoJet> set_pseudoJets_xyzm(ROOT::VecOps::RVec<float> px, 
						 ROOT::VecOps::RVec<float> py, 
						 ROOT::VecOps::RVec<float> pz, 
						 ROOT::VecOps::RVec<float> m);

  /** Get fastjet pseudoJet after reconstruction from FCCAnalyses jets*/
  ROOT::VecOps::RVec<fastjet::PseudoJet> get_pseudoJets(FCCAnalysesJet);

  /** Get fastjet constituents after reconstruction from FCCAnalyses jets*/
  std::vector<std::vector<int>> get_constituents(FCCAnalysesJet);
  
  /** Get jet px. Details. */
  ROOT::VecOps::RVec<float> get_px(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
  
  /** Get jet py. Details. */
  ROOT::VecOps::RVec<float> get_py(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
  
  /** Get jet pz. Details. */
  ROOT::VecOps::RVec<float> get_pz(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
  
  /** Get jet energy. Details. */
  ROOT::VecOps::RVec<float> get_e(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
  
  /** Get jet pt. Details. */
  ROOT::VecOps::RVec<float> get_pt(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
  
  /** Get jet p. Details. */
  ROOT::VecOps::RVec<float> get_p(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
  
  /** Get jet mass. Details. */
  ROOT::VecOps::RVec<float> get_m(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
  
  /** Get jet eta. Details. */
  ROOT::VecOps::RVec<float> get_eta(ROOT::VecOps::RVec<fastjet::PseudoJet> in);

  /** Get jet phi. Details. */
  ROOT::VecOps::RVec<float> get_phi(ROOT::VecOps::RVec<fastjet::PseudoJet> in);
  
  /** Get jet theta. Details. */
  ROOT::VecOps::RVec<float> get_theta(ROOT::VecOps::RVec<fastjet::PseudoJet> in);


  ///Internal methods
  FCCAnalysesJet initialise_FCCAnalysesJet();

  FCCAnalysesJet build_FCCAnalysesJet(std::vector<fastjet::PseudoJet> in);
  
  std::vector<fastjet::PseudoJet> build_jets(fastjet::ClusterSequence & cs, int exclusive, float cut, int sorted);

  bool check(unsigned int n, int exclusive, float cut);

  fastjet::RecombinationScheme recomb_scheme(int recombination);
  
  ///@}
}


#endif
