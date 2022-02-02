#ifndef  JETTAGGINGUTILS_ANALYZERS_H
#define  JETTAGGINGUTILS_ANALYZERS_H

#include <vector>
#include "Math/Vector4D.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/MCParticleData.h"
#include "fastjet/JetDefinition.hh"
#include "TRandom3.h"
#include "TLorentzVector.h"

/** Jet tagging utilities interface.
This represents a set functions and utilities to perfom jet tagging from a list of jets.
*/

namespace JetTaggingUtils{

  /** @name JetTaggingUtils
   *  Jet tagging interface utilities.
  */

  //Get flavour association of jet 
  ROOT::VecOps::RVec<int> get_flavour(ROOT::VecOps::RVec<fastjet::PseudoJet> in, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin);
  //Get flavour association of jet 
  ROOT::VecOps::RVec<int> get_flavour_qqbar(ROOT::VecOps::RVec<fastjet::PseudoJet> in, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin);
  //Get flavour association of jet via ghost matching (manual)
  ROOT::VecOps::RVec<int> get_flavour_gm(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<float> pdg_gm);
  //Get flavour association of jet via ghost matching (automated)
  ROOT::VecOps::RVec<int> get_flavour_gm_auto(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin, std::vector<fastjet::PseudoJet> PJin);
  //Get flavour association of jet via ghost matching (automated status 71-79)
  ROOT::VecOps::RVec<int> get_flavour_gm7x_auto(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin, std::vector<fastjet::PseudoJet> PJin);
  //Get flavour association of jet via ghost matching (automated status 71-79, with cut on |p|)
  ROOT::VecOps::RVec<int> get_flavour_gm_pcut(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin, std::vector<fastjet::PseudoJet> PJin, float p_cut);
  //Get b-tags with an efficiency applied
  ROOT::VecOps::RVec<int> get_btag(ROOT::VecOps::RVec<int> in, float efficiency);
  //Get c-tags with an efficiency applied
  ROOT::VecOps::RVec<int> get_ctag(ROOT::VecOps::RVec<int> in, float efficiency);

    ///@}                                                                                                                                                                              
}


#endif
