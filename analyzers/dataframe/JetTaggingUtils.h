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
  //Get b-tags with an efficiency applied
  ROOT::VecOps::RVec<int> get_btag(ROOT::VecOps::RVec<int> in, float efficiency, float mistag_c=0., float mistag_l=0., float mistag_g=0.);
  //Get c-tags with an efficiency applied
  ROOT::VecOps::RVec<int> get_ctag(ROOT::VecOps::RVec<int> in, float efficiency, float mistag_b=0., float mistag_l=0., float mistag_g=0.);
  //Get l-tags with an efficiency applied
  ROOT::VecOps::RVec<int> get_ltag(ROOT::VecOps::RVec<int> in, float efficiency, float mistag_b=0., float mistag_c=0., float mistag_g=0.);
  //Get g-tags with an efficiency applied
  ROOT::VecOps::RVec<int> get_gtag(ROOT::VecOps::RVec<int> in, float efficiency, float mistag_b=0., float mistag_c=0., float mistag_l=0.);

  /// select a list of jets depending on the status of a certain boolean flag (corresponding to its tagging state)
  struct sel_tag {
    bool m_pass; // if pass is true, select tagged jets. Otherwise select anti-tagged ones
    sel_tag(bool arg_pass);
    ROOT::VecOps::RVec<fastjet::PseudoJet> operator() (ROOT::VecOps::RVec<bool> tags, ROOT::VecOps::RVec<fastjet::PseudoJet> in);
  };

  //Get flavour association of jet via ghost matching (manual)
  ROOT::VecOps::RVec<int> get_flavour_gm_manual(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<float> pdg_gm);
  //Get flavour association of jet via ghost matching (automated)
  ROOT::VecOps::RVec<int> get_flavour_gm_auto(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin, std::vector<fastjet::PseudoJet> PJin);
  //Get flavour association of jet via ghost matching (automated status 71-79)
  ROOT::VecOps::RVec<int> get_flavour_gm7x_auto(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin, std::vector<fastjet::PseudoJet> PJin);
  //Get flavour association of jet via ghost matching (automated status 71-79, with cut on |p|)
  ROOT::VecOps::RVec<int> get_flavour_gm_pcut(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin, std::vector<fastjet::PseudoJet> PJin, float p_cut);
  
  //Get flavour association of jet via ghost matching with cut on |p|
  std::vector<std::vector<int>> get_flavour_gm(ROOT::VecOps::RVec<fastjet::PseudoJet> in, std::vector<std::vector<int>> inJC, ROOT::VecOps::RVec<edm4hep::MCParticleData> MCin, std::vector<fastjet::PseudoJet> PJin, int statCode, float p_cut);
  // in: jets that need to be assigned a flavour
  // inJC: jet constituents after clustering (including ghosts)
  // MCin: MC particle collection
  // PJin: jet constituents before clustering (excluding ghosts)
  // statCode: 0->23, 1->71-79 (choice of parton status codes)
  // p_cut: choose 0 for no cut


  ///@}
}

#endif
